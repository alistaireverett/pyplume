#!/usr/bin/python

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import ode
import matplotlib
from scipy.interpolate import interp1d
import gsw
import const
#import ambient

from IPython import embed

class InitAmbient:
    
    def __init__(self,z_max,salinity,temperature,pressure=None,depth=None):

        if depth is None and pressure is None:
            raise ValueError("Must pass either pressure or depth as an argument to InitAmbient")
        elif depth is None and pressure is not None:
            self.pressure = pressure
            self.depth = pressure/(1027.*9.81*1.e-4)
        elif depth is not None and pressure is None:
            self.pressure = depth*(1027.*9.81*1.e-4)
            self.depth = depth
        else:
            if np.any(depth != pressure/(1027.*9.81*1.e-4)):
			    print "UserWarning: You specified a pressure and depth profile which may not match assumptions used elsewhere in this model"
            self.pressure = pressure
            self.depth = depth

        self.z_max = z_max
        self.z = z_max - self.depth 
        self.salinity = salinity
        self.temperature = temperature
        self.rho = gsw.rho(salinity, temperature, self.pressure)
        
        # internal functions to return interpolated values at arbitrary depths
        self.__f_sal_d = interp1d(self.depth, salinity)
        self.__f_temp_d = interp1d(self.depth, temperature)
        self.__f_rho_d = interp1d(self.depth, self.rho)
        self.__f_pres_d = interp1d(self.depth, self.pressure)

    def get_sal_d(self,x):
        return self.__f_sal_d(x)[()]

    def get_temp_d(self,x):
        return self.__f_temp_d(x)[()]

    def get_rho_d(self,x):
        return self.__f_rho_d(x)[()]

    def get_pres_d(self,x):
        return self.__f_pres_d(x)[()]

    def get_sal_z(self,x):
        return self.__f_sal_d(self.z_max - x)[()]

    def get_temp_z(self,x):
        return self.__f_temp_d(self.z_max - x)[()]
                    
    def get_rho_z(self,x):
        return self.__f_rho_d(self.z_max - x)[()]
                    
    def get_pres_z(self,x):
            return self.__f_pres_d(self.z_max - x)[()]
                    

def inlet(h_i,h_w,q):
    
    n_eff = (const.RHO_I*const.G*h_i-const.RHO_W*const.G*h_w)
    assert n_eff > 0., "Terminus is floating - need shallower water or more ice!"
    factor = (const.C1/(const.C2*const.C3**2.*n_eff**const.GLEN_N))**(2./7.)
    
    area = factor*q**(6./7.)
    radius = (2*area/math.pi)**0.5
    velocity = q/area

    return radius, velocity


def get_melt(u_inf,t_inf,s_inf,pressure):

    # solve the three equation melt-parameterisation quadratically
    aa = const.A*(const.GAM_T*const.C_W - const.GAM_S*const.C_I)
    bb = (const.GAM_S*const.C_I*
             (const.A*s_inf - const.B - const.C*pressure + const.T_I 
                - (const.L/const.C_I)) 
             - const.GAM_T*const.C_W*(t_inf - const.B - const.C*pressure))
    cc = const.GAM_S*s_inf*(const.C_I*(const.B + const.C*pressure - const.T_I)                             + const.L)

    s_b = (1./(2.*aa))*(-bb-((bb**2.-4.*aa*cc)**0.5))
    t_b = const.A*s_b+const.B+const.C*pressure
    mdot = const.GAM_S*(const.C_D**0.5)*u_inf*(s_inf-s_b)/s_b

    return t_b, s_b, mdot

def plume_geometry(A,X):

	# A in: plume cross-sectional area
	# X in: distance of center from wall
	# C out: perimeter exposed to ambient
	# W out: perimeter exposed to wall

	B = math.sqrt(A/math.pi)
	
	if (B < X):

		#PLUME IS FULLY DETACHED    

		C = 2.*math.pi*B
		W = 0.

	else:

		#PLUME IS PARTLY OR NON-DETACHED (half-conical)
      	#FIND SHAPE CONSISTENT WITH X AND A

		tol = 1.e-12
		itermax = 300

		for i in range(1,itermax):
			f = B*B * (math.pi-math.acos(X/B)) + X*(B*B-X*X)**0.5 - A
			fp = 2*B*(math.pi - math.acos(X/B)) - X/(1 - X*X/(B*B))**0.5 + (X*B)/(B*B - X*X)**0.5
			ADJ = f/fp
			B = B - ADJ
			if (B < X):
				B = X + .001
			if abs(ADJ) < tol: break
		
		if (i == itermax):
			print "didn't work"

		W = 2 * math.sqrt(B*B-X*X)
		C = 2 * (math.pi-math.acos(X/B)) * B

	return C,W

def Cone(Z,Y,ambient,z_max,MELT=True):

    # The governing set of equations to be solved for the plume
    if Z > z_max:
            return None

    # unpack args
    #[g,z_max,LAMBDA1,LAMBDA2,LAMBDA3,GAMT,GAMS,C_W,C_I,RHO_REF,ICETEMP,L] = args

    # initialise array for output
    YDOT = np.zeros(Y.shape)
    '''
    if MELT:
        t_b, s_b, mdot = get_melt(Y[1],Y[2],Y[3],ambient.get_pres_z(Z))
    else:
        t_b = s_b = mdot = 0.
    '''
    t_amb = ambient.get_temp_z(Z)
    s_amb = ambient.get_sal_z(Z)
    rho_0 = ambient.get_rho_z(Z)
    # approximate pressure at the current depth in Pa 
    pressure = ambient.get_pres_z(Z)

    # calculate current plume density (needs pressure in decibar)
    # gives density in kg/m3
    print pressure
    rho_1 = gsw.rho(Y[3],Y[2],pressure)
    #rho_0 = gsw.rho(Sambient,Tambient,pressure*1.e-4)

    if rho_1 > rho_0:
        Y[0] = 0.
        Y[1] = 0.

    #Cplume = mat
    #print pressure.pi * Y[0]
    #Cplume,Wplume = plume_geometry(Y[0],Y[9])

    YDOT[0] = 2.*const.E_0 - Y[0]*const.G*(rho_0-rho_1)/(2.*Y[1]*Y[1]*const.RHO_REF)

    YDOT[1] = -2.*const.E_0*Y[1]/Y[0] + const.G*(rho_0-rho_1)/(Y[1]*const.RHO_REF)

    YDOT[2] = 2.*const.E_0*(t_amb-Y[2])/Y[0]

    YDOT[3] = 2.*const.E_0*(s_amb-Y[3])/Y[0]

    YDOT[4] = 0.#Wplume

    YDOT[5] = t_amb-Y[5]
    YDOT[6] = s_amb-Y[6]

    YDOT[7] = 0.

    #YDOT[8] = (Y[1]*Y[1]*YDOT[1]+2*Y[1]*Y[0]*YDOT[2])/(Y[1]*Y[1]*Y[0])/math.tan(Y[8])
    #YDOT[9] = 1.0/(math.tan(Y[8]))-Y[9]

    return YDOT

def detachCone(Z,Y,ambient,z_max,MELT=True):

    # The governing set of equations to be solved for the plume
    if Z > z_max:
            return None

    # unpack args
    #[g,z_max,LAMBDA1,LAMBDA2,LAMBDA3,GAMT,GAMS,C_W,C_I,RHO_REF,ICETEMP,L] = args

    # initialise array for output
    YDOT = np.zeros(Y.shape)

    if MELT:
        t_b, s_b, mdot = get_melt(Y[1],Y[2],Y[3],ambient.get_pres_z(Z))
    else:
        t_b = s_b = mdot = 0.

    t_amb = ambient.get_temp_z(Z)
    s_amb = ambient.get_sal_z(Z)
    rho_0 = ambient.get_rho_z(Z)
    # approximate pressure at the current depth in Pa 
    pressure = ambient.get_pres_z(Z)

    # calculate current plume density (needs pressure in decibar)
    # gives density in kg/m3
    rho_1 = gsw.rho(Y[3],Y[2],pressure)
    #rho_0 = gsw.rho(Sambient,Tambient,pressure*1.e-4)

    if rho_1 > rho_0:
        Y[0] = 0.
        Y[1] = 0.

    Cplume = math.pi * Y[0]
    Cplume,Wplume = plume_geometry(Y[0],Y[9])

    YDOT[0] = 2*const.E_0*Cplume/math.sin(Y[8])+2*Wplume*mdot/(math.sin(Y[8])*Y[1])-Y[0]*const.G*(rho_0-rho_1)/(Y[1]*Y[1]*const.RHO_REF)+const.CD*Wplume

    YDOT[1] = -const.E_0*Cplume*Y[1]/Y[0]/math.sin(Y[8])-Wplume*mdot/(Y[0]*math.sin(Y[8]))+const.G*(rho_0-rho_1)/(Y[1]*const.RHO_REF)-const.CD*Wplume*Y[1]/Y[0]

    YDOT[2] = (const.E_0*Cplume*(Tambient-Y[2])/Y[0]+Wplume*mdot*(Tb-Y[2])/(Y[0]*Y[1])-const.GAMT*math.sqrt(const.CD)*Wplume*(Y[2]-Tb)/(Y[0]))/math.sin(Y[8])

    YDOT[3] = (const.E_0*Cplume*(Sambient-Y[3])/Y[0]+Wplume*mdot*(Sb-Y[3])/(Y[0]*Y[1])-const.GAMS*math.sqrt(const.CD)*Wplume*(Y[3]-Sb)/(Y[0]))/math.sin(Y[8])

    YDiOT[4] = Wplume

    YDOT[5] = Tambient-Y[5]
    YDOT[6] = Sambient-Y[6]

    YDOT[7] = mdot-Y[7]

    YDOT[8] = (Y[1]*Y[1]*YDOT[1]+2*Y[1]*Y[0]*YDOT[2])/(Y[1]*Y[1]*Y[0])/math.tan(Y[8])
    YDOT[9] = 1.0/(math.tan(Y[8]))-Y[9]

    return YDOT

def halfCone(Z,Y,ambient,z_max,MELT=True):

	# The governing set of equations to be solved for the plume
    if Z > z_max:
        return None
    # unpack args
    #[g,z_max,LAMBDA1,LAMBDA2,LAMBDA3,GAMT,GAMS,C_W,C_I,RHO_REF,ICETEMP,L] = args

    # initialise array for output
    YDOT = np.zeros(Y.shape)
    #print Z,ambient.get_sal_z(Z)
    if MELT:
        t_b, s_b, mdot = get_melt(Y[1],Y[2],Y[3],ambient.get_pres_z(Z))
    else:
        t_b = s_b = mdot = 0.

    t_amb = ambient.get_temp_z(Z)
    s_amb = ambient.get_sal_z(Z)
    rho_0 = ambient.get_rho_z(Z)
    #Tambient, Sambient, rho_0 = get_TS(Z,ambient.temp,ambient.sal,ambient.depths,ambient.rho) 
    # approximate pressure at the current depth in Pa 
    pressure = ambient.get_pres_z(Z)

    #print "compare pressure: ", 1027.*const.G*(z_max-Z), pressure

    # calculate current plume density (needs pressure in decibar)
    # gives density in kg/m3
    rho_1 = gsw.rho(Y[3],Y[2],pressure)
    #print rho_1,rho_0 
    # check if Neutral Buoyancy is reached
     
    if rho_1 > rho_0:
        Y[0] = 0.
        Y[1] = 0.
        
    
    # Solve the plume equations and store in YDOT
    YDOT[0] = (2.*const.E_0 + 4.*mdot/(math.pi*Y[1])
               - Y[0]*const.G*(rho_0-rho_1)/(2.*Y[1]*Y[1]*const.RHO_REF)
               + 2. * (const.C_D/math.pi))

    YDOT[1] = (- 2.*const.E_0*Y[1]/Y[0] - 4.*mdot/(math.pi*Y[0])
               + const.G*(rho_0-rho_1)/(Y[1]*const.RHO_REF)
               - 4.*const.C_D*Y[1]/(math.pi*Y[0]))

    YDOT[2] = (2.*const.E_0*(t_amb-Y[2])/Y[0]
               + 4.*mdot*(t_b-Y[2])/(math.pi*Y[0]*Y[1])
               - 4.*const.GAM_T*(const.C_D**0.5)*(Y[2]-t_b)/(math.pi*Y[0]))

    YDOT[3] = (2.*const.E_0*(s_amb-Y[3])/Y[0] 
               + 4.*mdot*(s_b-Y[3])/(math.pi*Y[0]*Y[1])
               - 4.*const.GAM_S*(const.C_D**0.5)*(Y[3]-s_b)/(math.pi*Y[0]))

    YDOT[4] = 2.*Y[0]

    YDOT[5] = t_amb - Y[5]
    YDOT[6] = s_amb - Y[6]

    YDOT[7] = mdot - Y[7]

    return YDOT

############################################
'''
# Variables:

T_sg = 1.0e-3    # inlet temperature
S_sg = 1.0e-3    # inlet salinity
w_sg = 1.#.434411779623        # inlet velocity
r_sg = 1.#1.#1.2105682097        # inlet radius
z_max = 1000.

# Detached P
theta_sg = math.pi/2.

#Tambient = 5.       # ambient salinity
#Sambient = 33.      # ambient salinity

detached = False

# Ambient T/S profile

ambient.depths = ambient.pressure/(1027.*const.G*1.e-4)

ambient.rho = np.array([gsw.rho(Samb,Tamb,pressure) for Samb,Tamb,pressure in 
								zip(ambient.sal,ambient.temp,ambient.pressure)])

'''
#get_sal = interp1d(


#------------------------------------------
'''
# Constants:

E_0             = 0.05
ICETEMP  		= -0.
RHO_REF 		= 1020.0
g 				= 9.81
C_W 			= 3994.
C_I				= 2009.
L     			= 334.e3
LAMBDA1 		= -0.0573
LAMBDA2			= 0.0832
LAMBDA3			= 0.000761
GAMT   			= 0.022
GAMS   			= 0.00062
CD     			= 0.0025
backgroundVel 	= 0.001
'''

############################################

def calc_plume(u_sg, r_sg, h_w, ambient, 
               T_sg = 1.0e-3,S_sg = 1.0e-3, 
               detached=False,MELT=True):

    # depth = 0 at surface, increasing down
    # z = z_max at surface, decreasing down (z = z_max - depth)

    #r_sg, u_sg = inlet(h_i,h_w,q)
    #print u_sg, r_sg
    # temporary work around:
    #if z_max > ambient.depth.max(): z_max = ambient.depth.max()
    #print z_max

    #     Y is input/output vector for DLSODE
    #       Y(0) = plume thickness/radius
    #       Y(1) = plume velocity
    #       Y(2) = plume temperature
    #       Y(3) = plume salinity
    #       Y(4) = plume area
    #		Y(5) = ambient temperature
    #		Y(6) = ambient salinity
    #       Y(7) = area integrated melt
    #       Y(8) = angle of plume (detatched plume only)
    #       Y(9) = distance of plume from ice (detatched plume only)

    Y = np.zeros((10,1))

    # Initial conditions
    Y[0] = r_sg			# initial plume thickness
    Y[1] = u_sg         # initial vertical velocity
    Y[2] = T_sg         # initial temperature
    Y[3] = S_sg         # initial salinity
    Y[4] = 0.0          # integrated contact area
    Y[7] = 0.0          # integrated melt rate

    Y[5] = ambient.get_temp_z(0)
    Y[6] = ambient.get_sal_z(0)

    #Y[5], Y[6],_ = get_TS(0,ambient.temp,ambient.sal,ambient.pressure,ambient.rho)

    # detached plume only
    #if detached:
	#   Y[0] = math.pi * r_sg**2
    #    Y[8] = theta_sg
    #    Y[9] = 0.
    #    #Y[9] = delta_y
    #    Y[0] = 0.5 * math.pi * r_sg**2
        

    # set up the depths to iterate over
    z_step = .1
    z_range = np.arange(0,h_w+z_step,z_step)
    
    # Create output arrays to be populated
    rProfPlume = np.zeros(len(z_range))*np.nan
    wProfPlume = np.zeros(len(z_range))*np.nan
    tProfPlume = np.zeros(len(z_range))*np.nan
    sProfPlume = np.zeros(len(z_range))*np.nan
    aProfPlume = np.zeros(len(z_range))*np.nan
    tProfAmb = np.zeros(len(z_range))*np.nan
    sProfAmb = np.zeros(len(z_range))*np.nan
    mIntProfPlume = np.zeros(len(z_range))*np.nan
    thetaProfPlume = np.zeros(len(z_range))*np.nan
    distanceProfPlume = np.zeros(len(z_range))*np.nan

    # Set intial conditions
    rProfPlume[0] = Y[0]
    wProfPlume[0] = Y[1]
    tProfPlume[0] = Y[2]
    sProfPlume[0] = Y[3]
    aProfPlume[0] = Y[4]
    tProfAmb[0] = Y[5]
    sProfAmb[0] = Y[6]
    mIntProfPlume[0] = Y[7]

    #if detached:
    #    thetaProfPlume[0] = Y[8]
    #    distanceProfPlume[0] = Y[9]
    #out = odeint(halfCone,Y,depths)
    
    #embed()



    # Select the solveri
    if detached:
        solver = ode(Cone).set_integrator('lsoda', atol=1.e-5, rtol=1.e-5)
    else:
        solver = ode(halfCone).set_integrator('lsoda', atol=1.e-5, rtol=1.e-5)#,nsteps=1000)#,min_step=0.1)#lsoda
    
    # set initial conidtions in solver with Y array, z=0 and args
    solver.set_initial_value(Y,0).set_f_params(ambient,h_w,MELT)
#[const.G,z_max,const.LAMBDA1,const.LAMBDA2,const.LAMBDA3,const.GAMT,const.GAMS,const.C_W,const.C_I,const.RHO_REF,const.ICETEMP,const.L])




    # iterate over depths while the solver completes successfully
    # and we are below the maximum height
    i = 0
    while solver.successful() and solver.t < h_w and solver.y[0]>0.:
        i += 1
        #print solver.t,z_range[i]
        solver.integrate(z_range[i])
        #solver.integrate(solver.t+z_step)
        #print solver.y
        # populate output arrays with solver result for current step
        rProfPlume[i] = solver.y[0]
        wProfPlume[i] = solver.y[1]
        tProfPlume[i] = solver.y[2]
        sProfPlume[i] = solver.y[3]
        aProfPlume[i] = solver.y[4]
        tProfAmb[i] = solver.y[5]
        sProfAmb[i] = solver.y[6]
        mIntProfPlume[i] = solver.y[7]

        thetaProfPlume[i] = solver.y[8]
        distanceProfPlume[i] = solver.y[9]
        #print tProfPlume[i],sProfPlume[i]
    #if gsw

    return [z_range,rProfPlume,wProfPlume,tProfPlume,sProfPlume,aProfPlume,tProfAmb,sProfAmb,
            mIntProfPlume,thetaProfPlume,distanceProfPlume]
    
'''
[d,r_p,w_p,t_p,s_p,a_p,t_a,s_a,melt,theta,distance] = calc_plume(w_sg, r_sg, z_max)

# plot some things
f,ax = plt.subplots(1,5,sharey=True, figsize=(10,5))

print "Terminal Height: ", max(d[np.isnan(r_p)])

#ax[0].plot(rProfPlume,depths,color='r')
ax[0].plot(r_p,d,color='r')
ax[0].set_title('Radius')
ax[0].set_xlim(0, np.nanmax(r_p)*1.1)
ax[0].set_ylim(0,100)
ax[0].set_xlabel('(m)')

ax[0].plot([0, np.nanmax(r_p)*1.1], [max(d[~np.isnan(r_p)]),max(d[~np.isnan(r_p)])])
ax[0].annotate('NB = %sm'%max(d[~np.isnan(r_p)]),xy=(np.nanmax(r_p)*0.05,max(d[~np.isnan(r_p)])+1))

ax[1].plot(w_p,d,color='b')
ax[1].set_title('Velocity')
ax[1].set_xlim(0, np.nanmax(w_p)*1.1)
ax[1].set_xlabel('(m/s)')

ax[2].plot(t_p,d,color='g')
ax[2].plot(t_p,d,ls='--',color='g')
ax[2].set_title('Temp')
ax[2].set_xlim(0, np.nanmax(t_p)*1.1)
ax[2].set_xlabel('($^o$C)')

ax[3].plot(s_p,d,color='k')
ax[3].plot(s_p,d,ls='--',color='k')
ax[3].set_title('Sal')
ax[3].set_xlim(0, np.nanmax(s_p)*1.1)
ax[3].set_xlabel('(g/kg)')

ax[4].plot(melt*3600.*24.,d,color='k')
#ax[4].plot(2.*rProfPlume*mIntProfPlume,depths,ls='--',color='k')
ax[4].set_title('Melt')
ax[4].set_xlim(0, np.nanmax(melt*3600.*24.)*1.1)
ax[4].set_xlabel('(m/d)')

if detached:
    ax[5].plot(distance,d,color='k')
    ax[5].set_title('"offline" melt')


#############################################   

# calculate an anlytical solution for the velocity from the same initial conditions

pressure = 1027.*const.G*500.
rho_p = gsw.rho(0.,0.,pressure*1.e-4)
rho_a = gsw.rho(34.,4.,pressure*1.e-4)
g_dash = -const.G*(rho_p-rho_a)/const.RHO_REF

B=g_dash*w_sg*math.pi*r_sg**2.

w = (5./(6.*const.E_0))*((9./10.)*const.E_0*B)**(1./3.)*math.pi**(-1./3.)*d**(-1./3.)

ax[1].plot(w,d,'--')


############################################

plt.savefig('iceplume.jpg')

plt.show()


'''

