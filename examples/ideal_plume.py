#/usr/bin/env python

from pyplume import pyplume, const
import math
import numpy as np
import matplotlib.pyplot as plt
#import const
import gsw

from IPython import embed

# Variables:

# to keep

h_w = 700.          # water depth
h_i = 780.          # ice depth

Tambient = 3.       # ambient salinity
Sambient = 33.      # ambient salinity

q = 10000.            # discharge

# Ambient T/S profile

pressure = np.arange(0,h_w+10.,1.)
sal_array = np.linspace(30.,33.,len(pressure))

new_amb = pyplume.InitAmbient(h_w,sal_array,#np.ones(pressure.shape)*Sambient,
                                  np.ones(pressure.shape)*Tambient,
                                  pressure)


for q in np.array([1,]):#2,5,10,20,50,100,200,500,1000]):#range(100,1000,100):
    out = pyplume.calc_plume(q, h_i, h_w, new_amb,MELT=True)

    [z,r_p,w_p,t_p,s_p,a_p,t_a,s_a,melt,theta,distance] = out

    plt.figure(10)
    plt.plot(t_a,z)
    plt.figure(11)
    plt.plot(s_a,z)

    if np.any(np.isnan(r_p)):
        nb_depth = max(z[~np.isnan(r_p)])
        print "Terminal Height: ", nb_depth
        print "plume density: ",gsw.rho(s_p[z==nb_depth],t_p[z==nb_depth],                                            new_amb.get_pres_z(nb_depth))
        print "ambient density: ",new_amb.get_rho_z(nb_depth)

    else:
        plume_density = gsw.rho(s_p[-1],t_p[-1],
								new_amb.get_pres_z(max(z)))
        amb_density = new_amb.get_rho_z(max(z))
        g_dash = -const.G*((plume_density-amb_density)/const.RHO_REF)
        print "plume density: ",plume_density
        print "ambient density: ",amb_density
		
        print "ref density: ",const.RHO_REF
        print "Excess buoyancy: ", g_dash
        nb_depth = np.nan

    plt.figure(1)
    plt.scatter(q,g_dash)
	
    plt.figure(2)
    plt.scatter(q,r_p[-1])

    plt.figure(3)
    plt.scatter(q,s_p[-1])
    
    plt.figure(4)
    plt.scatter(q,math.pi*r_p[-1]**2.*w_p[-1]**2.)

    plt.figure(5)
    plt.scatter(q,math.pi*r_p[-1]**2.*w_p[-1]*g_dash)

#plt.plot(s_a,d,'--')
#plt.show()
plt.figure()
plt.plot(s_a,z)

# plot some things
f,ax = plt.subplots(1,5,sharey=True, figsize=(10,5))

if np.any(np.isnan(r_p)):
    nb_depth = max(z[~np.isnan(r_p)])
    print "Terminal Height: ", nb_depth
    print "plume density: ",gsw.rho(s_p[z==nb_depth],t_p[z==nb_depth],                                            new_amb.get_pres_z(nb_depth))
    print "ambient density: ",new_amb.get_rho_z(nb_depth)

else:
    plume_density = gsw.rho(s_p[-1],t_p[-1],
                            new_amb.get_pres_z(max(z)))
    amb_density = new_amb.get_rho_z(max(z))
    print "plume density: ",plume_density
    print "ambient density: ",amb_density
    print "ref density: ",const.RHO_REF
    print "Excess buoyancy: ", const.G*((plume_density-amb_density)/const.RHO_REF)
    nb_depth = np.nan
    
#ax[0].plot(rProfPlume,depths,color='r')
ax[0].plot(r_p,z,color='r')
ax[0].set_title('Radius')
ax[0].set_xlim(0, np.nanmax(r_p)*1.1)
ax[0].set_ylim(0,h_w)
ax[0].set_xlabel('(m)')

ax[0].plot([0, np.nanmax(r_p)*1.1], [max(z[~np.isnan(r_p)]),max(z[~np.isnan(r_p)])])
ax[0].annotate('NB = %sm'%max(z[~np.isnan(r_p)]),xy=(np.nanmax(r_p)*0.05,max(z[~np.isnan(r_p)])+1))

ax[1].plot(w_p,z,color='b')
ax[1].set_title('Velocity')
ax[1].set_xlim(0, np.nanmax(w_p)*1.1)
ax[1].set_xlabel('(m/s)')

ax[2].plot(t_p,z,color='g')
ax[2].plot(t_a,z,ls='--',color='g')
ax[2].set_title('Temp')
ax[2].set_xlim(0, np.nanmax(t_p)*1.1)
ax[2].set_xlabel('($^o$C)')

ax[3].plot(s_p,z,color='k')
ax[3].plot(s_a,z,ls='--',color='k')
ax[3].set_title('Sal')
ax[3].set_xlim(0, np.nanmax(s_p)*1.1)
ax[3].set_xlabel('(g/kg)')

ax[4].plot(melt*3600.*24.,z,color='k')
#ax[4].plot(2.*rProfPlume*mIntProfPlume,depths,ls='--',color='k')
ax[4].set_title('Melt')
ax[4].set_xlim(0, np.nanmax(melt*3600.*24.)*1.1)
ax[4].set_xlabel('(m/d)')



#############################################   

# calculate an anlytical solution for the velocity from the same initial conditions

r_sg, u_sg = pyplume.inlet(h_i,h_w,q)

pressure = new_amb.get_pres_z(0)#1027.*const.G*
print "zeros pressure:",pressure,new_amb.get_pres_d(h_w)
rho_p = gsw.rho(0.,0.,pressure)
rho_a = new_amb.get_rho_z(0)#sw.rho(Sambient,Tambient,pressure)
print rho_a, new_amb.get_rho_z(0),const.RHO_REF
g_dash = -const.G*(rho_p-rho_a)/const.RHO_REF

B=g_dash*u_sg*math.pi*r_sg**2.

w = (5./(6.*const.E_0))*((9./10.)*const.E_0*B)**(1./3.)*math.pi**(-1./3.)*z**(-1./3.)

ax[1].plot(w,z,'--')


############################################

plt.savefig('iceplume.jpg')

plt.show()




