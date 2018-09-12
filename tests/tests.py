#/usr/bin/env python

from pyplume import pyplume
#import ambient#test_data
import math
import numpy as np
import matplotlib.pyplot as plt
#import const
import gsw

from IPython import embed

# Variables:

# to move
#T_sg = 1.0e-3    # inlet temperature
#S_sg = 1.0e-3    # inlet salinity


#w_sg = .1#.434411779623        # inlet velocity
#r_sg = .1#1.#1.2105682097        # inlet radius

# to keep

h_w = 70.          # water depth
h_i = 100.          # ice depth

Tambient = 3.       # ambient salinity
Sambient = 33.      # ambient salinity

q = 100.            # discharge

# Ambient T/S profile

#new_amb = pyplume.InitAmbient(h_w,ambient.sal,ambient.temp,pressure=ambient.pressure)

pressure = np.arange(0,h_w+10.,1.)
new_amb = pyplume.InitAmbient(h_w,np.ones(pressure.shape)*Sambient,
                                  np.ones(pressure.shape)*Tambient,
                                  pressure=pressure)

r, u = pyplume.inlet(h_i,h_w,q)

out = pyplume.calc_plume(u,r, h_w, new_amb,MELT=False)

[z,r_p,w_p,t_p,s_p,a_p,t_a,s_a,melt,theta,distance] = out

#plt.plot(s_a,d,'--')
#plt.show()
plt.figure()
plt.plot(s_a,z)

# plot some things
f,ax = plt.subplots(1,5,sharey=True, figsize=(10,5))
embed()

if np.any(np.isnan(r_p)):
    nb_depth = max(z[~np.isnan(r_p)])
    print "Terminal Height: ", nb_depth
    print "plume density: ",gsw.rho(s_p[z==nb_depth],t_p[z==nb_depth],                                            new_amb.get_pres_z(nb_depth))
    print "ambient density: ",new_amb.get_rho_z(nb_depth)

else:
    plume_density = gsw.rho(s_p[-1],t_p[-1],
                            new_amb.get_pres_z(max(z)))
    amb_density = new_amb.get_rho_z(max(z))
    print "Excess buoyancy: ", pyplume.const.G*((plume_density-amb_density)/pyplume.const.RHO_REF)
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
print r_sg, u_sg
pressure = new_amb.get_pres_z(0)#1027.*const.G*
print "zeros pressure:",pressure,new_amb.get_pres_d(h_w)
rho_p = gsw.rho(0.,0.,pressure)
rho_a = new_amb.get_rho_z(0)#sw.rho(Sambient,Tambient,pressure)
print rho_a, new_amb.get_rho_z(0),pyplume.const.RHO_REF
g_dash = -pyplume.const.G*(rho_p-rho_a)/pyplume.const.RHO_REF
print len(rho_p),len(rho_a)
B=g_dash*u_sg*math.pi*r_sg**2.

w = (5./(6.*pyplume.const.E_0))*((9./10.)*pyplume.const.E_0*B)**(1./3.)*math.pi**(-1./3.)*z**(-1./3.)

ax[1].plot(w,z,'--')


############################################

plt.savefig('iceplume.jpg')

plt.show()




