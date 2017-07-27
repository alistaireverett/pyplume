#/usr/bin/env python

import pyplume
import ambient#test_data
import math
import numpy as np
import matplotlib.pyplot as plt
import const
import gsw

# Variables:

T_sg = 1.0e-3    # inlet temperature
S_sg = 1.0e-3    # inlet salinity
w_sg = .1#.434411779623        # inlet velocity
r_sg = .1#1.#1.2105682097        # inlet radius
z_max = 1000.

# Detached P
theta_sg = math.pi/2.

Tambient = 5.       # ambient salinity
Sambient = 33.      # ambient salinity

detached = False

# Ambient T/S profile

z_max=100.

#new_amb = pyplume.InitAmbient(z_max,ambient.sal,ambient.temp,pressure=ambient.pressure)

pressure = np.arange(0,z_max+10,1)
new_amb = pyplume.InitAmbient(z_max,np.ones(pressure.shape)*Sambient,np.ones(pressure.shape)*Tambient,pressure)

'''
plt.figure()
#plt.plot(ambient.sal,ambient.pressure)
plt.plot(new_amb.salinity,new_amb.z)
for i in range(0,int(z_max)+1,1):
    plt.scatter(new_amb.get_sal_d(i),z_max-i,marker='+')
for i in range(0,int(z_max)+1,1):
    plt.scatter(new_amb.get_sal_z(i),i,marker='x')

plt.plot([30,35],[z_max,z_max])
'''
print z_max,len(ambient.sal_const),len(ambient.temp_const),len(ambient.pressure)
#print new_amb.get_pres_z(91.)
[d,r_p,w_p,t_p,s_p,a_p,t_a,s_a,melt,theta,distance] = pyplume.calc_plume(w_sg, r_sg, z_max, new_amb)

#plt.plot(s_a,d,'--')
#plt.show()
plt.figure()
plt.plot(s_a,d)

# plot some things
f,ax = plt.subplots(1,5,sharey=True, figsize=(10,5))

try: 
    print "Terminal Height: ", max(d[np.isnan(r_p)])
except:
    pass 

#ax[0].plot(rProfPlume,depths,color='r')
ax[0].plot(r_p,d,color='r')
ax[0].set_title('Radius')
ax[0].set_xlim(0, np.nanmax(r_p)*1.1)
ax[0].set_ylim(0,z_max)
ax[0].set_xlabel('(m)')

ax[0].plot([0, np.nanmax(r_p)*1.1], [max(d[~np.isnan(r_p)]),max(d[~np.isnan(r_p)])])
ax[0].annotate('NB = %sm'%max(d[~np.isnan(r_p)]),xy=(np.nanmax(r_p)*0.05,max(d[~np.isnan(r_p)])+1))

ax[1].plot(w_p,d,color='b')
ax[1].set_title('Velocity')
ax[1].set_xlim(0, np.nanmax(w_p)*1.1)
ax[1].set_xlabel('(m/s)')

ax[2].plot(t_p,d,color='g')
ax[2].plot(t_a,d,ls='--',color='g')
ax[2].set_title('Temp')
ax[2].set_xlim(0, np.nanmax(t_p)*1.1)
ax[2].set_xlabel('($^o$C)')

ax[3].plot(s_p,d,color='k')
ax[3].plot(s_a,d,ls='--',color='k')
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

pressure = new_amb.get_pres_z(0)#1027.*const.G*
rho_p = gsw.rho(0.,0.,pressure)
rho_a = gsw.rho(Sambient,Tambient,pressure)
print rho_a, new_amb.get_rho_z(0),const.RHO_REF
g_dash = -const.G*(rho_p-rho_a)/const.RHO_REF

B=g_dash*w_sg*math.pi*r_sg**2.

w = (5./(6.*const.E_0))*((9./10.)*const.E_0*B)**(1./3.)*math.pi**(-1./3.)*d**(-1./3.)

ax[1].plot(w,d,'--')


############################################

plt.savefig('iceplume.jpg')

plt.show()




