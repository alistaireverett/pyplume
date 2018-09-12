#/usr/bin/env python

import pyplume
import math
import numpy as np
import matplotlib.pyplot as plt
import const
import gsw
import ambient
import csv
import datetime as dt

from IPython import embed

def mlab2datetime(num):
    day = dt.datetime.fromordinal(int(num))
    dayfrac = dt.timedelta(days=num%1) - dt.timedelta(days=366)
    return day + dayfrac

# Variables:

# to keep

h_w = 70.          # water depth
h_i = 100.          # ice depth

# Ambient T/S profile

new_amb = pyplume.InitAmbient(h_w,ambient.sal,ambient.temp,
                                pressure=ambient.pressure)

# load discharges
times = []
disch_flux = []

with open('discharge_kronebreen.csv') as csvfile:
    readcsv = csv.reader(csvfile,delimiter=';')
    for row in readcsv:
        times.append(float(row[0]))
        disch_flux.append(float(row[1]))

times = np.array(times)
disch_flux = np.array(disch_flux)

str_times = [mlab2datetime(t).isoformat(' ') for t in times]

g_dash = []
nb_depth =[]
vol_flux = []
momentum = []
buoyancy = []
surface = []
ent_flux = []

count=0

cut = 10000

for t,q in zip(times,disch_flux):

    print "count: ",count,"time: ",t
    count+=1
    
    if count>cut: break
    if q < .1:
        g_dash.append(0)
        nb_depth.append(0)
        vol_flux.append(0)
        momentum.append(0)
        buoyancy.append(0)
        surface.append(0)
        ent_flux.append(0)
        continue

    r, u = pyplume.inlet(h_i,h_w,q)
    out = pyplume.calc_plume(u,r, h_w, new_amb,MELT=True)
    
    [z,r_p,w_p,t_p,s_p,a_p,t_a,s_a,melt,theta,distance] = out

    this_nb_depth = max(z[~np.isnan(r_p)])

    if this_nb_depth == max(z):
        plume_density = gsw.rho(s_p[-1],t_p[-1],
                                new_amb.get_pres_z(max(z)))
        amb_density = new_amb.get_rho_z(max(z))
        this_vol_flux = 0.5*math.pi*r_p[-1]**2.*w_p[-1]
        this_g_dash = -const.G*((plume_density-amb_density)/const.RHO_REF)
        this_momentum = 0.5*math.pi*r_p[-1]**2.*w_p[-1]**2.
        this_buoyancy = 0.5*math.pi*r_p[-1]**2.*w_p[-1]*this_g_dash
        this_surface = 1
        this_ent_flux = this_vol_flux - q
    else:
        this_vol_flux = 0.
        this_g_dash = 0.
        this_momentum = 0.
        this_buoyancy = 0.
        this_surface = 0
        this_ent_flux = 0.

    g_dash.append(this_g_dash)	
    nb_depth.append(this_nb_depth)   
    vol_flux.append(this_vol_flux)
    momentum.append(this_momentum)
    buoyancy.append(this_buoyancy)
    surface.append(this_surface)
    ent_flux.append(this_ent_flux)
'''
if vol_flux > 1.e-3:
    entrained_flux = np.array(vol_flux) - disch_flux[:cut]
else:
    entrained_flux = 0.
'''



# disharge is m3/s for a 3 hour period
entrained_vol = np.array(ent_flux) * 60. * 60. * 3.
#cumul_ent_vol = [] # np.cumsum(entrained_vol)
#for yr in np.unique(str_times.year):
#   cumul_ent_vol = np.concatenate((cumul_disch_vol,np.cumsum(this[this[:,0]==yr,1])))
#for yr in range(2010,2017,1):
#    this_yr = dt.datetime.toordinal(dt.date(yr,1,1))
#    cumul_ent_vol = np.concatenate((cumul_ent_vol,np.cumsum(entrained_vol[times[:cut]<this_yr])))

disch_vol = disch_flux[:cut] * 60. * 60. * 3.
#cumul_disch_vol = [] # np.cumsum(disch_vol)
#for yr in range(2010,2017,1):
#    this_yr = dt.datetime.toordinal(dt.date(yr,1,1))
#    cumul_disch_vol = np.concatenate((cumul_disch_vol,np.cumsum(disch_vol[times[:cut]<this_yr])))

np.savez('KbPlume',[str_times,disch_flux,surface,ent_flux,
                    disch_vol,entrained_vol,g_dash,nb_depth,vol_flux,momentum,buoyancy,surface])
 

with open('kb_upwelling.csv', 'wb') as csvout:
    print 'writing'
    csvwriter = csv.writer(csvout, delimiter=',')
    csvwriter.writerow(['Time','Discharge Flux (m3/s)','Surface?','Entrained Flux (m3/s)',
                        'Runoff Vol (m3)', 'Entrained Vol (m3)',
                        #'Cumulative Runoff (m3)','Cumulative Entrained Vol (m3)'
                        ])
    for row in zip(str_times,disch_flux,surface,ent_flux,
                    disch_vol,entrained_vol,
                    #cumul_disch_vol,cumul_ent_vol
                    ):
        csvwriter.writerow(row)

plt.figure()
plt.title('Discharge vs entrainment')
plt.plot(disch_flux,ent_flux)

plt.figure()
plt.title('Discharge')
plt.plot(times[:cut],disch_flux[:cut])

plt.figure()
plt.title('Volume Flux')
plt.plot(times[:cut],vol_flux[:cut])

plt.figure()
plt.title("g'")
plt.plot(times[:cut],g_dash)

plt.figure()
plt.title('Buoyancy')
plt.plot(times[:cut],buoyancy)

plt.figure()
plt.title('Momentum')
plt.plot(times[:cut],momentum)

plt.figure()
plt.title('Neutral Buoyancy Depth')
plt.plot(times[:cut],nb_depth)

plt.figure()
plt.title('Discharge')
plt.plot(times[:cut],disch_flux[:cut])

plt.figure()
plt.title('Reached the surface?')
plt.plot(times[:cut],surface)

plt.show()
#plt.plot(s_a,d,'--')
#plt.show()
#plt.figure()
#plt.plot(s_a,z)

# plot some things
f,ax = plt.subplots(1,5,sharey=True, figsize=(10,5))
'''
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
''' 
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




############################################

plt.savefig('iceplume.jpg')

plt.show()




