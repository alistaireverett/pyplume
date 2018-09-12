#/usr/bin/env python

import math
import numpy as np
import matplotlib.pyplot as plt

import pyplume
import const
import ambient


h_w = 70.          # water depth
h_i = 100.          # ice depth

Tambient = 3.8       # ambient salinity
Sambient = 34.      # ambient salinity


#pressure = np.arange(0,h_w+10.,1.)
#sal_array = np.ones(pressure.shape)*Sambient#
#sal_array = np.linspace(33.,34.,len(pressure))
#temp_array = np.ones(pressure.shape)*Tambient

#start_depth = h_w[:]

nb_depth = 0
while nb_depth < 70.:
    
    # Ambient T/S profile

    new_amb = pyplume.InitAmbient(h_w-nb_depth,
                                #sal_array,temp_array,pressure=pressure,
                                ambient.sal,ambient.temp,pressure=ambient.pressure
                                )

    out = pyplume.calc_plume(1.e-6, 1.e-6, h_w-nb_depth, new_amb,MELT=True)

    [z,r_p,w_p,t_p,s_p,a_p,t_a,s_a,melt,theta,distance] = out

    plt.figure(1)
    plt.plot(w_p,z+nb_depth)

    plt.figure(2)
    plt.plot(s_a,z+nb_depth)

    plt.figure(3)
    plt.plot(t_a,z+nb_depth)

    nb_depth = nb_depth + max(z[~np.isnan(r_p)])

    print nb_depth

    #h_w = h_w-nb_depth          # water depth
    #h_i = h_i-nb_depth          # ice depth

    print h_w,h_i



'''
new_amb = pyplume.InitAmbient(h_w,ambient.sal,ambient.temp,
                                pressure=ambient.pressure)

out = pyplume.calc_plume(1.e-6,1.e-6, h_w, new_amb,MELT=True)

[z2,r_p2,w_p2,t_p2,s_p2,a_p2,t_a2,s_a2,melt2,theta2,distance2] = out

nb_depth = max(z[~np.isnan(r_p)])

print nb_depth

plt.plot(w_p,z)
plt.plot(w_p2,z2+nb_depth)
'''
plt.show()

