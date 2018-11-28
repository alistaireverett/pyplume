#/usr/bin/env python

############################################################################################
#### !!!! THIS CODE IS PROBABLY BROKEN AT THE MOMENT, BUT WILL BE UPDATED SHORTLY !!!! #####
############################################################################################

import math
import numpy as np
import pandas as pd
import datetime as dt

import gsw
from pyplume import pyplume

from IPython import embed

# Variables:

h_w = 70.	# water depth at the terminus (ie. how high the plume will rise)
h_i = 100.	# ice thickness (used to calculate effective pressure for inlet conditions)

# Ambient T/S profile

# load a csv containing columns for 
# z (depth - metres) 
# SA (Absolute Salinity - Kg/m3)
# CT (Conservative temperature - C)

# This one was generated using the script in the ctdtools repo here:
# https://github.com/alistaireverett/ctdtools/blob/master/examples/proc_cnv_short.py

ambient = pd.read_csv('ambient.csv',';')

##########################
# start hack
##########################

# This is a short fix because our ambient profile doesn't start at zero
# this should be moved into the ctdtools preprocessing or into the InitAmbient class
# ignore this if your ambient profile already starts from zero

# TODO z should really be renamed to depth (at the start)

# copy the top row
shallowest_data = pd.DataFrame(ambient.iloc[0])

# set the depth to zero, include a catch 'if' to check whether z is a column or index
if 'z' in shallowest_data.index:
    shallowest_data.loc['z'] = 0.
else:
    shallowest_data.columns = [0.]

# append back onto ambient data frame
ambient = pd.concat([pd.DataFrame(shallowest_data).T,ambient],ignore_index=True)

##########################
# end hack
##########################

# create an ambient profile object from the temp/sal data
new_amb = pyplume.Ambient(h_w, ambient['SA'], ambient['CT'], depth=ambient['z'])


# Load discharges

# These are from the model presented in Pramanik et al. (2018) DOI: 10.1017/jog.2018.80
# Contains the following columns:
# Year, Month, Day, Hour, Discharge (m3/s)

disch_data = pd.read_csv('kronebreen_runoff.csv',';')

disch_data = disch_data[:100]

def to_dt(x):
    """
    Function to convert Y/M/D/H columns to python datetime object

    """
    return dt.datetime(int(x['year']),int(x['month']),int(x['day']),int(x['hour']))

# Add pandas datetime column using date conversion function to DataFrame
disch_data['date'] = disch_data.apply(to_dt,axis=1)

# Calculate inlet radius and velocity from discharge (for every discharge)
# Method follows Slater et al. (2015) / Schoof (2010)
rad, vel = pyplume.inlet(h_i,h_w,disch_data['discharge'])

# Add to new columns of DataFrame
disch_data['src_radius'] = rad
disch_data['src_vel'] = vel

# Getting to the good bit
# Now we can loop through all of the inlet properties and calculate a steady state
# solution for each given discharge (now defined by inlet radius and velocity)
res_data = []
for r, u in disch_data[["src_radius","src_vel"]].values:
	# The key line -  get steady state plume solution for given inputs
    out = pyplume.calc_plume(u, r, h_w, new_amb, MELT=True)

    '''
	# build dictionary from results - 
	# can (and should) be moved into the calc_plume func
    results_dict = {'z_range': out[0],
                    'rProfPlume': out[1],
                    'wProfPlume': out[2],
                    'tProfPlume': out[3],
                    'sProfPlume': out[4],
                    'aProfPlume': out[5],
                    'tProfAmb': out[6],
                    'sProfAmb': out[7],
                    'mIntProfPlume': out[8],
                    'thetaProfPlume': out[9],
                    'distanceProfPlume': out[10]}
    '''
    res_data.append(out)

### Now post process the results

# convert results to dataframe
#res_df = pd.DataFrame(res_data, index = disch_data.index)

# append to discharge dataframe
disch_data['results'] = res_data#res_df

# option to save the results now
disch_data.to_pickle('temp_save.pkl')

# calculate plume density profile from temp, sal and depth
disch_data['plume_density'] = disch_data['results'].apply(
                                lambda x: gsw.rho(x['s_p'],
                                                  x['t_p'],
                                                  new_amb.get_pres_z(x['z'])))
# Define a function to extract volume a given 
def vol_flux_at_depth(data, height_above_source):
    vol_flux = 0.5 * math.pi * data['b_p']**2. *data['u_p']
    vol_flux_out = vol_flux[data['z']==height_above_source]
    return vol_flux_out[0]

# apply function to get vol flux 1 metre below the surface
disch_data['vol_flux'] = disch_data['results'].apply(vol_flux_at_depth, 
													 height_above_source = h_w - 1.)

# calculate difference between surface flux and discharge as entrained flux
disch_data['ent_flux'] = disch_data['vol_flux'] - disch_data['discharge']

# calc_plume returns NaN if the plume is above the neutral buoyancy height
# so extract nb_height as the max height where the plume radius is not NaN
disch_data['nb_depth'] = disch_data['results'].apply(lambda x: max(x['z'][~np.isnan(x['b_p'])]))

### Write results out to a .csv

# Pull relevant columns from disch_data dataframe
csv_out = disch_data[['year','month','day','hour','discharge','ent_flux','nb_depth']]

# Rename columns to something a bit friendlier
csv_out.rename(columns={'discharge':	'Discharge Flux (m3/s)',
						'ent_flux':		'Entrained Flux (m3/s)',
						'nb_depth': 	'Plume Height (m)'},
						inplace=True)

# Add total volumes (remember, inputs are at 6-hour intervals) 
csv_out['Discharge Vol (m3)'] = csv_out['Discharge Flux (m3/s)'] * 60. * 60. * 6.

csv_out['Entrained Vol (m3)'] = csv_out['Entrained Flux (m3/s)'] * 60. * 60. * 6.

csv_out.fillna(0.,inplace=True)

csv_out.to_csv('upwelling_results.csv', ';')

# save data in full for future analysis
disch_data.to_pickle('upwelling.pkl')

