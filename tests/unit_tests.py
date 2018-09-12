#!/usr/bin/env python

import numpy as np
from pyplume import pyplume
import gsw
import unittest

def ambient_profile_test():
    
    z_max = 70.
    salinity = np.array([30.,34.])
    temperature = np.array([0.,2.])
    pressure = np.array([0.,100.])

    ambient = pyplume.InitAmbient(z_max,salinity,temperature,pressure=pressure)

    assert ambient.rho[0] == gsw.rho(salinity[0],temperature[0],pressure[0]), "Density not calculated correctly"
    assert ambient.rho[1] == gsw.rho(salinity[1],temperature[1],pressure[1]), "Density not calculated correctly"

    mid_rho = (ambient.rho[0]+ambient.rho[1])/2.
    mid_sal = (salinity[0]+salinity[1])/2.
    mid_temp = (temperature[0]+temperature[1])/2.
    mid_pressure = (pressure[0]+pressure[1])/2.

    mid_depth = mid_pressure/(1027.*9.81*1.e-4)


    assert mid_sal == ambient.get_sal_d(mid_depth), "Something wrong with salinity interpolation"


if __name__ == "__main__":

    ambient_profile_test()
