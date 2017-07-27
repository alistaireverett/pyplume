#!/usr/bin/env python

import numpy as np
import pyplume
import gsw
import unittest


class PyPlumeTests(unittest.TestCase):

    def setUp(self):
        self.salinity = np.array([30.,34.])
        self.temperature = np.array([0.,2.])
        self.pressure = np.array([0.,100.])
        self.depth = self.pressure/(1027.*9.81*1.e-4)
        self.z_max = 1.
        #self.ambient = pyplume.InitAmbient(self.salinity,self.temperature,depth=self.depth)


    def test_density_from_pressure(self):
        self.ambient = pyplume.InitAmbient(self.z_max,self.salinity,self.temperature,pressure=self.pressure)
        self.assertEqual(self.ambient.rho[0],gsw.rho(self.salinity[0],self.temperature[0],self.pressure[0]))
		#self.assertEqual(self.ambient.rho[1],gsw.rho(self.salinity[1],self.temperature[1],self.pressure[1]))
 
    def test_density_from_depth(self):
        self.ambient = pyplume.InitAmbient(self.z_max,self.salinity,self.temperature,depth=self.depth)
        self.assertEqual(self.ambient.rho[0],gsw.rho(self.salinity[0],self.temperature[0],self.pressure[0]))

    def test_profile_interpolation(self):
        self.ambient = pyplume.InitAmbient(self.z_max,self.salinity,self.temperature,pressure=self.pressure)

        mid_rho = (self.ambient.rho[0]+self.ambient.rho[1])/2.
        mid_sal = (self.salinity[0]+self.salinity[1])/2.
        mid_temp = (self.temperature[0]+self.temperature[1])/2.
        mid_pressure = (self.pressure[0]+self.pressure[1])/2.

        mid_depth = mid_pressure/(1027.*9.81*1.e-4)
        
        self.assertEqual(mid_sal, self.ambient.get_sal_d(mid_depth))

    def test_no_args(self):
		self.assertRaises(ValueError, pyplume.InitAmbient,self.z_max,self.salinity,self.temperature)

    def test_zmax_vs_d0(self):
        self.ambient = pyplume.InitAmbient(self.z_max,self.salinity,self.temperature,pressure=self.pressure)
        self.assertEqual(self.ambient.get_sal_z(self.z_max),self.ambient.get_sal_d(0))

    def test_z0_vs_dmax(self):
        self.ambient = pyplume.InitAmbient(self.z_max,self.salinity,self.temperature,pressure=self.pressure)
        self.assertEqual(self.ambient.get_sal_z(0),self.ambient.get_sal_d(self.z_max))

if __name__ == "__main__":

    unittest.main()
