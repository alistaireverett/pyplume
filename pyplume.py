#!/usr/bin/python

import numpy as np
import math
from scipy.integrate import ode
from scipy.interpolate import interp1d
import gsw
import const

class Ambient(object):
    """Class object to store ambient temperature and salinity profiles

    Used to store ambient temperature and salinity profiles.
    Automatically derives a density profile and either pressure or depth
    depending on which the object is initialised with.
    Also provides convenience functions to extract interpolated values from
    any depth in the profile by specifying either the depth or height in the profile.

    NB: You must specify either a pressure or a depth profile. 
    If you specify both, take note of the assumption used here to relate pressure
    and depth which other associated functions may also rely on.

    depth = 0 at surface, increasing down
    z = z_max at surface, decreasing down (z = z_max - depth)

    Args:
        z_max (float) :                     maximum depth/height (same difference!)
        salinity (np.array) :               salinity profile (kg/m3)
        temperature (np.array) :            temperature profile (deg C)
        pressure (np.array , optional) :    pressure profile (dbar)
        depth (np.array, optional) :        depth below surface (m)

    Attributes:
        depth :         depth below surface (m)
        z_max :         maximum depth/height (same difference!) (m)
        z :             depth profile (m)
        salinity :      salinity profile (kg/m3)
        temperature :   temperature profile (deg C)
        rho :           density profile determined using gsw TEOS-10 package (kg/m3)
	
    Raises:
        ValueError :    if neither pressure or depth is specified
        UserWarning :   if pressure and depth are specified to don't match model assumptions
    """

    
    def __init__(self, z_max, salinity, temperature, pressure=None, depth=None):

        if depth is None and pressure is None:
            raise ValueError("Must pass either pressure or depth as an argument to Ambient")
        elif depth is None and pressure is not None:
            if min(pressure) != 0.:
                raise ValueError("Profile must start from the surface for the " \
                "interpolation, try copying the minimum depth values to be the " \
                "surface values")
            self.pressure = pressure
            self.depth = pressure/(1027.*9.81*1.e-4)
        elif depth is not None and pressure is None:
            if min(depth) != 0:
                raise ValueError("Profile must start from the surface for the " \
                "interpolation, try copying the minimum depth values to be the " \
                "surface values")

            self.pressure = depth*(1027.*9.81*1.e-4)
            self.depth = depth
        else:
            if np.any(depth != pressure/(1027.*9.81*1.e-4)):
                print("UserWarning: You specified a pressure and depth profile "\
                "which may not match assumptions used elsewhere in this model")
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
        """Return salinity float for desired depth below surface"""
        return self.__f_sal_d(x)[()]

    def get_temp_d(self,x):
        """Return temperature float for desired depth below surface"""
        return self.__f_temp_d(x)[()]

    def get_rho_d(self,x):
        """Return density float for desired depth below surface"""
        return self.__f_rho_d(x)[()]

    def get_pres_d(self,x):
        """Return pressure float for desired depth below surface"""
        return self.__f_pres_d(x)[()]

    def get_sal_z(self,x):
        """Return salinity float for desired height above botttom"""
        return self.__f_sal_d(self.z_max - x)[()]

    def get_temp_z(self,x):
        """Return temperature float for desired height above botttom"""
        return self.__f_temp_d(self.z_max - x)[()]
                    
    def get_rho_z(self,x):
        """Return density float for desired height above botttom"""
        return self.__f_rho_d(self.z_max - x)[()]
                    
    def get_pres_z(self,x):
        """Return pressure float for desired height above botttom"""
        return self.__f_pres_d(self.z_max - x)[()]
                    

def inlet(h_i, h_w, q):
    """Calculate radius and velocity at source

    Function to calculate inlet velocity and radius from effective pressure
    Assumes semi-circular cross-section
    See: Eqn S4, Everett et al. (2018) DOI: 10.1038/s41598-018-31875-8

    Args:
        h_i (float) :   ice thickness (m)
        h_w (float) :   water depth (m)
        q (float) :     discharge (m3/s)

    Returns:
        radius (float), velocity (float)

    Raises:
        AssertionError: if values produce effective pressure greater than zero
    """

    n_eff = (const.RHO_I*const.G*h_i-const.RHO_W*const.G*h_w)
    assert n_eff > 0., "Terminus is floating - need shallower water or more ice!"
    factor = (const.C1/(const.C2*const.C3**2.*n_eff**const.GLEN_N))**(2./7.)
    
    area = factor*q**(6./7.)
    radius = (2*area/math.pi)**0.5
    velocity = q/area

    return radius, velocity


def get_melt(v_inf,t_inf,s_inf,pressure):
    """Quadratic solution to the three equation melt parameterisation

    Args:
        v_inf (float) :     farfield velocity (m/s)
        t_inf (float) :     farfield temperature (deg C)
        s_inf (float) :     farfield salinity (kg/m3)
        pressure (float) :  pressure (Pa)

    Returns:
        t_b (float) :   boundary temperature (deg C)
        s_b (float) :   boundary salinity (kg/m3)
        mdot (float) :  melt rate (m/s)
    """
    # solve the three equation melt-parameterisation quadratically
    aa = const.A*(const.GAM_T*const.C_W - const.GAM_S*const.C_I)
    bb = (const.GAM_S*const.C_I*
             (const.A*s_inf - const.B - const.C*pressure + const.T_I 
                - (const.L/const.C_I)) 
             - const.GAM_T*const.C_W*(t_inf - const.B - const.C*pressure))
    cc = const.GAM_S*s_inf*(const.C_I*(const.B + const.C*pressure - const.T_I)                             + const.L)

    s_b = (1./(2.*aa))*(-bb-((bb**2.-4.*aa*cc)**0.5))
    t_b = const.A*s_b+const.B+const.C*pressure
    mdot = const.GAM_S*(const.C_D**0.5)*v_inf*(s_inf-s_b)/s_b

    return t_b, s_b, mdot

def wallPlume(z, y, ambient, z_max, MELT=True):
    """Solve the equations for a wallPlume

    wallPlume formulation (halfCone in Cowton et al.)
    See: Cowton et al. (2015) DOI: 10.1002/2014JC010324
    """

    # this was a safety check at some point - is it still needed?!
    if z > z_max:
        return None

    # initialise array for output
    ydot = np.zeros(y.shape)
    
    # calculate melt rate if required
    if MELT:
        t_b, s_b, mdot = get_melt(y[1], y[2], y[3], ambient.get_pres_z(z))
    else:
        t_b = 0.
        s_b = 0.
        mdot = 0.

    # get ambient conditions at whatever depth we're at
    t_amb = ambient.get_temp_z(z)
    s_amb = ambient.get_sal_z(z)
    rho_a = ambient.get_rho_z(z)
    
    # approximate pressure at the current depth in dbar 
    pressure = ambient.get_pres_z(z)

    # calculate current plume density (needs pressure in decibar)
    # gives density in kg/m3
    rho_p = gsw.rho(y[3], y[2], pressure)
    
    # check if Neutral Buoyancy is reached, if so this forces values to be nan 
    if rho_p > rho_a:
        y[0] = np.nan
        y[1] = np.nan
        
    # Solve the plume equations and store in ydot
    ydot[0] = (2.*const.E_0 + 4.*mdot/(math.pi*y[1])
               - y[0]*const.G*(rho_a-rho_p)/(2.*y[1]*y[1]*const.RHO_REF)
               + 2. * (const.C_D/math.pi))

    ydot[1] = (- 2.*const.E_0*y[1]/y[0] - 4.*mdot/(math.pi*y[0])
               + const.G*(rho_a-rho_p)/(y[1]*const.RHO_REF)
               - 4.*const.C_D*y[1]/(math.pi*y[0]))

    ydot[2] = (2.*const.E_0*(t_amb-y[2])/y[0]
               + 4.*mdot*(t_b-y[2])/(math.pi*y[0]*y[1])
               - 4.*const.GAM_T*(const.C_D**0.5)*(y[2]-t_b)/(math.pi*y[0]))

    ydot[3] = (2.*const.E_0*(s_amb-y[3])/y[0] 
               + 4.*mdot*(s_b-y[3])/(math.pi*y[0]*y[1])
               - 4.*const.GAM_S*(const.C_D**0.5)*(y[3]-s_b)/(math.pi*y[0]))

    ydot[4] = mdot - y[4]

    ydot[5] = t_amb - y[5]
    ydot[6] = s_amb - y[6]

    return ydot

def calc_plume(u_0, b_0, h_w, ambient, t_0 = 1.0e-3, s_0 = 1.0e-3, MELT=True):
    """Solve the plume equations for the specified initial conditions


    Args:
        u_0 (float) :               velocity at the source (m/s)
        b_0 (float) :               radius at the source (m)
        h_w (float) :               total water column depth (m)
        ambient (Ambient object) :  containing temperature and salinity profiles
        t_0 (float, optional) :     temperature at the plume source (deg C, default=1.0e-3)
        s_0 (float, optional) :     salinity at the plume source (kg/m3, default=1.0e-3)
        MELT (boolean, optional) :  whether to include the melt feedback (default=True)

    Returns:
        plume (dict) : dictionary containing:
                            'z'     height above source (m)
                            'b_p'   plume radius (m)
                            'w_p'   vertical velocity (m/s)
                            't_p'   temperature (deg C)
                            's_p'   salinity (kg/m3)
                            'm_p'   melt rate (m/s)
                            't_a'   ambient salinity (kg/m3)
                            's_a'   ambient temperature (deg C)
    """
    
    # TODO: The initialisation could still be made a bit tider/clearer

    plume_variables = ['b_p', 	# plume radius (m)
                       'w_p', 	# vertical velocity (m/s)
                       't_p', 	# temperature (deg C)
                       's_p', 	# salinity (kg/m3)
                       'm_p',   # melt rate (m/s)
                       't_a', 	# ambient salinity (kg/m3)
                       's_a', 	# ambient temperature (deg C)
                      ]

    # Create output dictionary containting all the keys in the plume variables
    # list
    plume = {key: [] for key in plume_variables}
    
    # Populate the output dict with the initial conditions
    plume['b_p'].append(b_0)
    plume['w_p'].append(u_0)
    plume['t_p'].append(t_0)
    plume['s_p'].append(s_0)
    plume['m_p'].append(0.)
    plume['t_a'].append(ambient.get_temp_z(0))
    plume['s_a'].append(ambient.get_sal_z(0))

    # extract the initial conditions as an input array for the ode solver
    y = np.array([plume[key][0] for key in plume_variables])

    # set up the depths to iterate over
    # depth = 0 at surface, increasing down
    # z = z_max at surface, decreasing down (z = z_max - depth)
    z_step = .1
    z_range = np.arange(0, h_w + z_step, z_step)

    # Set the solver and tolerances
    solver = ode(wallPlume).set_integrator('lsoda', atol=1.e-5, rtol=1.e-5)
    
    # set initial conidtions in solver with y array, z=0 and args
    solver.set_initial_value(y, 0).set_f_params(ambient, h_w, MELT)

    # iterate over depths while the solver completes successfully
    # and we are below the maximum height
    i = 0
    while solver.successful() and solver.t < h_w and ~np.isnan(solver.y[0]):
        i += 1

        # do the integration
        solver.integrate(z_range[i])
        
        # populate output dict with solver result for current step
        [plume[key].append(solver.y[k].item()) for k,key in enumerate(plume_variables)]

    # add z_range into output dict
    plume['z'] = z_range

    # populate rest of dict with nans if they solver didn't reach the surface
    for key in plume_variables:
        plume[key] += [np.nan]*(len(plume['z'])-len(plume[key]))

    return plume 
    

