#!/bin/python

import math

# Standard Constants

G = 9.81                    # acceleration due to gravity


# Ice/Water Properties

T_I = -0.                   # Ice temperature

RHO_I = 917.                # Ice density
RHO_W = 1000.               # Freshwater density
RHO_REF = 1020.0            # Reference density

C_W = 3994.                 # Specific heat capacity of seawater
C_I = 2009.                 # Specific heat capacity of ice
L = 334.e3                  # Latent heat of fusion
A = -0.0573                 # Liquidus slope
B = 0.0832                  # Liquidus intercept
C = 0.000761                # Liquidus pressure coefficient
GLEN_A = 6.e-24             # Flow law creep parameter
GLEN_N = 3.                 # Flow law creep exponent

# Plume/Turbulent Transfer Properties

E_0 = 0.09#832#0.08#1#05                  # Entrainment coefficient
GAM_T = 0.022               # Turbulent heat transfer coefficient
GAM_S = 0.00062             # Turbulent salinity transfer coefficient
C_D = 0.0025                # Drag coefficient
BACKGROUNDVEL = 0.001       # Background velocity

F = 0.1                     # Friction factor


# Derived constants

C1 = 1./(RHO_I*L)
C2 = 2.*GLEN_A*GLEN_N**(-GLEN_N)
C3 = 2.**0.25*(math.pi+2.)**0.5/(math.pi**0.25*(RHO_W*F)**0.5)
