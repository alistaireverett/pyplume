# pyplume v1.0

Pyplume calculates a numerical solution to the equations formulated for a
point source plume released at the base of a vertical wall. The code includes
capabilities for calculating the melt rate in the case where the vertical wall
is made of ice.

See the examples and docstrings for more details on how to use pyplume. The links below will launch interactive jupyter notebooks with Binder so you can see how the code
works and experiment with the examples.

### Simple plume

A minimal example demonstrating a simple case to use pyplume with linear ambient
temperature and salinity profiles and plot the results.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/alistaireverett/pyplume/master?filepath=examples%2Fsimple_plume.ipynb)

### Comparison to analytical solution

A slightly more detailed look at using pyplume and comparing the results to an
analytical solution.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/alistaireverett/pyplume/master?filepath=examples%2Fpyplume_vs_analytical.ipynb)

### Kronebreen Upwelling

An example looking at how to run pyplume using a year long time series of
discharges calculated using a surface mass balance model and ambient conditions
taken from CTD casts.

***To be added soon***

More examples and documentation will be added shortly.

The code was developed by Alistair Everett with inspiration from Tom Cowton's iceplume Fortran model (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014JC010324)

