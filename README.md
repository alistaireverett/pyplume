# pyplume v1.0

Pyplume calculates a numerical solution to the equations formulated for a point source plume released at the base of a vertical wall. The code includes capabilities for calculating the melt rate in the case where the vertical wall is made of ice.

See the examples and docstrings for more details on how to use pyplume. The links below will launch interactive jupyter notebooks with [Binder](https://mybinder.readthedocs.io/en/latest/) so you can see how the code works and experiment with the examples. [Here's a link](https://medium.com/codingthesmartway-com-blog/getting-started-with-jupyter-notebook-for-python-4e7082bd5d46#f8b4) explaining how to interact with jupyter notebooks if you haven't used them before.

## Examples:

### Simple plume

A minimal example demonstrating a simple case to use pyplume with linear ambient temperature and salinity profiles and plot the results.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/alistaireverett/pyplume/master?filepath=examples%2Fsimple_plume.ipynb)

### Comparison to analytical solution

A slightly more detailed look at using pyplume and comparing the results to an
analytical solution.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/alistaireverett/pyplume/master?filepath=examples%2Fpyplume_vs_analytical.ipynb)

### Kronebreen Upwelling

The code used in Halbach et al. (20??) *[add link]* to calculate entrainment at Kronebreen throughout 2017 using pyplume. The model is forced with a year long timeseries of discharges calculated using a surface mass balance model ([Pramanik et al., 2018](https://doi.org/10.1017/jog.2018.80)) and ambient conditions taken from CTD casts close to the terminus.

***To be added soon***

More examples and documentation will be added shortly.

The code was developed by Alistair Everett with inspiration from Tom Cowton's iceplume Fortran model (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014JC010324)

