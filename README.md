# pyplume

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2198959.svg)](https://doi.org/10.5281/zenodo.2198959)

Pyplume calculates a numerical solution to the equations formulated for a point
source plume released at the base of a vertical wall. The code includes
capabilities for calculating the melt rate in the case where the vertical wall
is made of ice.

See the examples and docstrings for more details on how to use pyplume. The
links below will launch interactive jupyter notebooks with
[Binder](https://mybinder.readthedocs.io/en/latest/) so you can see how the code
works and experiment with the examples. [Here's a
link](https://medium.com/codingthesmartway-com-blog/getting-started-with-jupyter-notebook-for-python-4e7082bd5d46#f8b4)
explaining how to interact with jupyter notebooks if you haven't used them
before.

The first release (v0.1.0) contains the code used in Halbach et al. (submitted).
Example 3 is a simplified, interactive verion of the code used in the paper,
with a synthetic discharge timeseries. The full code used in the paper is
available in the examples folder as
[kronebreen_upwelling.py](https://github.com/alistaireverett/pyplume/blob/master/examples/kronebreen_upwelling.py).
The full datasets as used in the paper are available upon request from the authors.

## Examples:

### 1. Simple plume

A minimal example demonstrating how to use pyplume with linear ambient
temperature and salinity profiles and plot the results.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/alistaireverett/pyplume/master?filepath=examples%2Fsimple_plume.ipynb)

### 2. Comparison to analytical solution

A slightly more detailed look at using pyplume and comparing the results to an
analytical solution assuming uniform ambient temperature and salinity, with the
inlet Richardson number set to $\Gamma_0=1$.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/alistaireverett/pyplume/master?filepath=examples%2Fpyplume_vs_analytical.ipynb)

### 3. Discharge upwelling timeseries

An example demonstrating how pyplume can be used with a discharge timeseries to
calculate how the plume behaves on seasonal timescales. The example also
demonstrates a method which can be used to estimate plume driven upwelling, as
described in Halbach et al. (submitted).

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/alistaireverett/pyplume/master?filepath=examples%2Fupwelling_timeseries.ipynb)

More examples and documentation will be added in future.

## Acknowledgements

The code was developed by Alistair Everett with inspiration from Tom Cowton's
iceplume Fortran model [(Cowton et al.,
2015)](https://doi.org/10.1002/2014JC010324). Some of this development work
was done as part of the TIGRIF (RCN project number: 243808/E40) and TW-ICE
research projects at the Norwegian Polar Institute.

This code relies upon many previous studies which derived the solutions and
equations built into this model, references have been included throughout the code.
If you use pyplume, in addition to citing this code, please also consider citing
some of the important papers on which this code is based.

