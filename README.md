[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![DOI](https://zenodo.org/badge/781898403.svg)](https://zenodo.org/doi/10.5281/zenodo.10925820)

# Py_Admittance

Admittance calculations for geophysics and gravity analyses.

## Description

**Py_Admittance** provides functions and example scripts for calculating global/localized admittances of gravity and topography. These can be used to invert for parameters such as the elastic thickness of the lithosphere and density of the crust and surface.

The provided admittance model is from [Broquet & Wieczorek (2019)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019JE005959).

## Documentation

The full documentation can be found at https://ab-ares.github.io/Py_Admittance/.

## Dependencies

Some of these functions rely heavily on the [pyshtools](https://shtools.github.io/SHTOOLS/) package of [Wieczorek & Meschede (2018)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018GC007529) that is used to perform the spherical harmonic transforms, spatio-spectral localization, and finite-amplitude calculations.

An example script shows how one can use my [DSP (Displacement_strain_planet)]( https://github.com/AB-Ares/Displacement_strain_planet) package to compute theoretical admittances. This package can be used to add more complexity to the load model (see [this notebook](https://ab-ares.github.io/Displacement_strain_planet/notebooks/Run_demo.html)) when compared to the simple transfer function described in Py_Admittance. DSP also considers tangential loads. 

## Methods
`ForwardGravity`  Compute the theortical gravity field given a set of input parameters.

`TransferTGz`  Compute theoretical transfer functions for calculating the theortical admittance and correlation given a set of input parameters.

`LocalAdmitCorr`  Compute the localized admittance of two functions.

`GlobalAdmitCorr` Compute the global admittance of two functions.

## Example scripts
`Mars_AdmitGlob`  A script to compute the global theoretical and observed admittance and correlation on Mars given a set of input parameters. The theortical admittance is computed following [Broquet & Wieczorek (2019)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019JE005959).

`Mars_AdmitLoc`  A script to compute the localized theoretical and observed admittance and correlation on Mars given a set of input parameters. The theortical admittance is computed following [Broquet & Wieczorek (2019)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019JE005959) and localized to Olympus Mons.

`Mars_DSPGlob`  A script to compute the global theoretical and observed admittance and correlation on Mars given a set of input parameters. The theortical admittance is computed using the [DSP (Displacement_strain_planet)]( https://github.com/AB-Ares/Displacement_strain_planet) package (see [Broquet & Andrews-Hanna, 2023](https://doi.org/10.1038/s41550-022-01836-3)).

`Mars_DSPLoc`  A script to compute the localized theoretical and observed admittance and correlation on Mars given a set of input parameters. The theortical admittance is computed using the [DSP (Displacement_strain_planet)]( https://github.com/AB-Ares/Displacement_strain_planet) package (see [Broquet & Andrews-Hanna, 2023](https://doi.org/10.1038/s41550-022-01836-3)) and localized to Olympus Mons.

## How to install and run Py_Admittance
If you would like to modify the source code, download the Displacement_strain_planet repository and install using pip (or pip3 depending on your installation).
```bash
    git clone https://github.com/AB-Ares/Py_Admittance.git
    cd Py_Admittance/
    pip install .
```

## To run the example scripts
```bash
    cd examples
    python Mars_AdmitLoc.py 
    python Mars_DSPLoc.py 
```

## Author
[Adrien Broquet](https://ab-ares.github.io/website/) (adrien.broquet@dlr.de)

## Cite
You can cite the latest release of the package as:
Adrien Broquet. AB-Ares/Py_Admittance: 0.1.0 (Version 0.1.0). Zenodo. http://doi.org/10.5281/zenodo.10925820
