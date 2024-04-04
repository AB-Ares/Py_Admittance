"""
Py_Admittance
============================
Py_Admittance provides functions and an example scripts for calculating global/localized admittances of gravity and topography. These can be used to invert for parameters such as the elastic thickness of the lithosphere and density of the crust and surface.

   ForwardGravity
      Compute the theortical gravity field given a set of input
      parameters.

   TransferTGz
      Compute theoretical transfer functions for calculating
      the admittance and correlation given a set of input
      parameters.

   LocalAdmitCorr
      Compute the localized admittance of two functions.

   GlobalAdmitCorr
      Compute the global admittance of two functions.
"""

from ._version import get_versions

from .Th_Admittance import ForwardGravity
from .Th_Admittance import TransferTGz
from .Th_Admittance import LocalAdmitCorr
from .Th_Admittance import GlobalAdmitCorr

del Th_Admittance

__version__ = get_versions()["version"]
del get_versions

__author__ = "Adrien Broquet"

__all__ = [
    "ForwardGravity",
    "TransferTGz",
    "LocalAdmitCorr",
    "GlobalAdmitCorr",
]
