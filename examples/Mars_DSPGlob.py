import numpy as np
import matplotlib.pyplot as plt
import pyshtools as pysh
from Py_Admittance import GlobalAdmitCorr
from Displacement_strain_planet import Thin_shell_matrix_nmax

#################################################################
# In this example, we show how to compute the observed and
# theoretical global admittance on Mars
# using our Displacement_strain_planet package (DSP).
# The DSP package allows to add more complexity to the admittance
# model and considers tangential loads. For more information, check
# https://github.com/AB-Ares/Displacement_strain_planet
#################################################################

lmax = 120  # Maximum spherical harmonic degree to perform all
# calculations
topo_lm = pysh.datasets.Mars.MOLA_shape(lmax=lmax).coeffs
pot_lm = pysh.datasets.Mars.GMM3(lmax=lmax)
R = topo_lm[0, 0, 0]  # Mean planetary radius planetary radius
pot_lm = pot_lm.change_ref(r0=R)  # Downward continue to Mean planetary radius

# Constants
G = pysh.constants.G.value  # Gravitational constant
gm = pot_lm.gm  # GM given in the gravity model file
mass = gm / G  # Mass of the planet
g0 = gm / R**2  # Mean gravitational attraction of the planet
rhobar = mass * 3.0 / 4.0 / np.pi / R**3

# Potential to free-air (mGal)
multipliers_l = np.arange(lmax + 1, dtype=float).reshape(1, -1, 1) + 1
FA_lm = pot_lm.coeffs * multipliers_l * g0 * 1e5

Tc = 50e3  # Crustal thickness
Te = 50e3  # Elastic thickness
rhoc = 2700  # Crustal density
rhol = 2700  # Surface load density
rhom = 3500  # Mantle density
E = 1e11  # Young's modulus
v = 0.25  # Poisson's ratio

args_param_m = (g0, R, Tc, Te, rhom, rhoc, rhol, lmax, E, v, mass)
G_lm = Thin_shell_matrix_nmax(
    *args_param_m,
    dc_lm=np.zeros((2, lmax + 1, lmax + 1)),
    drhom_lm=np.zeros((2, lmax + 1, lmax + 1)),
    H_lm=topo_lm
)[8]

# Geoid to free-air (mGal)
G_lm *= multipliers_l * g0 / R * 1e5

Admit_obs, Corr_obs, Admit_e_obs = GlobalAdmitCorr(topo_lm / 1e3, FA_lm)
Admit_Th, Corr_Th, Admit_e_Th = GlobalAdmitCorr(topo_lm / 1e3, G_lm)

f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
degrees = range(lmax + 1)
ax1.errorbar(
    degrees,
    Admit_obs,
    fmt="o",
    c="r",
    label="Observed admittance",
    yerr=Admit_e_obs,
)
ax1.plot(degrees, Admit_Th, "k", label="Theoretical admittance")
ax2.plot(degrees, Corr_Th, "k", label="Theoretical correlation")
ax2.plot(degrees, Corr_obs, "b", label="Observed correlation")
ax1.set_ylabel("Admittance, mGal/km")
ax1.set_xlabel("Degree")
ax2.set_ylabel("Correlation")
ax2.set_xlabel("Degree")
ax2.set_ylim(0, 1)
ax1.legend()
ax2.legend()
plt.show()
