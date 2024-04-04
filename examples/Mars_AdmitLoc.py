import numpy as np
import matplotlib.pyplot as plt
import pyshtools as pysh
from Py_Admittance import ForwardGravity, LocalAdmitCorr

#################################################################
# In this example, we show how to compute the observed and
# theoretical localized admittance for Olympus Mons on Mars.
#################################################################

lmax = 120  # Maximum spherical harmonic degree to perform all
# calculations
topo_lm = pysh.datasets.Mars.MOLA_shape(lmax=lmax).coeffs
pot_lm = pysh.datasets.Mars.GMM3(lmax=lmax)
R = topo_lm[0, 0, 0]  # Mean planetary radius planetary radius
pot_lm = pot_lm.change_ref(r0=R)  # Downward continue to Mean
# planetary radius

# Constants
G = pysh.constants.G.value  # Gravitational constant
gm = pot_lm.gm  # GM given in the gravity model file
mass = gm / G  # Mass of the planet
g0 = gm / R**2  # Mean gravitational attraction of the planet

# Potential to free-air (mGal)
multipliers_l = np.arange(lmax + 1, dtype=float) + 1
FA_lm = pot_lm.coeffs * multipliers_l.reshape(1, -1, 1) * g0 * 1e5

Tc = 50e3  # Crustal thickness
Te = 32e3  # Elastic thickness
rhoc = 2900  # Crustal density
rhol = 3200  # Surface load density
rhom = 3500  # Mantle density

# Localization parameters
lat = 18.5  # Latitude
lon = 226.0  # Longitude
theta = 15  # Localization window angular radius
lwin = 17  # Localization window bandwidth

############################# Note #############################
# Mars has a strong degree-2 topography, which is due to the
# rotational flattening of the planet. This should not be
# considered as a load. With the option_deg2, you can say that
# only a fraction of the degree-2 topography should be considered
# as a load with option_deg2 = X. This effect is important when
# considered finite-amplitude correction. Check the doc for other
# options!
#################################################################

# Only 5% of the degree-2 topography is considered as a load
option_deg2 = 5.0 / 100.0
G_lm = ForwardGravity(
    topo_lm, G, mass, g0, Tc, Te, rhoc, rhom, lmax, rhol=rhol, option_deg2=option_deg2
)

Admit_obs, Corr_obs, Admit_e_obs = LocalAdmitCorr(
    topo_lm / 1e3, FA_lm, lat, lon, theta, lwin
)
Admit_Th, Corr_Th, Admit_e_Th = LocalAdmitCorr(
    topo_lm / 1e3, G_lm, lat, lon, theta, lwin
)

f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
degrees_lwin = range(lwin, lmax - lwin + 1)
ax1.errorbar(
    degrees_lwin,
    Admit_obs[lwin:],
    fmt="o",
    c="r",
    label="Observed admittance",
    yerr=Admit_e_obs[lwin:],
)
ax1.plot(degrees_lwin, Admit_Th[lwin:], "k", label="Theoretical admittance")
ax2.plot(degrees_lwin, Corr_Th[lwin:], "k", label="Theoretical correlation")
ax2.plot(degrees_lwin, Corr_obs[lwin:], "b", label="Observed correlation")
f.suptitle("Admittance and correlation at Olympus Mons")
ax1.set_ylabel("Admittance, mGal/km")
ax1.set_xlabel("Degree")
ax2.set_ylabel("Correlation")
ax2.set_xlabel("Degree")
ax2.set_ylim(0, 1)
ax1.legend()
ax2.legend()
plt.show()
