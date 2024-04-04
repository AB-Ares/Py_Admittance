"""
Compute theoretical transfer functions and gravity field given the input parameters. These
can then be used to compute the localized/global admittance and correlation. 
"""

import numpy as np
import pyshtools as pysh


def TransferTGz(
    g0,
    R,
    Tc,
    Te,
    rhol,
    rhoc,
    rhom,
    rhobar,
    lmax,
    ratio_L=0,
    alpha_L=1,
    depth_L=50e3,
    E=100e9,
    v=0.25,
):
    """
    Compute theoretical transfer function given the input parameters, including the admittance and correlation. For more information, see Broquet & Wieczorek (2019).

    Returns
    -------
    T_s : array, dimension (lmax+1)
                    Transfer function for flexure due to surface loads.
    Qw_l : array, dimension (lmax+1)
                    Transfer function for flexure due to surface/internal loads.
    Qw_lz : array, dimension (lmax+1)
                    Transfer function for flexure du to internal load.
    Trsf_s : array, dimension (lmax+1)
                    Transfer function for gravity from surface load.
    Trsf_L : array, dimension (lmax+1)
                    Transfer for gravity from internal load.
    Correlation : array, dimension (lmax+1)
                    Theoretical global correlation for out of
                    phase surface/internal loads (alpha_L != 1).

    Parameters
    ----------
    g0 : float
                    Gravitational attraction at the surface.
    R : float
                    Mean radius of the planet.
    Tc : float
                    Average crustal thickness.
    Te : float
                    Elastic thickness of the lithosphere.
    rhol : float
                    Density of the surface topography.
    rhoc : float
                    Density of the crust.
    rhom : float
                    Density of the mantle.
    lmax : int
                    Maximum spherical harmonic degree for calculations.
    ratio_L : float, optional, default = 0
                    Ratio of the internal / surface load.
    alpha_L : float, optional, default = 1
                    Phase relationship for the internal / surface load.
                    This parameter is experimental.
    depth_L : float, optional, default = 50e3
                    Depth of the internal load.
    E : float, optional, default = 100e9
                    Young's modulus.
    v : float, optional, default = 0.25
                    Poisson's ratio.
    """

    Re = R - Te / 2.0  # Reference radius for displacement equations
    D = E * Te**3 / (12.0 * (1.0 - v**2))  # Shell's rigidity.
    drho = rhom - rhoc
    drho2 = rhoc - rhol

    rhocd = rhoc / drho
    rhold = rhol / drho
    rhocm = rhoc / rhom
    rhodd = drho2 / drho
    rhodm = drho / rhom
    rhodc = drho / rhoc
    rholc = rhol / rhoc
    rholm = rhol / rhom
    rhomd = rhom / drho

    g0_rho_c_bar = 3.0 * g0 * rhoc / rhobar
    g0_rho_l_bar = 3.0 * g0 * rhol / rhobar

    RTCR = (R - Tc) / R
    RZR = (R - depth_L) / R
    RTCRZ = (R - Tc) / (R - depth_L)
    RZRTC = (R - depth_L) / (R - Tc)

    RE4 = Re**4
    psi = 12.0 * Re**2 / Te**2

    gmoho = g0 * (1.0 + (RTCR**3 - 1.0) * rhoc / rhobar) / RTCR**2
    if depth_L <= Tc:
        gz = g0 * (1.0 + (RZR**3 - 1.0) * rhoc / rhobar) / RZR**2
    else:
        gz = g0 * (1.0 + (RZR**3 - 1) * rhom / rhobar) / RZR**2
    gmohog0 = gmoho / g0
    gzg0 = gz / g0

    degrees = np.arange(lmax + 1, dtype=float)
    degrees_2 = degrees + 2.0
    degrees_1 = degrees + 1.0
    degrees_2_1 = 2.0 * degrees + 1
    rhobconst = 3.0 / (rhobar * degrees_2_1)
    rhobconst_c = rhoc * rhobconst
    rhobconst_d = drho * rhobconst
    rhobconst_dl = drho2 * rhobconst
    Lapla = -degrees * degrees_1  # Laplacian identity.
    epsi = np.where(
        degrees != 1,
        (
            -RE4
            / D
            * (
                (Lapla + 1.0 - v)
                / (
                    1e-12  # Avoid Division by zero error
                    + Lapla**3
                    + 4.0 * Lapla**2
                    + 4.0 * Lapla
                    + psi * (1.0 - v**2) * (Lapla + 2.0)
                )
            )
        ),
        np.inf,
    )

    Cs_bar = (1.0 - rhobconst_c - (rhobconst_d * RTCR**degrees)) / (
        gmohog0
        + rhodd
        - 1.0 / (epsi * g0 * drho)
        - rhobconst_c * (RTCR**degrees_2 + rhodd)
        - rhobconst_d * RTCR
        - rhobconst_dl * RTCR**degrees
    )

    if depth_L > Tc:
        # In the mantle
        phy_2 = (
            gzg0
            - (rhobconst_c * RZR**degrees_2)
            - (rhobconst_d * RZR * RZRTC**degrees_1)
        )
    else:
        # In the crust
        phy_2 = (
            gzg0 - (rhobconst_c * RZR**degrees_2) - (rhobconst_d * RZR * RTCRZ**degrees)
        )

    phy_3 = (
        rhocd
        + gmohog0
        - (1.0 / (epsi * g0 * drho))
        - rhobconst_c * (rhocd + RTCR**degrees_2 + rhodc * RTCR + RTCR**degrees)
    )

    Cz = (rhomd * phy_2) / phy_3

    phy_1 = 1.0 + rhold * Cs_bar
    Cs = rhomd * (Cs_bar / phy_1)
    Xs = (Cs * rholm) / (Cs * rholm - 1.0)

    Qs = (g0_rho_l_bar / degrees_2_1) * (
        1.0 - Cs_bar * rhodd - Cs_bar * RTCR**degrees_2
    )

    Qz = (g0_rho_c_bar / degrees_2_1) * (
        (RZR**degrees_2) - (Cz * rhocm) - (Cz * rhodm * (RTCR**degrees_2))
    )

    Czlm = Cz * rholm * (1.0 - Xs) * ratio_L
    Qzlc = Qz * rholc * (1.0 - Xs) * ratio_L

    T_s = 1.0 - Czlm

    Qw_l = (Xs - Czlm) / T_s

    Qw_lz = -Czlm / T_s

    Trsf_s = (Qs + Qzlc) / T_s

    Trsf_L = (
        (g0_rho_l_bar / degrees_2_1) * (RZR**degrees_2) * (1.0 - Xs) * ratio_L
    ) / T_s

    # From Wieczorek 2007, adapted for difference between rhol and rhoc.
    # From these, one can get the theoretical correlation for loads
    # that are not in phase. This parameter is experimental.
    if alpha_L != 1:
        Shh = 1.0 + Czlm**2 - alpha_L * 2.0 * Czlm
        Sgg = (degrees_1 / R) ** 2 * (Qs**2 + Qzlc**2 + alpha_L * (2.0 * Qs * Qzlc))
        Sgh = (degrees_1 / R) * (Qs - Czlm * Qzlc + alpha_L * (Qzlc - Qs * Czlm))
        Correlation = np.where((Shh * Sgg) != 0, Sgh / np.sqrt(Shh * Sgg), 0)
        Trsf_s = (Sgh / Shh) / (degrees_1 / R)
        alpha_correction = Trsf_s - ((Qs + Qzlc) / T_s)
    else:
        Correlation = np.ones_like(degrees)
        alpha_correction = 0.0

    return T_s, Qw_l, Qw_lz, Trsf_s, Trsf_L + alpha_correction, Correlation


def ForwardGravity(
    topo,
    G,
    mass,
    g0,
    Tc,
    Te,
    rhoc,
    rhom,
    lmax,
    ratio_L=0,
    alpha_L=1,
    depth_L=50e3,
    E=100e9,
    v=0.25,
    nmax=5,
    R_ref=None,
    rhol=None,
    lmaxgrid=None,
    option_deg1=None,
    option_deg2=None,
    return_corr=False,
):
    """
    Compute the theoretical gravity field given the input parameters.
    For more information, see Broquet & Wieczorek (2019).

    Returns
    -------
    grav_Th : array, dimension (2, lmax+1, lmax+1)
                    Theoretical global gravity field in mGal.
    Corr_Th : array, dimension (lmax+1)
                    If return_corr is True, returns the theoretical
                    correlation for surface/internal loads that are
                    out of phase (alpha_L != 1).

    Parameters
    ----------
    topo : array, dimension (2,lmax+1,lmax+1)
                    Array with the spherical harmonic coefficients
                    of the surface topography.
    G : float
                    Gravitational constant.
    mass : float
                    Mass of the planet.
    g0 : float
                    Gravitational attraction at the surface.
    Tc : float
                    Average crustal thickness.
    Te : float
                    Elastic thickness of the lithosphere.
    rhoc : float
                    Density of the crust.
    rhom : float
                    Density of the mantle.
    lmax : int
                    Maximum spherical harmonic degree for calculations.
    ratio_L : float, optional, default = 0
                    Ratio of the internal / surface load.
    alpha_L : float, optional, default = 1
                    Phase relationship for the internal / surface load.
                    This parameter is experimental.
    depth_L : float, optional, default = 50e3
                    Depth of the internal load.
    E : float, optional, default = 100e9
                    Young's modulus.
    v : float, optional, default = 0.25
                    Poisson's ratio.
    nmax : int, optional, default = 5
                    Order of the finite-amplitude correction.
    R_ref : float, optional, default = None
                    Reference radius for gravity field calculations.
                    If None, this parameter is set to the mean radius of
                    the topography file.
    rhol : float, optional, default = None
                    Density of the surface topography. If None, this
                    parameter is set to rhoc.
    lmaxgrid : int, optional, default = None
                    Resolution of the input grid for the finite-amplitude correction
                    routines. If None, this parameter is set to 3*lmax.
                    For accurate results, this parameter should be about
                    3 times lmax, though this should be verified for each application.
                    Lowering this parameter significantly increases speed.
    option_deg1 : string, optional, default = None
                    How to treat degree-1 displacement. If set to "Zero",
                    the degree-1 displacement is zeroed out. If set to "Airy",
                    the degree-1 topography is assumed to be Airy compensated.
                    If anything else, no special treatment is applied to
                    degree-1.
    option_deg2 : string/float, optional, default = None
                    How to treat degree-2 topography. If set to "Zero",
                    the C20 topography is zeroed out. If set to "Flat",
                    the C20 topography is not considered as a load, but
                    is added back for finite-amplitude calculations.
                    If set to a float, only option_deg2 * topography is used as a load.
                    If anything else, no special treatment is applied to degree-2.
    return_corr : string, optional, default = False
                    If set to True, return the theoretical global correlation.
    """

    if lmaxgrid is None:
        lmaxgrid = 3 * lmax
    if rhol is None:
        rhol = rhoc

    kw_exp = {"sampling": 2, "lmax": lmaxgrid, "lmax_calc": lmax}
    topo_clm = topo.copy()
    R = topo_clm[0, 0, 0]
    if R_ref is None:
        R_ref = R

    if option_deg2 == "Zero":
        topo_clm[0, 2, 0] = 0.0
        C20 = 0.0
    elif option_deg2 == "Flat":
        C20 = topo_clm[0, 2, 0]
        topo_clm[0, 2, 0] = 0.0
    elif isinstance(option_deg2, float):
        C20 = (1.0 - option_deg2) * topo_clm[0, 2, 0]
        topo_clm[0, 2, 0] *= option_deg2
    else:
        C20 = 0.0

    rhobar = mass * 3.0 / 4.0 / np.pi / R**3

    # Compute transfer functions
    T_s, Qw_l, Qw_lz, Trsf_s, Trsf_L, Corr_Th = TransferTGz(
        g0,
        R,
        Tc,
        Te,
        rhol,
        rhoc,
        rhom,
        rhobar,
        lmax,
        ratio_L=ratio_L,
        alpha_L=alpha_L,
        depth_L=depth_L,
        E=E,
        v=v,
    )
    Trsf_L *= 1.0e5
    Trsf_s *= 1.0e5
    multipliers_l = (np.arange(lmax + 1, dtype=float) + 1).reshape(1, -1, 1)

    if nmax == 1:
        grav_Th = (
            topo_clm
            * Trsf_s.reshape(1, -1, 1)
            * (multipliers_l / R_ref)
            * (R / R_ref) ** (multipliers_l - 1)
        )
        grav_Th[0, 0, 0] = 0.0
        return grav_Th

    W_lm = topo_clm * Qw_l.reshape(1, -1, 1)  # Both surface/internal loads
    W_lmz = topo_clm * Qw_lz.reshape(1, -1, 1)  # Just internal load
    topo_s = topo_clm / T_s.reshape(1, -1, 1)  # No internal load

    # Define degree-0 relief
    W_lm[0, 0, 0] = R
    W_lmz[0, 0, 0] = R
    topo_s[0, 0, 0] = R

    # Add back degree-2 relief for finite-amplitude calculations
    W_lm[0, 2, 0] += C20
    W_lmz[0, 2, 0] += C20
    topo_s[0, 2, 0] += C20

    if ratio_L != 0:
        # Deflection of the surface due to the internal load (not considered in grid_topo_s)
        grid_Wz = pysh.expand.MakeGridDH(W_lmz, **kw_exp)
        B_wl, r0_deflecz = pysh.gravmag.CilmPlusDH(grid_Wz, nmax, mass, rhol, lmax=lmax)
        B_load = (
            topo_clm
            * Trsf_L.reshape(1, -1, 1)
            * (multipliers_l / R_ref)
            * (R / R_ref) ** multipliers_l
        )
        B_load[0, 0, 0] = 0.0
    else:
        B_wl = 0.0
        B_load = 0.0

    # Deg-1 Airy at the surface
    if option_deg1 == "Airy":
        if rhol == rhoc:
            W_lm[:, 1, :2] = 0.0
        else:
            W_lm[:, 1, :2] = -((topo_clm[:, 1, : 1 + 1] * rhol) / (rhoc - rhol))
    elif option_deg1 == "Zero":
        W_lm[:, 1, :2] = 0.0

    grid_W = pysh.expand.MakeGridDH(W_lm, **kw_exp)
    grid_topo_s = pysh.expand.MakeGridDH(topo_s, **kw_exp)

    # Gravity contribution of topography caused by the surface load
    B_tl, r0 = pysh.gravmag.CilmPlusDH(grid_topo_s, nmax, mass, rhol, lmax=lmax)

    if rhoc != rhol:
        # Gravity contribution Topo/Crust relief due to the surface/internal load
        B_wcl, r0_deflec = pysh.gravmag.CilmPlusDH(
            grid_W, nmax, mass, rhoc - rhol, lmax=lmax
        )
    else:
        B_wcl = 0.0
        r0_deflec = 0.0

    # Deg-1 Airy at the moho
    if option_deg1 == "Airy":
        W_lm[:, 1, :2] = -((topo_clm[:, 1, :2] * rhoc) / (rhom - rhoc)) * (
            g0 / (g0 * (1.0 - Tc / R))
        )
        grid_W = pysh.expand.MakeGridDH(W_lm, **kw_exp)
    elif option_deg1 == "Zero":
        W_lm[:, 1, :2] = 0.0

    # Gravity contribution of Crust/Mantle relief due to the surface/internal load
    B_wmc, r0_deflec_Tc = pysh.gravmag.CilmPlusDH(
        grid_W - Tc, nmax, mass, rhom - rhoc, lmax=lmax
    )

    # Potential to Free-air
    B_tl *= multipliers_l * (r0 / R_ref) ** (multipliers_l - 1)
    if ratio_L != 0:
        B_wl *= multipliers_l * (r0_deflecz / R_ref) ** (multipliers_l - 1)
    B_wcl *= multipliers_l * (r0_deflec / R_ref) ** (multipliers_l - 1)
    B_wmc *= multipliers_l * (r0_deflec_Tc / R_ref) ** (multipliers_l - 1)

    grav_Th = (B_tl + B_wcl + B_wmc + B_wl) * 1.0e5 * (G * mass / R_ref**2) + B_load
    grav_Th[0, 0, 0] = 0.0

    if not return_corr:
        return grav_Th

    return grav_Th, Corr_Th


def LocalAdmitCorr(topo, grav, lat, lon, theta, lwin, lmax=None, quiet=True):
    """
    Compute the localized admittance and correlation functions from the
    input gravity and topography. For more information, see Broquet &
    Wieczorek (2019).

    Returns
    -------
    Admittance : array, dimension (lmax-lwin+1)
                    Localized admittance function in mGal/km.
    Correlation : array, dimension (lmax-lwin+1)
                    Localized correlation function.
    Admit_error : array, dimension (lmax-lwin+1)
                    Localized admittance uncertainty in mGal/km.

    Parameters
    ----------
    topo : dimension (2, lmax+1, lmax+1)
                    Spherical harmonic coefficients of the topography (km).
    grav : dimension (2, lmax+1, lmax+1)
                    Spherical harmonic coefficients of the gravity field (mGal).
    lat : float
                    Central latitude (°) of the localization window.
    lon : float
                    Central longitude (°) of the localization window.
    theta : float
                    Angular radius (°) of the localization window.
    lwin : int
                    Bandwidth of the localization window.
    lmax : int, optional, default = None
                    Maximum degree at which the admittance and correlation
                    are computed. If None, lmax = min(lmax_topo, lmax_grav).
                    lmax must be <= min(lmax_topo, lmax_grav).
    quiet : string, optional, default = True
                    If False, the function will provide information regarding the
                    spatio-spectral concentration of the localization window.
    """
    if not quiet:
        capwin = pysh.SHWindow.from_cap(theta=theta, lwin=lwin)
        print("Best concentrated window energy %s" % (capwin.number_concentrated(0.99)))
        print(
            "Number of optimally (>99%%) concentrated windows %s"
            % (capwin.number_concentrated(0.99))
        )

    if lmax is None:
        lmax = np.min([np.shape(topo)[1] - 1, np.shape(grav)[1] - 1])

    taps, _, taper_order = pysh.spectralanalysis.SHReturnTapers(
        theta * np.pi / 180.0, lwin
    )
    dj_lwin = pysh.rotate.djpi2(lwin)
    SH_tapers = pysh.SHCoeffs.from_zeros(lwin)
    SH_tapers.coeffs[:, : lwin + 1, taper_order[0]] = taps[: lwin + 1, 0]
    SHtaperRot = pysh.rotate.SHRotateRealCoef(
        SH_tapers.coeffs,
        np.asarray([0.0, -(90.0 - lat) * np.pi / 180.0, -lon * np.pi / 180.0]),
        dj_lwin,
    )
    grid_taper = pysh.expand.MakeGridDH(SHtaperRot, lmax=lmax, extend=0)
    grav_grid = pysh.expand.MakeGridDH(grav, lmax=lmax, extend=0)
    topo_grid = pysh.expand.MakeGridDH(topo, lmax=lmax, extend=0)
    grav_clm_loc = pysh.expand.SHExpandDH(grid_taper * grav_grid, lmax_calc=lmax - lwin)
    topo_clm_loc = pysh.expand.SHExpandDH(grid_taper * topo_grid, lmax_calc=lmax - lwin)
    admit_loc, error_loc, corr_loc = pysh.spectralanalysis.SHAdmitCorr(
        grav_clm_loc, topo_clm_loc
    )

    return admit_loc, corr_loc, error_loc


def GlobalAdmitCorr(topo, grav, lmax=None):
    """
    Compute the global admittance and correlation functions from the
    input gravity and topography. For more information, see Broquet &
    Wieczorek (2019).

    Returns
    -------
    Admittance : array, dimension (lmax+1)
                    Global admittance function in mGal/km.
    Correlation : array, dimension (lmax+1)
                    Global correlation function.
    Admit_error : array, dimension (lmax+1)
                    Global admittance uncertainty in mGal/km.

    Parameters
    ----------
    topo : dimension (2, lmax+1, lmax+1)
                    Spherical harmonic coefficients of the topography (km).
    grav : dimension (2, lmax+1, lmax+1)
                    Spherical harmonic coefficients of the gravity field (mGal).
    lmax : int, optional, default = None
                    Maximum degree at which the admittance and correlation
                    are computed. If None, lmax = min(lmax_topo, lmax_grav).
                    lmax must be <= min(lmax_topo, lmax_grav).
    """

    admit, error, corr = pysh.spectralanalysis.SHAdmitCorr(grav, topo, lmax=lmax)

    return admit, corr, error
