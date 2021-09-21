import math

import numpy
from parcels import ParcelsRandom
import settings
import scipy.optimize
import numpy as np
import os
from copy import deepcopy
from utils.file_utils import print_statement


def anti_beach_nudging(particle, fieldset, time):
    """
    If a particle is within 0.5km of the nearest coastline (determined by sampling
    the distance2shore field), then it gets nudged away back out to sea. I have
    fields for border currents, so that the particle gets nudged in the right
    direction with a speed of 1 - 1.414 m s^{-1}.

    With dt=10 minutes a particle gets displaced by 600 - 848 m back out to sea.
    """
    d1 = particle.depth
    if fieldset.distance2shore[time, d1, particle.lat, particle.lon] < 0.5:
        borUab, borVab = fieldset.borU[time, d1, particle.lat, particle.lon], fieldset.borV[
            time, d1, particle.lat, particle.lon]
        particle.lon -= borUab * particle.dt
        particle.lat -= borVab * particle.dt


def floating_advection_rk4(particle, fieldset, time):
    """Advection of particles using fourth-order Runge-Kutta integration.

    Function needs to be converted to Kernel object before execution

    A particle only moves if it has not beached (rather obviously)
    """
    if particle.beach == 0:
        # particle.distance = fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
        d2 = particle.depth
        if particle.lon > 180:
            particle.lon -= 360
        if particle.lon < -180:
            particle.lon += 360
        (u1, v1) = fieldset.UV[time, d2, particle.lat, particle.lon]
        (uS1, vS1) = fieldset.Ust[time, d2, particle.lat, particle.lon], fieldset.Vst[
            time, d2, particle.lat, particle.lon]
        # lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
        lon1, lat1 = (particle.lon + (u1 + uS1) * .5 * particle.dt, particle.lat + (v1 + vS1) * .5 * particle.dt)

        if lon1 > 180:
            lon1 -= 360
        if lon1 < -180:
            lon1 += 360
        (u2, v2) = fieldset.UV[time + .5 * particle.dt, d2, lat1, lon1]
        (uS2, vS2) = fieldset.Ust[time + .5 * particle.dt, d2, lat1, lon1], fieldset.Vst[
            time + .5 * particle.dt, d2, lat1, lon1]
        # lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
        lon2, lat2 = (particle.lon + (u2 + uS2) * .5 * particle.dt, particle.lat + (v2 + vS2) * .5 * particle.dt)

        if lon2 > 180:
            lon2 -= 360
        if lon2 < -180:
            lon2 += 360
        (u3, v3) = fieldset.UV[time + .5 * particle.dt, d2, lat2, lon2]
        (uS3, vS3) = fieldset.Ust[time + .5 * particle.dt, d2, lat2, lon2], fieldset.Vst[
            time + .5 * particle.dt, d2, lat2, lon2]
        # lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
        lon3, lat3 = (particle.lon + (u3 + uS3) * particle.dt, particle.lat + (v3 + vS3) * particle.dt)

        if lon3 > 180:
            lon3 -= 360
        if lon3 < -180:
            lon3 += 360
        (u4, v4) = fieldset.UV[time + particle.dt, d2, lat3, lon3]
        (uS4, vS4) = fieldset.Ust[time + particle.dt, d2, lat3, lon3], fieldset.Vst[time + particle.dt, d2, lat3, lon3]

        # particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lon += ((u1 + uS1) + 2 * (u2 + uS2) + 2 * (u3 + uS3) + (u4 + uS4)) / 6. * particle.dt
        if particle.lon > 180:
            particle.lon -= 360
        if particle.lon < -180:
            particle.lon += 360
        # particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        particle.lat += ((v1 + vS1) + 2 * (v2 + vS2) + 2 * (v3 + vS3) + (v4 + vS4)) / 6. * particle.dt


def floating_2d_brownian_motion(particle, fieldset, time):
    """Kernel for simple Brownian particle diffusion in zonal and meridional direction.
    Assumes that fieldset has fields Kh_zonal and Kh_meridional
    we don't want particles to jump on land and thereby beach"""
    if particle.beach == 0:
        # Wiener increment with zero mean and std of sqrt(dt)
        dWx = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
        dWy = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)

        bx = math.sqrt(2 * fieldset.Kh_zonal[time, particle.depth, particle.lat, particle.lon])
        by = math.sqrt(2 * fieldset.Kh_meridional[time, particle.depth, particle.lat, particle.lon])

        particle.lon += bx * dWx
        particle.lat += by * dWy


def floating_AdvectionRK4DiffusionEM_stokes_depth(particle, fieldset, time):
    """Kernel for 2D advection-diffusion,  with advection solved
    using fourth order Runge-Kutta (RK4) and diffusion using the
    Euler-Maruyama scheme (EM). Using the RK4 scheme for diffusion is
    only advantageous in areas where the contribution from diffusion
    is negligible.
    Assumes that fieldset has fields `Kh_zonal` and `Kh_meridional`
    and variable `fieldset.dres`, setting the resolution for the central difference
    gradient approximation. This should be at least an order of magnitude
    less than the typical grid resolution.
    The Wiener increment `dW` should be normally distributed with zero
    mean and a standard deviation of sqrt(dt). Instead, here a uniform
    distribution with the same mean and std is used for efficiency and
    to keep random increments bounded. This substitution is valid for
    the convergence of particle distributions. If convergence of
    individual particle paths is required, use normally distributed
    random increments instead. See Gräwe et al (2012)
    doi.org/10.1007/s10236-012-0523-y for more information.

    Comment 16/11/2020: Changes I have made from the original parcels kernel is that the advection and diffusion only
    occurs when a particle is afloat. We also add in stokes drift separately, and include depth dependence for the
    magnitude of the Stokes drift.
    depth dependence: Breivik et al. (2016) https://doi.org/10.1016/j.ocemod.2016.01.005
    """
    if particle.beach == 0:
        t, d, la, lo, dt = time, particle.depth, particle.lat, particle.lon, particle.dt
        # Stokes depth dependence term
        w_p = (2 * math.pi / fieldset.WP[t, d, la, lo])  # Peak period

        k_p = w_p ** 2 / 9.81  # peak wave number
        z_correc = max(d - fieldset.SURF_Z, 0)  # depth correction
        st_z = min(max(math.exp(-2 * k_p * z_correc) - math.sqrt(2 * math.pi * k_p * z_correc) * math.erfc(2 * k_p * z_correc), 0), 1)


        # RK4 terms
        (u1, v1) = fieldset.UV[t, d, la, lo]
        (uS1, vS1) = fieldset.Ust[t, d, la, lo], fieldset.Vst[t, d, la, lo]
        lon1, lat1 = (lo + (u1 + uS1 * st_z) * .5 * particle.dt, la + (v1 + vS1 * st_z) * .5 * dt)

        (u2, v2) = fieldset.UV[t + .5 * dt, d, lat1, lon1]
        (uS2, vS2) = fieldset.Ust[t + .5 * dt, d, lat1, lon1], fieldset.Vst[t + .5 * dt, d, lat1, lon1]
        lon2, lat2 = (lo + (u2 + uS2 * st_z) * .5 * dt, la + (v2 + vS2 * st_z) * .5 * dt)

        (u3, v3) = fieldset.UV[t + .5 * dt, d, lat2, lon2]
        (uS3, vS3) = fieldset.Ust[t + .5 * dt, d, lat2, lon2], fieldset.Vst[t + .5 * dt, d, lat2, lon2]
        lon3, lat3 = (lo + (u3 + uS3 * st_z) * dt, la + (v3 + vS3 * st_z) * dt)

        (u4, v4) = fieldset.UV[t + dt, d, lat3, lon3]
        (uS4, vS4) = fieldset.Ust[t + dt, d, lat3, lon3], fieldset.Vst[t + dt, d, lat3, lon3]

        # Wiener increment with zero mean and std of sqrt(dt)
        dWx = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
        dWy = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)

        Kxp1 = fieldset.Kh_zonal[t, d, la, lo + fieldset.dres]
        Kxm1 = fieldset.Kh_zonal[t, d, la, lo - fieldset.dres]
        dKdx = (Kxp1 - Kxm1) / (2 * fieldset.dres)
        bx = math.sqrt(2 * fieldset.Kh_zonal[t, d, la, lo])

        Kyp1 = fieldset.Kh_meridional[t, d, la + fieldset.dres, lo]
        Kym1 = fieldset.Kh_meridional[t, d, la - fieldset.dres, lo]
        dKdy = (Kyp1 - Kym1) / (2 * fieldset.dres)
        by = math.sqrt(2 * fieldset.Kh_meridional[t, d, la, lo])

        # Particle positions are updated only after evaluating all terms.
        particle.lon += (((u1 + uS1 * st_z) + 2 * (u2 + uS2 * st_z) + 2 * (u3 + uS3 * st_z) + (u4 + uS4 * st_z))
                         / 6. + dKdx) * particle.dt + bx * dWx
        particle.lat += (((v1 + vS1 * st_z) + 2 * (v2 + vS2 * st_z) + 2 * (v3 + vS3 * st_z) + (v4 + vS4 * st_z))
                         / 6. + dKdy) * particle.dt + by * dWy


def initial_input(particle, fieldset, time):
    """
    Since we have many instances that particles start at the very same position,
    when a particle is first added to the simulation it will get a random kick
    that moves it slightly away from the initial position, and so with multiple
    particles we will see a sort of star pattern around the central position.
    However, the particle shouldn't be put on land though, so a particle will
    have 100000 attempts to be placed in the simulation within a cell that is not
    land. If it doesn't work after 100000 attempts, then the particle just ends up
    starting from the unchanged position. We also set it so that the initial
    position of the particle is always closer to land than the original position

    Note: Tests show that at most we get at most around 7 or 8 particles per
    release getting placed immediately on land. this varies a bit

    """
    if particle.age == 0:
        check = 0
        distCur = fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
        d_lon = fieldset.dlon[time, particle.depth, particle.lat, particle.lon]
        d_lat = fieldset.dlat[time, particle.depth, particle.lat, particle.lon]
        while check < 100000:
            potential_lat = particle.lat + ParcelsRandom.uniform(-d_lat, d_lat)
            potential_lon = particle.lon + ParcelsRandom.uniform(-d_lon, d_lon)
            potential_land = math.floor(fieldset.landID[time, particle.depth, particle.lat, particle.lon])
            potential_dist = fieldset.distance2shore[time, particle.depth, potential_lat, potential_lon]
            if potential_land == 0 and potential_dist <= distCur:
                check += 100001
                particle.lat = potential_lat
                particle.lon = potential_lon
            check += 1


def delete_particle(particle, fieldset, time):
    # This delete particle format from Philippe Delandmeter
    # https://github.com/OceanParcels/Parcelsv2.0PaperNorthSeaScripts/blob/master/northsea_mp_kernels.py
    print("Particle [%d] lost !! (%g %g %g %g)" % (
    particle.id, particle.lon, particle.lat, particle.depth, particle.time))
    particle.delete()


def get_kinematic_viscosity(particle, fieldset, time):
    # Using equations 25 - 29 from Kooi et al. 2017
    # Assuming that the fieldset has a salinity field called abs_salinity and a temperature field called
    # cons_temperature.
    # Salinity units: PSU or g/kg
    # Temperature units: temperature
    # Salinity and Temperature at the particle position, where salinity is converted from g/kg -> kg/kg
    Sz = fieldset.abs_salinity[time, particle.depth, particle.lat, particle.lon] / 1000
    Tz = fieldset.cons_temperature[time, particle.depth, particle.lat, particle.lon]
    # The constants A and B
    A = 1.541 + 1.998 * 10 ** -2 * Tz - 9.52 * 10 ** -5 * math.pow(Tz, 2)
    B = 7.974 - 7.561 * 10 ** -2 * Tz + 4.724 * 10 ** -4 * math.pow(Tz, 2)
    # Calculating the water dynamic viscosity
    mu_wz = 4.2844 * 10 ** -5 + math.pow(0.156 * math.pow(Tz + 64.993, 2) - 91.296, -1)
    # Calculating the sea water kinematic viscosity
    particle.kinematic_viscosity = mu_wz * (1 + A * Sz + B * math.pow(Sz, 2)) / particle.density


def kooi_rising_velocity(particle, fieldset, time):
    """
    Kernel to compute the vertical velocity (Vs) of particles due to their different sizes and densities
    This is very heavily based on the Kernel "Kooi_no_biofouling" written by Delphine Lobelle
    https://github.com/dlobelle/TOPIOS/blob/master/scripts/Kooi%2BNEMO_3D_nobiofoul.py
    """

    # ------ Profiles from MEDUSA or Kooi theoretical profiles -----
    z = particle.depth  # [m]
    kin_visc = particle.kinematic_viscosity  # kinematic viscosity[m2 s-1]
    rho_sw = particle.density  # seawater rho[kg m-3]
    rise = particle.rise_velocity  # vertical velocity[m s-1]

    # ------ Constants -----
    g = 7.32e10 / (86400. ** 2.)  # gravitational acceleration (m d-2), now [s-2]

    # ------ Diffusivity -----
    r_tot = particle.size  # total radius [m]
    rho_tot = (particle.size ** 3. * particle.rho_plastic) / (particle.size) ** 3.  # total rho [kg m-3]

    dn = 2. * r_tot  # equivalent spherical diameter [m]
    delta_rho = (
                        rho_tot - rho_sw) / rho_sw  # normalised difference in rho between total plastic+bf and seawater[-]
    dstar = ((rho_tot - rho_sw) * g * dn ** 3.) / (rho_sw * kin_visc ** 2.)  # dimensional diameter[-]

    # Getting the dimensionless settling velocity
    if dstar > 5e9:
        w = 1000.
    elif dstar < 0.05:
        w = (dstar ** 2.) * 1.71E-4
    else:
        w = 10. ** (-3.76715 + (1.92944 * math.log10(dstar)) - (0.09815 * math.log10(dstar) ** 2.) - (
                0.00575 * math.log10(dstar) ** 3.) + (0.00056 * math.log10(dstar) ** 4.))
    # ------ Settling of particle -----

    if delta_rho > 0:  # sinks
        vs = (g * kin_visc * w * delta_rho) ** (1. / 3.)
    else:  # rises
        a_del_rho = delta_rho * -1.
        vs = -1. * (g * kin_visc * w * a_del_rho) ** (1. / 3.)  # m s-1

    particle.rise_velocity = vs


def PolyTEOS10_bsq(particle, fieldset, time):
    '''
    calculates rho based on the polyTEOS10-bsq algorithm from Appendix A.2 of
    https://www.sciencedirect.com/science/article/pii/S1463500315000566
    requires fieldset.abs_salinity and fieldset.cons_temperature Fields in the fieldset
    and a particle.density Variable in the ParticleSet
    References:
    Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate
    polynomial expressions for the rho and specific volume of
    seawater using the TEOS-10 standard. Ocean Modelling.
    McDougall, T. J., D. R. Jackett, D. G. Wright and R. Feistel, 2003:
    Accurate and computationally efficient algorithms for potential
    temperature and rho of seawater.  Journal of Atmospheric and
    Oceanic Technology, 20, 730-741.
    '''

    Z = - particle.depth  # note: use negative depths!
    SA = fieldset.abs_salinity[time, particle.depth, particle.lat, particle.lon]
    CT = fieldset.cons_temperature[time, particle.depth, particle.lat, particle.lon]
    SA_SURF = fieldset.abs_salinity[time, fieldset.SURF_Z, particle.lat, particle.lon]
    CT_SURF = fieldset.cons_temperature[time, fieldset.SURF_Z, particle.lat, particle.lon]

    SAu = 40 * 35.16504 / 35
    CTu = 40
    Zu = 1e4
    deltaS = 32
    R000 = 8.0189615746e+02
    R100 = 8.6672408165e+02
    R200 = -1.7864682637e+03
    R300 = 2.0375295546e+03
    R400 = -1.2849161071e+03
    R500 = 4.3227585684e+02
    R600 = -6.0579916612e+01
    R010 = 2.6010145068e+01
    R110 = -6.5281885265e+01
    R210 = 8.1770425108e+01
    R310 = -5.6888046321e+01
    R410 = 1.7681814114e+01
    R510 = -1.9193502195e+00
    R020 = -3.7074170417e+01
    R120 = 6.1548258127e+01
    R220 = -6.0362551501e+01
    R320 = 2.9130021253e+01
    R420 = -5.4723692739e+00
    R030 = 2.1661789529e+01
    R130 = -3.3449108469e+01
    R230 = 1.9717078466e+01
    R330 = -3.1742946532e+00
    R040 = -8.3627885467e+00
    R140 = 1.1311538584e+01
    R240 = -5.3563304045e+00
    R050 = 5.4048723791e-01
    R150 = 4.8169980163e-01
    R060 = -1.9083568888e-01
    R001 = 1.9681925209e+01
    R101 = -4.2549998214e+01
    R201 = 5.0774768218e+01
    R301 = -3.0938076334e+01
    R401 = 6.6051753097e+00
    R011 = -1.3336301113e+01
    R111 = -4.4870114575e+00
    R211 = 5.0042598061e+00
    R311 = -6.5399043664e-01
    R021 = 6.7080479603e+00
    R121 = 3.5063081279e+00
    R221 = -1.8795372996e+00
    R031 = -2.4649669534e+00
    R131 = -5.5077101279e-01
    R041 = 5.5927935970e-01
    R002 = 2.0660924175e+00
    R102 = -4.9527603989e+00
    R202 = 2.5019633244e+00
    R012 = 2.0564311499e+00
    R112 = -2.1311365518e-01
    R022 = -1.2419983026e+00
    R003 = -2.3342758797e-02
    R103 = -1.8507636718e-02
    R013 = 3.7969820455e-01
    ss = math.sqrt((SA + deltaS) / SAu)
    tt = CT / CTu
    zz = -Z / Zu
    rz3 = R013 * tt + R103 * ss + R003
    rz2 = (R022 * tt + R112 * ss + R012) * tt + (R202 * ss + R102) * ss + R002
    rz1 = (((R041 * tt + R131 * ss + R031) * tt + (R221 * ss + R121) * ss + R021) * tt + (
                (R311 * ss + R211) * ss + R111) * ss + R011) * tt + (
                      ((R401 * ss + R301) * ss + R201) * ss + R101) * ss + R001
    rz0 = (((((R060 * tt + R150 * ss + R050) * tt + (R240 * ss + R140) * ss + R040) * tt + (
                (R330 * ss + R230) * ss + R130) * ss + R030) * tt + (
                        ((R420 * ss + R320) * ss + R220) * ss + R120) * ss + R020) * tt + (
                       (((R510 * ss + R410) * ss + R310) * ss + R210) * ss + R110) * ss + R010) * tt + (
                      ((((R600 * ss + R500) * ss + R400) * ss + R300) * ss + R200) * ss + R100) * ss + R000
    particle.density = ((rz3 * zz + rz2) * zz + rz1) * zz + rz0
    ss = math.sqrt((SA_SURF + deltaS) / SAu)
    tt = CT_SURF / CTu
    rz0 = (((((R060 * tt + R150 * ss + R050) * tt + (R240 * ss + R140) * ss + R040) * tt + (
                (R330 * ss + R230) * ss + R130) * ss + R030) * tt + (
                        ((R420 * ss + R320) * ss + R220) * ss + R120) * ss + R020) * tt + (
                       (((R510 * ss + R410) * ss + R310) * ss + R210) * ss + R110) * ss + R010) * tt + (
                      ((((R600 * ss + R500) * ss + R400) * ss + R300) * ss + R200) * ss + R100) * ss + R000
    particle.surface_density = rz0


def KPP_wind_mixing(particle, fieldset, time):
    """
    Markov-0 implementation of wind mixing based on the KPP wind mixing parametrization. For more exact details, look
    at the full wind mixing code at https://github.com/VictorOnink/Wind-Mixing-Diffusion/ or at the paper draft
    """
    if particle.beach == 0:
        # Loading the mixed layer depth
        mld = fieldset.MLD[time, fieldset.SURF_Z, particle.lat, particle.lon]

        # Below the MLD there is no wind-driven turbulent diffusion according to KPP theory
        if particle.depth > mld:
            Kz = 0
            dKz = 0
        # Within the MLD we compute the vertical diffusion according to Boufadel et al. (2020)
        else:
            # Wind speed
            if particle.lon < 0:
                temp_lon = particle.lon + 360.
            else:
                temp_lon = particle.lon
            w_10 = math.sqrt(fieldset.u10[time, fieldset.SURF_Z, particle.lat, temp_lon] ** 2 + \
                             fieldset.v10[time, fieldset.SURF_Z, particle.lat, temp_lon] ** 2)
            # Drag coefficient
            C_D = min(max(1.2E-3, 1.0E-3 * (0.49 + 0.065 * w_10)), 2.12E-3)
            # wind stress
            tau = C_D * fieldset.RHO_A * w_10 ** 2
            # Frictional velocity of water at the ocean surface
            U_W = math.sqrt(tau / particle.surface_density)
            # Surface roughness z0 following Zhao & Li (2019)
            z0 = 3.5153e-5 * fieldset.BETA ** (-0.42) * w_10 ** 2 / fieldset.G
            # The corrected particle depth, since the depth is not always zero for the surface circulation data
            z_correct = particle.depth - fieldset.SURF_Z
            # The diffusion gradient at particle.depth
            alpha = (fieldset.VK * U_W) / (fieldset.PHI * mld ** 2)
            dKz = alpha * (mld - z_correct) * (mld - 3 * z_correct - 2 * z0)
            # The KPP profile vertical diffusion, at a depth corrected for the vertical gradient in Kz, and including the
            # bulk diffusivity
            alpha = (fieldset.VK * U_W) / fieldset.PHI
            Kz = alpha * (z_correct + z0) * math.pow(1 - z_correct / mld, 2)
        # The Markov-0 vertical transport from Grawe et al. (2012)
        gradient = dKz * particle.dt
        R = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3) * math.sqrt(2 * Kz)
        rise = particle.rise_velocity * particle.dt

        # The ocean surface acts as a lid off, and if a particle goes above the ocean surface it is placed back at the
        # ocean surface (so at fieldset.SURF_Z)
        particle.potential = particle.depth + gradient + R + rise


def internal_tide_mixing(particle, fieldset, time):
    """
    Here we add vertical stochastic transport according to the internal tide mixing
    """
    if particle.beach == 0:
        # The grid of the tidal mixing isn't an exact match with the CMEMS data, so in regions where Kz_tidal is either 0
        # or very very small, we take the Waterhouse et al. (2014) estimate of the diapycnal diffusion below the MLD
        Kz_tidal = fieldset.TIDAL_Kz[time, particle.depth, particle.lat, particle.lon]
        dKz_tidal = fieldset.TIDAL_dKz[time, particle.depth, particle.lat, particle.lon]
        if Kz_tidal < fieldset.K_Z_BULK:
            Kz_tidal = fieldset.K_Z_BULK
            dKz_tidal = 0.0

        # Implementing Markov-0 style random walk
        tidal_gradient = dKz_tidal * particle.dt
        tidal_random = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3) * math.sqrt(2 * Kz_tidal)

        particle.potential = particle.depth + tidal_gradient + tidal_random


def KPP_TIDAL_mixing(particle, fieldset, time):
    """
    Combination of KPP and internal tide mixing. By running only one kernel I hope to make everything a bit more
    efficient computationally
    :param particle:
    :param fieldset:
    :param time:
    :return:
    """
    if particle.beach == 0:
        # Loading the mixed layer depth
        mld = fieldset.MLD[time, particle.depth, particle.lat, particle.lon]

        if particle.depth > mld:
            Kz = 0
            dKz = 0
        # Within the MLD we compute the vertical diffusion according to Boufadel et al. (2020)
        else:
            # Wind speed
            if particle.lon < 0:
                temp_lon = particle.lon + 360.
            else:
                temp_lon = particle.lon
            w_10 = math.sqrt(fieldset.u10[time, fieldset.SURF_Z, particle.lat, temp_lon] ** 2 + \
                             fieldset.v10[time, fieldset.SURF_Z, particle.lat, temp_lon] ** 2)
            # Drag coefficient
            C_D = min(max(1.2E-3, 1.0E-3 * (0.49 + 0.065 * w_10)), 2.12E-3)
            # wind stress
            tau = C_D * fieldset.RHO_A * w_10 ** 2
            # Frictional velocity of water at the ocean surface
            U_W = math.sqrt(tau / particle.surface_density)
            # Surface roughness z0 following Zhao & Li (2019)
            z0 = 3.5153e-5 * fieldset.BETA ** (-0.42) * w_10 ** 2 / fieldset.G
            # The corrected particle depth, since the depth is not always zero for the surface circulation data
            z_correct = particle.depth - fieldset.SURF_Z
            # The diffusion gradient at particle.depth
            alpha = (fieldset.VK * U_W) / (fieldset.PHI * mld ** 2)
            dKz = alpha * (mld - z_correct) * (mld - 3 * z_correct - 2 * z0)
            # The KPP profile vertical diffusion, at a depth corrected for the vertical gradient in Kz, and including
            # the bulk diffusivity
            alpha = (fieldset.VK * U_W) / fieldset.PHI
            Kz = alpha * (z_correct + z0) * math.pow(1 - z_correct / mld, 2)

        # Adding in the tidal mixing
        Kz += fieldset.TIDAL_Kz[time, particle.depth, particle.lat, particle.lon]
        dKz += fieldset.TIDAL_dKz[time, particle.depth, particle.lat, particle.lon]

        # The Markov-0 vertical transport from Grawe et al. (2012)
        gradient = dKz * particle.dt
        R = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3) * math.sqrt(2 * Kz)
        rise = particle.rise_velocity * particle.dt

        # The ocean surface acts as a lid off, and if a particle goes above the ocean surface it is placed back at the
        # ocean surface (so at fieldset.SURF_Z)
        particle.potential = particle.depth + gradient + R + rise


def vertical_reflecting_boundary(particle, fieldset, time):
    """
    Somewhat in contrast with the name, the reflecting boundary is only applied at the sea bottom. At the surface, we
    have that a particle is set at fieldset.SURF_Z if it crosses the surface boundary, as physically you would expect
    that a particle just stays at the surface if it is rising (instead of shooting farther up or being reflected down).

    This is meant to be used in conjuction with a vertical transport kernel (e.g. KPP_TIDAL_mixing), where the kernel
    calculates the potential depth in the next time step before the boundary condition is applied.
    """
    local_bathymetry = fieldset.bathymetry[time, fieldset.SURF_Z, particle.lat, particle.lon]
    if particle.potential < fieldset.SURF_Z:
        particle.depth = fieldset.SURF_Z
    elif particle.potential > local_bathymetry:
        overshoot = math.fabs(local_bathymetry - particle.potential)
        particle.depth = local_bathymetry - overshoot
    else:
        particle.depth = particle.potential
    if particle.depth < fieldset.SURF_Z:
        print('This is where it goes wrong')
        a=fieldset.SURF_Z
        b=particle.potential
        c=particle.depth
        print(a)
        print(b)
        print(c)
        print(local_bathymetry)


def sticky_seafloor_boundary(particle, fieldset, time):
    """
    If the particle crosses the bottom boundary (e.g. it is deeper than the bathymetry), then the particle is stuck
    and we set particle.beach = 3.
    """
    local_bathymetry = fieldset.bathymetry[time, fieldset.SURF_Z, particle.lat, particle.lon]
    if particle.potential < fieldset.SURF_Z:
        particle.depth = fieldset.SURF_Z
    elif particle.potential > local_bathymetry:
        particle.beach = 3
    else:
        particle.depth = particle.potential


def initial_estimate_particle_rise_velocity(L=settings.INIT_SIZE, print_rise=False):
    """
    Here we make an initial estimate of what the particle rise velocity is based on Enders et al. (2015)
    https://doi.org/10.1016/j.marpolbul.2015.09.027
    :param L: the particle size in meters (diameter, not radius)
    :return:
    """
    if type(L) in [float, int, np.float64, np.float32, np.int64, np.int32]:
        def to_optimize(w_rise):
            rho_p = settings.INIT_DENSITY  # Density of plastic particle
            rho_w = 1027  # rho sea water (kg/m^3)
            nu = 1.1e-6  # kinematic viscosity of sea water (Enders et al., 2015)
            left = (1. - rho_p / rho_w) * 8. / 3. * L * settings.G
            Re = 2. * L * np.abs(w_rise) / nu
            right = np.square(w_rise) * (24. / Re + 5. / np.sqrt(Re) + 2. / 5.)
            return np.abs(left - right)

        w_rise = scipy.optimize.minimize_scalar(to_optimize, bounds=[-100, 0], method='bounded').x
        print_statement("The rise velocity is for a particle with size {} is {}".format(L, w_rise), to_print=print_rise)
        return w_rise
    elif type(L) in [np.ndarray, list]:
        w_rise = deepcopy(L)
        for index_L, L_size in enumerate(L):
            def to_optimize(w_rise):
                rho_p = settings.INIT_DENSITY  # Density of plastic particle
                rho_w = 1027  # rho sea water (kg/m^3)
                nu = 1.1e-6  # kinematic viscosity of sea water (Enders et al., 2015)
                left = (1. - rho_p / rho_w) * 8. / 3. * L_size * settings.G
                Re = 2. * L_size * np.abs(w_rise) / nu
                right = np.square(w_rise) * (24. / Re + 5. / np.sqrt(Re) + 2. / 5.)
                return np.abs(left - right)
            w_rise[index_L] = scipy.optimize.minimize_scalar(to_optimize, bounds=[-100, 0], method='bounded').x
        return w_rise
    else:
        print_statement("the input type is {}, and the value is {}".format(type(L), L), to_print=True)


def initial_reynolds_number(L=settings.INIT_SIZE, rise_velocity=None):
    if rise_velocity is None:
        rise_velocity = initial_estimate_particle_rise_velocity(L=L)
    nu = 1.1e-6  # kinematic viscosity of sea water (Enders et al., 2015)
    return L * np.abs(rise_velocity) / nu


def get_resuspension_timescale(L=settings.INIT_SIZE, w_rise=None, print_size=False):
    """
    This follows equation 9 from Hinata et al. (2017)
    https://doi.org/10.1016/j.marpolbul.2017.05.012
    :param L:
    :param print_size:
    :return:
    """
    if w_rise is None:
        w_rise = initial_estimate_particle_rise_velocity(L=L)
    lambda_R = 2.6e2 * np.abs(w_rise) + 7.1
    print_statement("The resuspension timescale for a particle of size {} is {:.6f} days".format(settings.INIT_SIZE, np.mean(lambda_R)), to_print=print_size)
    return lambda_R


def resuspension_probability(w_rise):
    """
    This computes the resuspension probability for a particle with size L
    :param L:
    :return:
    """
    lambda_resus = get_resuspension_timescale(L=0, w_rise=w_rise)
    prob_resus = np.exp(-settings.TIME_STEP.total_seconds() / (np.array(lambda_resus) * 86400.))
    return prob_resus


def get_reynolds_number(particle, fieldset, time):
        """
        Calculating the reynolds number (https://en.wikipedia.org/wiki/Reynolds_number)
        """
        w_b = math.fabs(particle.rise_velocity)
        particle.reynolds = particle.size * w_b / particle.kinematic_viscosity


def get_rising_velocity(particle, fieldset, time):
    """
    Calculating the rise velocity of a particle following the Enders et al. (2015) paper
    https://doi.org/10.1016/j.marpolbul.2015.09.027
    :return:
    """
    rho_sw = particle.density  # sea water rho (kg m^-3)
    rho_p = particle.rho_plastic  # plastic particle rho (kg m^-3)
    left = (1. - rho_p / rho_sw) * 8. / 3. * particle.size * fieldset.G
    right = 24. / particle.reynolds + 5. / math.sqrt(particle.reynolds) + 2. / 5.
    if left > 0:
        particle.rise_velocity = - 1 * math.sqrt(left / right)
    else:
        particle.rise_velocity = math.sqrt(math.fabs(left) / right)


def TotalDistance(particle, fieldset, time):
    """
    Calculating the cumulative distance travelled by the particle in vertical and horizontal directions
    """
    # Calculate the distance in latitudinal direction (using 1.11e2 kilometer per degree latitude)
    lat_dist = (particle.lat - particle.prev_lat) * 1.11e2
    # Calculate the distance in longitudinal direction, using cosine(latitude) - spherical earth
    lon_dist = (particle.lon - particle.prev_lon) * 1.11e2 * math.cos(particle.lat * math.pi / 180)
    # Calculate the total Euclidean distance travelled by the particle
    particle.distance_horizontal += math.sqrt(math.pow(lon_dist, 2) + math.pow(lat_dist, 2))
    particle.distance_vertical += math.fabs(particle.depth - particle.prev_depth)

    particle.prev_lon = particle.lon  # Set the stored values for next iteration.
    particle.prev_lat = particle.lat
    particle.prev_depth = particle.depth


def mass_per_size_class(k, f, p=settings.P_FRAG):
    gamma_ratio = math.gamma(k + f) / (math.gamma(k + 1) * math.gamma(f))
    return gamma_ratio * p ** k * (1 - p) ** f


def particle_number_per_size_class(k, f=1.0, p=settings.P_FRAG, Dn=settings.DN):
    return mass_per_size_class(k, f, p) * 2 ** (Dn * k)


def seabed_resuspension_beaching_kernel(particle, fieldset, time):
    """
    The beaching and resuspension kernels for beaching on the coastline follows the procedure outlined in Onink et
    al. (2021).
    Onink et al. (2021) = https://doi.org/10.1088/1748-9326/abecbd

    For particles on the seabed, we follow the resuspension procedure outlined in Carvajalino-Fernandez et al.
    (2020), where a particle at the sea bed gets resuspended if the estimated sea floor sea stress is greater than a
    critical threshold. Particles can get stuck on the seabed if the potential depth due to KPP or internal tide
    mixing is below the bathymetry depth

    The bottom sea stress is calculated using a quadratic drag extrapolation according to Warner et al. (2008)
    Carvajalino-Fernandez et al. (2020) = https://doi.org/10.1016/j.marpolbul.2020.111685
    Warner et al. (2008) = https://doi.org/10.1016/j.cageo.2008.02.012

    Note: This beaching kernel was developed for 3D size dependent transport runs, but given the high uncertainty in how
    to represent the seabed resuspension it was ultimately not used.
    """
    # First, the beaching of particles on the coastline
    if particle.beach == 0:
        dist = fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
        if dist < fieldset.Coastal_Boundary:
            if ParcelsRandom.uniform(0, 1) > fieldset.p_beach:
                particle.beach = 1
    # Next the resuspension of particles on the coastline
    elif particle.beach == 1:
        if ParcelsRandom.uniform(0, 1) > fieldset.p_resus:
            particle.beach = 0
    # Finally, the resuspension of particles on the seabed
    elif particle.beach == 3:
        dWx = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
        dWy = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)

        bx = math.sqrt(2 * fieldset.SEABED_KH)

        # Getting the current strength at the particle position at the sea bed, and converting it to m/s
        U_bed, V_bed = fieldset.U[time, particle.depth, particle.lat, particle.lon], fieldset.V[time, particle.depth, particle.lat, particle.lon]
        U_bed, V_bed = U_bed * 1852. * 60. * math.cos(40. * math.pi / 180.), V_bed * 1852. * 60.
        U_bed, V_bed = U_bed + bx * dWx, V_bed + bx * dWy
        # Getting the bottom shear stress
        tau_bss = 0.003 * (math.pow(U_bed, 2) + math.pow(V_bed, 2))
        # if tau_bss is greater than fieldset.SEABED_CRIT, then the particle gets resuspended
        if tau_bss > fieldset.SEABED_CRIT:
            particle.beach = 0
    # Update the age of the particle
    particle.age += particle.dt