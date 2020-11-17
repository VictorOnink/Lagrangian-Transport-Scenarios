import math
from parcels import rng as random


def _anti_beach_nudging(particle, fieldset, time):
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


def _floating_advection_rk4(particle, fieldset, time):
    """Advection of particles using fourth-order Runge-Kutta integration.

    Function needs to be converted to Kernel object before execution

    A particle only moves if it has not beached (rather obviously)
    """
    if particle.beach == 0:
        particle.distance = fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
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


def _floating_2d_brownian_motion(particle, fieldset, time):
    """Kernel for simple Brownian particle diffusion in zonal and meridional direction.
    Assumes that fieldset has fields Kh_zonal and Kh_meridional
    we don't want particles to jump on land and thereby beach"""
    if particle.beach == 0:
        # Wiener increment with zero mean and std of sqrt(dt)
        dWx = random.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
        dWy = random.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)

        bx = math.sqrt(2 * fieldset.Kh_zonal[time, particle.depth, particle.lat, particle.lon])
        by = math.sqrt(2 * fieldset.Kh_meridional[time, particle.depth, particle.lat, particle.lon])

        particle.lon += bx * dWx
        particle.lat += by * dWy


def _floating_AdvectionRK4DiffusionEM_stokes_depth(particle, fieldset, time):
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
        w_p = 2 * math.pi / fieldset.WP[t, d, la, lo]  # Peak period
        k_p = w_p ** 2 / 9.81  # peak wave number
        z_correc = max(d - 1.472102, 0)  # correction since the surface in the CMEMS Mediterranean data is not at 0
        st_z = math.exp(-2 * k_p * z_correc) - math.sqrt(2 * math.pi * k_p * z_correc) * math.erfc(2 * k_p * z_correc)

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
        dWx = random.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
        dWy = random.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)

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


def _initial_input(particle, fieldset, time):
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
            potential_lat = particle.lat + random.uniform(-d_lat, d_lat)
            potential_lon = particle.lon + random.uniform(-d_lon, d_lon)
            potential_land = math.floor(fieldset.landID[time, particle.depth, particle.lat, particle.lon])
            potential_dist = fieldset.distance2shore[time, particle.depth, potential_lat, potential_lon]
            if potential_land == 0 and potential_dist <= distCur:
                check += 100001
                particle.lat = potential_lat
                particle.lon = potential_lon
            check += 1


def _delete_particle(particle, fieldset, time):
    # This delete particle format from Philippe Delandmeter
    # https://github.com/OceanParcels/Parcelsv2.0PaperNorthSeaScripts/blob/master/northsea_mp_kernels.py
    print("Particle [%d] lost !! (%g %g %g %g)" % (
    particle.id, particle.lon, particle.lat, particle.depth, particle.time))
    particle.delete()


def _get_kinematic_viscosity(particle, fieldset, time):
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
    rho_sw = particle.density  # seawater density[kg m-3]
    rise = particle.rise_velocity  # vertical velocity[m s-1]

    # ------ Constants -----
    g = 7.32e10 / (86400. ** 2.)  # gravitational acceleration (m d-2), now [s-2]

    # ------ Diffusivity -----
    r_tot = particle.size  # total radius [m]
    rho_tot = (particle.size ** 3. * particle.rho_plastic) / (particle.size) ** 3.  # total density [kg m-3]

    dn = 2. * r_tot  # equivalent spherical diameter [m]
    delta_rho = (
                            rho_tot - rho_sw) / rho_sw  # normalised difference in density between total plastic+bf and seawater[-]
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

    # z0 = z + vs * particle.dt
    # # # 1.472102
    # if z0 <= 1.472102 or z0 >= fieldset.bathymetry[time, particle.depth, particle.lat, particle.lon]:
    #     vs = 0
    # else:
    #     particle.depth = z0

    particle.rise_velocity = vs