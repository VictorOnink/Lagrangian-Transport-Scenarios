from parcels import FieldSet, ParticleSet
import numpy as np
import settings as settings
import scenarios.base_scenario as base_scenario
import factories.fieldset_factory as fieldset_factory
from advection_scenarios import advection_files
import utils as utils
from datetime import datetime, timedelta
import os
from parcels import ParcelsRandom
import math


class Turrell_Resuspension(base_scenario.BaseScenario):
    """Stochastic beaching and shore dependent resuspension """

    def __init__(self, server, stokes):
        """Constructor for coastal_proximity"""
        super().__init__(server, stokes)
        self.prefix = "Turrell"
        self.input_dir = utils._get_input_directory(server=self.server)
        self.output_dir = utils._get_output_directory(server=self.server)
        if settings.RESTART == 0:
            self.repeat_dt = timedelta(days=31)
        else:
            self.repeat_dt = None
        if settings.SUBMISSION == 'simulation':
            advection_scenario = advection_files.AdvectionFiles(server=self.server, stokes=self.stokes,
                                                                advection_scenario=settings.ADVECTION_DATA,
                                                                repeat_dt=self.repeat_dt)
            self.file_dict = advection_scenario.file_names
            self.field_set = self.create_fieldset()

    var_list = ['lon', 'lat', 'beach', 'age', 'weight']

    def create_fieldset(self) -> FieldSet:
        os.system('echo "Creating the fieldset"')
        fieldset = fieldset_factory.FieldSetFactory().create_fieldset(file_dict=self.file_dict, stokes=self.stokes,
                                                                      border_current=True, diffusion=True, landID=True,
                                                                      distance=True,beach_timescale=True,wind=True,
                                                                      sea_elev=True, wind_min=True
                                                                      )
        return fieldset

    def _get_pset(self, fieldset: FieldSet, particle_type: utils.BaseParticle, var_dict: dict,
                  start_time: datetime, repeat_dt: timedelta):
        """
        :return:
        """
        os.system('echo "Creating the particle set"')
        pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                           lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
                           age=var_dict['age'], weights=var_dict['weight'],
                           time=start_time, repeatdt=repeat_dt)
        return pset

    def _get_pclass(self):
        os.system('echo "Creating the particle class"')
        particle_type = utils.BaseParticle
        utils._add_var_particle(particle_type, 'distance', dtype=np.float32, set_initial=False)
        return particle_type

    def _file_names(self, new: bool = False, run: int = settings.RUN, restart: int = settings.RESTART):
        odirec = self.output_dir+"Turrell/st_"+str(settings.SHORE_TIME)+"_W_"+str(settings.WMIN)+\
                 "_e_"+str(settings.ENSEMBLE)+"/"
        if new == True:
            os.system('echo "Set the output file name"')
            return odirec + self.prefix +'_{}'.format(settings.ADVECTION_DATA) + "_Wmin="+str(settings.WMIN)+\
                   "_st="+str(settings.SHORE_TIME)+"_y=" +str(settings.START_YEAR)+"_I="+str(settings.INPUT)+\
                   "_r="+str(restart)+"_run="+str(run)+".nc"
        else:
            os.system('echo "Set the restart file name"')
            return odirec + self.prefix +'_{}'.format(settings.ADVECTION_DATA) + "_Wmin="+str(settings.WMIN)+\
                   "_st="+str(settings.SHORE_TIME)+"_y=" +str(settings.START_YEAR)+"_I="+str(settings.INPUT)+\
                   "_r="+str(restart)+"_run="+str(run)+".nc"

    def _beaching_kernel(particle, fieldset, time):
        """
            Beaching is implemented the same way as in the stochastic and shore dependent resuspension
            scenarios.

            Resuspension is based on Turrell 2018 & 2020. Resuspension is possible when
            water levels are at the same level as that of the beached particle. Then,
            only when the offshore wind component is greater than the threshold Wmin
            will the particle actually be resuspended
            """
        t, d, la, lo = time, particle.depth, particle.lat, particle.lon
        # Beaching
        if particle.beach == 0:
            dist = fieldset.distance2shore[t, d, la, lo]
            if dist < fieldset.Coastal_Boundary:
                if ParcelsRandom.uniform(0, 1) > fieldset.p_beach:
                    particle.beach = 1
                    particle.depth = fieldset.eta[t, d, la, lo]
        # Resuspension
        elif particle.beach == 1:
            sea_elev = fieldset.eta[t, d, la, lo]
            # If particles are beached above sea level, then they will remain there
            if particle.depth < sea_elev:
                # particles will get pushed up to the present water level
                particle.depth = sea_elev
                # Now, we need to get the offshore wind component, for which we first get the direction
                # of the border current. We convert the border current back to m/s and multiply by -1 to reverse
                #the sign such that the current is directed offshore
                bU, bV = fieldset.borU[t, d, la, lo] * -1 * 1852 * 60 * math.cos(la * math.pi / 180), fieldset.borV[
                    t, d, la, lo] * -1 * 1852 * 60
                wU, wV = fieldset.U[t, d, la, lo] * 1852 * 60 * math.cos(la * math.pi / 180), fieldset.V[
                    t, d, la, lo] * 1852 * 60
                # magnitude of the b and w vectors, and then dot product between them
                mB, mW = math.sqrt(bU ** 2 + bV ** 2), math.sqrt(wU ** 2 + wV ** 2)
                dot = bU * wU + bV * wV
                # Angle between b and w
                alpha = math.acos(dot / (mB * mW))
                # If the angle between thexs wind and border current is <90 degrees,
                # then the wind is directed offshore
                if alpha < math.pi / 2:
                    # If the offshore component of wind is greater than Wmin, then
                    # the particle gets resuspended
                    if mW * math.cos(alpha) < fieldset.Wmin:
                        particle.beach = 0
        # Update the age of the particle
        particle.age += particle.dt

    def _get_particle_behavior(self, pset: ParticleSet):
        os.system('echo "Setting the particle behavior"')
        base_behavior = pset.Kernel(utils._initial_input) + pset.Kernel(utils._floating_advection_rk4) + \
                        pset.Kernel(utils._floating_2d_brownian_motion)
        total_behavior = base_behavior + pset.Kernel(utils._anti_beach_nudging) + pset.Kernel(self._beaching_kernel)
        return total_behavior
