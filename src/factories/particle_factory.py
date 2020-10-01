import os
import src.settings as settings
from netCDF4 import Dataset
import numpy as np
from parcels import GeographicPolar, Geographic, FieldSet, Field, JITParticle


class ParticleFactory:
    """A factory class for the particles"""
    @classmethod
    def create_particle(cls) -> JITParticle:
        """

        :rtype: object
        """
        particle = _create_base_particle()
        return particle


def _create_base_particle() -> JITParticle:
    pass




