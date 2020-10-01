from parcels import FieldSet, JITParticle

from src.scenarios.base_scenario import BaseScenario
from src.factories.fieldset_factory import FieldSetFactory
from src.factories.particle_factory import ParticleFactory


class CoastalProximity(BaseScenario):
    """Coastal proximity scenario"""

    def __init__(self, server, stokes):
        """Constructor for coastal_proximity"""
        super().__init__(server, stokes)

    def create_fieldset(self) -> FieldSet:
        fieldset = FieldSetFactory().create_fieldset(server=self.server, stokes=self.stokes,
                                                     stokes_drift=False, border_current=True, diffusion=True)
        return fieldset

    def create_particle(self) -> JITParticle:
        particle = ParticleFactory().create_particle()
        return particle

    def run(self) -> None:
        pass







