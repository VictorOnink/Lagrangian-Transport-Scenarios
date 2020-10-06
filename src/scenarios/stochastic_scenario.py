from parcels import ParticleSet, JITParticle, FieldSet

from scenarios.base_scenario import BaseScenario


class StochasticSeneratio(BaseScenario):
    """A stocastic scenario"""

    def __init__(self):
        """Constructor for StochasticSeneratio"""

    def create_fieldset(self) -> FieldSet:
        pass

    def create_particle_set(self) -> JITParticle:
        pass

    def _get_pset(self) -> ParticleSet:
        pass

    def _file_names(self, input_dir: str, new: bool) -> str:
        pass

    def _get_start_var(self) -> dict:
        pass

    def _beaching_kernel(self) -> ParticleSet:
        pass

    def _get_particle_behavior(self):
        pass


