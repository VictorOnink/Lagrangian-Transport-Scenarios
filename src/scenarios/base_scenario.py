from abc import ABC, abstractmethod

from parcels import FieldSet, JITParticle


class BaseScenario(ABC):
    server: int
    stokes: int
    field_set: FieldSet
    particle: JITParticle

    """A base class for the different scenarios"""
    def __init__(self, server, stokes):
        self.server = server
        self.stokes = stokes
        self.field_set = self.create_fieldset()
        self.particle = self.create_particle()

    @abstractmethod
    def create_fieldset(self) -> FieldSet:
        pass

    @abstractmethod
    def create_particle(self) -> JITParticle:
        pass

    @abstractmethod
    def run(self):
        pass







