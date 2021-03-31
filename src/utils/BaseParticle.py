from operator import attrgetter
import numpy as np
from parcels import JITParticle, Variable


class BaseParticle(JITParticle):
    # 0=open ocean, 1=beached
    beach = Variable('beach', dtype=np.int32, initial=attrgetter('beach'))
    # Finally, I want to keep track of the age of the particle
    age = Variable('age', dtype=np.int32, initial=attrgetter('age'))
    # # Weight of the particle in tons
    # weights = Variable('weights', dtype=np.float32, initial=attrgetter('weights'))
