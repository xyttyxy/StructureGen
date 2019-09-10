from catkit.gen.surface import SlabGenerator
from ase.build import bulk
from ase.visualize import view
from ase import Atom, Atoms
import numpy as np

atoms_bulk = bulk(name='Cu', crystalstructure='fcc', a=3.6302862146117354, cubic=True)
gen = SlabGenerator(atoms_bulk,
                    miller_index = (2,1,0),
                    layers = 8,
                    fixed = 4,
                    vacuum = 10,
                    standardize_bulk = True,
                    layer_type = "trim")

primitive = gen.get_slab()
#supercell = primitive.repeat(2,3)
supercell_1 = gen.get_slab( size=2 )
supercell_2 = gen.get_slab( size=(2,2) )
supercell_3 = gen.get_slab( size=np.array([[2,0],[0,2]]))
view(primitive)
view(supercell_3)