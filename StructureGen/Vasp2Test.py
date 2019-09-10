from ase.calculators.vasp import Vasp, Vasp2
from catkit.gen.surface import SlabGenerator
from ase.build import bulk
from ase.visualize import view
from ase import Atom, Atoms
import copy

#d = 1.1
#co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)], cell=[10,10,100],pbc=True)
#calc = Vasp(encut=400,
#            ismear=2,
#            sigma=0.2,
#            gga='PE',
#            xc='PBE')
#label = "CarbonMonoxide"

#co2 = copy.deepcopy(co)
#calc2 = Vasp2(label=label,
#              encut=400,
#              ismear=2,
#              sigma=0.2,
#              gga='PE',
#              xc='PBE')

#co.set_calculator(calc)
#co.get_potential_energy()

atoms_bulk = bulk(name='Cu', crystalstructure='fcc', a=3.6302862146117354, cubic=True)
gen = SlabGenerator(atoms_bulk,
                    miller_index = [2,1,0],
                    layers = 8,
                    fixed = 4,
                    vacuum = 10,
                    standardize_bulk = True,
                    layer_type = "trim")
bare = gen.get_slab(size=[2,2])
view(bare)