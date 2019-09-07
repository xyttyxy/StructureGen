from ase.calculators.vasp import Vasp, Vasp2
from ase.build import bulk
from ase.visualize import view
from ase import Atom, Atoms
from ase.calculators.vasp import Vasp
import copy

d = 1.1
co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)], cell=[10,10,100],pbc=True)
calc = Vasp(encut=400,
            ismear=2,
            sigma=0.2,
            gga='PE',
            xc='PBE')
label = "CarbonMonoxide"

co2 = copy.deepcopy(co)
calc2 = Vasp2(label=label,
              encut=400,
              ismear=2,
              sigma=0.2,
              gga='PE',
              xc='PBE')

co.set_calculator(calc)
co.get_potential_energy()
