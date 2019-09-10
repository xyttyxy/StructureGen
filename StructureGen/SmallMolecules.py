from ase.build import molecule
from ase import Atoms
from ase.visualize import view
from ase.calculators.vasp import Vasp2
label = ['CH4', 'NH3', 'H2O', 'HF', 'SiH4', 'PH3', 'SH2', 'HCl']
molecules = []

for i in label:
	temp = molecule(name=i,cell=[10,11,12],pbc=True)
	temp.center()
	calc = Vasp2(label=i,
			  sigma=0,
			  ismear=0.01,
			  ispin=2,
			  encut=400,
			  ibrion=1,
			  nsw=100,
			  ediff=1E-4,
			  ediffg=-2e-2)

	temp.set_calc(calc)
	molecules.append(temp)


for i in molecules:
	energy = i.get_potential_energy()
	print(energy)

