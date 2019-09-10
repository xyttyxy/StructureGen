#from catkit.gen.surface import SlabGenerator
from ase.build import fcc100, fcc111_root, add_adsorbate, make_supercell
from ase.visualize import view
from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp2
import copy, os, subprocess

ads_height = 2
def slabgen(termination, size, adsorbate, position):
    if termination == '100':
        prim = fcc100('Cu', a=3.6302862146117354, size=(1,1,5),vacuum=15, orthogonal=True, periodic=True)
    elif termination == '111':
        prim = fcc111_root('Cu', root=3, a=3.6302862146117354, size=(1,1,5),vacuum=15)

    super = make_supercell(prim, size)
    add_adsorbate(slab=super, adsorbate=adsorbate, height=ads_height, position=position)
    constr = FixAtoms(indices=[atom.index for atom in super if atom.position[2] < 19])
    super.set_constraint(constr)
    return super

modifiers = ['C', 'Si', 'O', 'S', 
             'F', 'Cl', 'N', 'P']

terminations = ['100', '100', '100', '100',
               '100', '100', '111', '111']

positions = ['hollow', 'hollow', 'hollow', 'hollow', 
             'hollow', 'hollow', (1.853,1.925), (1.853,1.925)]
sizes = [[[2,0,0],[0,2,0],[0,0,1]], [[2,0,0],[0,2,0],[0,0,1]],
         [[2,0,0],[0,1,0],[0,0,1]], [[2,0,0],[0,1,0],[0,0,1]],
         [[1,0,0],[0,1,0],[0,0,1]], [[1,0,0],[0,1,0],[0,0,1]],
         [[1,0,0],[0,1,0],[0,0,1]], [[1,0,0],[0,1,0],[0,0,1]]]
kpts = [[7,7,1], [7,7,1], [7,13,1], [7,13,1], 
        [13,13,1], [13,13,1], [9,9,1], [9,9,1]]
slabs = []
calc_template = Vasp2(encut=400,
             ismear=2,
             sigma=0.2,
             gga='PE',
             xc='PBE',
             ibrion=1,
             algo='FAST',
             potim=0.5,
             nsw=200,
             ediff=1E-6,
             ediffg=-2E-2,
             lreal='AUTO',
             kpts=[7,7,1])

# write input files
for idx, mod in enumerate(modifiers):
    calc = copy.deepcopy(calc_template)
    workdir = mod+'/'
    calc.set_label(workdir)
    calc.set(kpts=kpts[idx])
    slab = slabgen(terminations[idx],sizes[idx], mod, positions[idx])
    slab.set_calculator(calc)
    calc.write_input(slab)
    slabs.append(slab)
    #print(slabs[-1].get_calculator().label)

# execute
cwd = os.getcwd()
for slab in slabs:
    workdir = slab.get_calculator().label
    os.chdir(workdir)
#    subprocess.run("qsub ~/scripts/job_scripts/direct/job_dc_24_ex.sh")
    print(os.getcwd())
    os.chdir(cwd)