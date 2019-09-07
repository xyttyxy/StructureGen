from catkit.gen.surface import SlabGenerator
from ase.build import bulk
from ase.visualize import view
from ase import Atom, Atoms
from ase.calculators.vasp import Vasp, Vasp2
import numpy as np
import copy

# rules: 
# add first, if distance less than eps, reject the structure
def next_ads(picked_struct, size, cell, sites):
    returned = []
    for i, c in enumerate(sites):
        for a in range(size[0]):
            for b in range(size[1]):
                candidate = copy.deepcopy(picked_struct)
                candidate += Atom('N', c+[0,0,2]+a*cell[0]+b*cell[1])
                distances = candidate.get_all_distances(mic=True)[-1]
                if min(distances[0:-1]) < 2:
                    pass
                else: 
                 returned.append(candidate)
    return returned

atoms_bulk = bulk(name='Cu', crystalstructure='fcc', a=3.6302862146117354, cubic=True)

gen = SlabGenerator(atoms_bulk,
                    miller_index=(1,1,1),
                    layers=4,
                    fixed=2,
                    vacuum=4,
                    standardize_bulk=True,
                    layer_type="trim")
size = [3,4]
primitive = gen.get_slab()
bare = gen.get_slab(size=size)

coordinates, connectivity= gen.adsorption_sites(primitive)
cell_a = primitive.get_cell()[0]
cell_b = primitive.get_cell()[1]

# adsorption of 1 nitrogen, different because only 1st primitive cell need to be considered
Cu111_1N = next_ads(bare, [1,1], primitive.get_cell(), coordinates)
print("1N Adsorption: %d structures" % len(Cu111_1N))

for idx, atoms in enumerate(Cu111_1N):
    label = 'Cu_N_ads/111/1N/s'+str(idx+1)+'/'
    calc = Vasp2(label=label,
                 encut=400,
                 ismear=2,
                 sigma=0.2,
                 gga='PE',
                 xc='PBE',
                 algo='FAST',
                 ispin=2,
                 ediff=1E-6,
                 ediffg=-2E-2,
                 ibrion=1,
                 nsw=200,
                 lreal='AUTO',
                 kpts=[5,5,1])
    calc.write_input(atoms)
    
# for each of these atoms object generate a folder

## eliminating top and bridge sites
#coordinates = coordinates[2:3]
#connectivity = connectivity[2:4]

#######################################################################
# pretend you calculated and found lowest energy
# now add the second nitrogen
#Cu111_1N_pick = Cu111_1N[3]
#Cu111_2N = next_ads(Cu_N_1_pick, size, primitive.get_cell(), coordinates)
## view(Cu_N_2)
#print("2N Adsorption: %d structures" % len(Cu111_2N))

# TODO: PRESERVE THE MS PROJECT FOLDER STRUCTURE
# NAMELY:
# Cu_N_ads/<surface termination>/<number of nitrogen adsorbed>/<configuration>
# Use the FileIOCalculator. Learn based on Gaussian

########################################################################
## pretend you calculated and found lowest energy
## now add the third nitrogen
#Cu_N_2_pick = Cu_N_2[0]
#Cu_N_3 = next_ads(Cu_N_2_pick, size, primitive.get_cell(), coordinates)
## view(Cu_N_3)
#print("3N Adsorption: %d structures" % len(Cu_N_3))

########################################################################
## pretend you calculated and found lowest energy
## now add the third nitrogen
#Cu_N_3_pick = Cu_N_3[0]
#Cu_N_4 = next_ads(Cu_N_3_pick, size, primitive.get_cell(), coordinates)
#view(Cu_N_4)
#print("4N Adsorption: %d structures" % len(Cu_N_4))

#print("Testing update 2")