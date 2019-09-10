from catkit.gen.surface import SlabGenerator
from ase.build import bulk, make_supercell
from ase.visualize import view
from ase import Atom, Atoms
from ase.calculators.vasp import Vasp, Vasp2
import numpy as np
import copy, os, subprocess

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

terminations   = ['100', '110', '111', '210', '211', '221', '311']
miller_indices = [(1,0,0),(1,1,0),(1,1,1),(2,1,0),(2,1,1), (2,2,1), (3,1,1)]
layers = [4, 5, 4, 8, 8, 10, 6]
fixed  = [2, 3, 2, 4, 4, 5, 3]
kpoints = [[5,3,1],[5,3,1],[5,5,1],[5,3,1],[7,3,1],[7,1,1],[5,4,1]]
sizes = [(3,4,1), (3,3,1), (3,3,1), (2,2,1),(2,2,1), (2,2,1), (3,2,1)]
vacuum = 14

primitive = []
bare = []
coordinates = []
connectivity = []


Cu_N_ads = {}

for idx, term in enumerate(terminations):
    gen = SlabGenerator(atoms_bulk, 
                        miller_index = miller_indices[idx],
                        layers = layers[idx],
                        fixed = fixed[idx],
                        vacuum = 10,
                        standardize_bulk = True,
                        layer_type = "trim")

    prim = gen.get_slab()
    supercell = prim.repeat( sizes[idx] )
    #print(supercell.get_cell())
    primitive.append( prim )
    bare.append( supercell )
    coords_tmp, connect_tmp = gen.adsorption_sites( prim )
    coordinates.append( coords_tmp )
    connectivity.append( connect_tmp )

    list_ads = next_ads(bare[-1], [1,1], primitive[-1].get_cell(), coordinates[-1])
    Cu_N_ads.update( {term : {'1N': list_ads}} )
    print("Termination (" + term + ") 1N Adsorption: %d structures" % len(Cu_N_ads[term]['1N']))

# Data structure
# Cu_N_ads = {'111' : {'1N' : [s1, s2, s3, ...]}}

calc_template = Vasp2(encut=400,
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

for termination, value_1 in Cu_N_ads.items():
    for num_ads, value_2 in value_1.items():
        t = 0
        for i, slab in enumerate(value_2):
            label = 'Cu_N_ads/'+termination+'/'+num_ads+'/t'+str(t)+'/'
            calc = copy.deepcopy(calc_template)
            calc.set_label(label)
            calc.set(kpts= kpoints[terminations.index(termination)])
            slab.set_calculator(calc)
            calc.write_input(slab)
            # do stuff here
            t += 1

# submit
cwd = os.getcwd()
for termination, value_1 in Cu_N_ads.items():
    for num_ads, value_2 in value_1.items():
        t = 0
        for i, slab in enumerate(value_2):
            workdir = 'Cu_N_ads/'+termination+'/'+num_ads+'/t'+str(t)+'/'
            os.chdir(workdir)
            #subprocess.run("qsub ~/scripts/job_scripts/direct/job_dc_24_ex.sh")
            #print(os.getcwd())
            os.chdir(cwd)
            # do stuff here
            t += 1

    
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