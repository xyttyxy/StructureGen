from catkit.gen.surface import SlabGenerator
from ase.build import bulk
from ase.visualize import view

atoms_bulk = bulk('Cu', 'fcc', cubic=True)
gen = SlabGenerator(atoms_bulk,
                    miller_index=(1,1,1),
                    layers=4,
                    fixed=2,
                    vacuum=4,
                    standardize_bulk=True,
                    layer_type="trim")

primitive = gen.get_slab()
supercell = gen.get_slab(size=(3,4))
view(supercell,viewer="RASMOL")