from ase.io import read
import numpy as np
from ase.io.pov import write_pov # , get_bondpairs
from pathlib import Path

atoms = read('./POSCAR-1', format='vasp')
atoms.wrap()
# atoms = atoms.repeat([3, 3, 1])
pov_file_name = 'POSCAR.pov'    # suffix must be .pov
ini_file_name = Path(pov_file_name).with_suffix('.ini')
png_file_name = 'POSCAR.png'
self_define_colors = False
atoms_show_bonds_model = 3    # 1, all atoms show bonds; 2, no atom show bonds; 3, only selected atoms show bonds
bondlinewidth = 0.1

def _self_define_get_bondpairs(atoms, indices_of_interest:list, radius=1.2):
    from ase.data import covalent_radii
    from ase.neighborlist import NeighborList
    cutoffs = radius * covalent_radii[atoms.numbers]
    nl = NeighborList(cutoffs=cutoffs, self_interaction=False)
    nl.update(atoms)
    bondpairs = []
    for a in indices_of_interest:
        indices, offsets = nl.get_neighbors(a)
        bondpairs.extend([(a, a2, offset)                           # 除了a要在列表里, a2也要在列表里
                          for a2, offset in zip(indices, offsets) if a2 in indices_of_interest])
    return bondpairs

bond_atom_index = [i for i, atom in enumerate(atoms) if atom.symbol in ['H','O'] if i not in [305, 386,387]]
bondatoms_list =  _self_define_get_bondpairs(atoms, bond_atom_index, radius=0.50) # []

plotting_var_settings = {'maxwidth': 50000000000000, "show_unit_cell": 0, "rotation": '-0z, -0x, 0y'}
if atoms_show_bonds_model == 1:
    radii_scale = 0.57          # float, radii = covalent_radii*radii_scale
elif atoms_show_bonds_model == 2:
    radii_scale = None          # None, radii = covalent_radii
elif atoms_show_bonds_model == 3:
    from ase.data import covalent_radii     # a numpy array
    radii_scale = 0.57          # list, radii = self define list
    radii_list = [bondlinewidth # covalent_radii[atomic_number]*radii_scale or bondlinewidth
                  if index in bond_atom_index else covalent_radii[atomic_number] 
                  for index, atomic_number in enumerate(atoms.get_atomic_numbers())]
    radii_scale = radii_list
else:
    raise ValueError(f"atoms_show_bonds_model must be 1, 2 or 3, while you set {atoms_show_bonds_model}")
plotting_var_settings['radii'] = radii_scale


povray_settings = {'canvas_width':1000, "display":True, "transparent":True, "textures":['jmol'] * len(atoms), 
                   "bondlinewidth":bondlinewidth, "bondatoms":bondatoms_list}
""" canvas_width is image width, one canot set image height and width at the same time; 
Display the image while rendering, default is False;
Transparent background, default is True;
textures default is ase3;
bondatoms default is [], one can set to get_bondpairs(atoms) """


transmittances = np.zeros(len(atoms))  # 初始化所有原子的透明度为0
h_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'H']
o_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
# 设置除了索引为 307, 388, 389, 390 之外的 H 和 O 原子的透明度
for i in h_indices + o_indices:
    if i not in [305, 386,387]:
        transmittances[i] = 0.01
povray_settings["transmittances"] = transmittances  # 将 transmittances 赋值给 povray_settings


if self_define_colors:
    from ase.data.colors import jmol_colors
    ncolors = len(jmol_colors)
    colors = jmol_colors[atoms.get_atomic_numbers().clip(max=ncolors - 1)]
    # colors[1] = np.array([0.122, 0.941, 0.122])         # Change the color of the 2th atom to green
    povray_settings["colors"] = colors


#=======================Method 1: Use pre setting=======================
pov_ini_put = write_pov(pov_file_name, atoms, povray_settings=povray_settings, **plotting_var_settings)
# For function write_pov, param povray_settings is param for Class POVRAY.
'''povray_settings include:
cell, cell_vertices, positions, diameters, colors,
image_width, image_height, constraints=tuple(), isosurfaces=[],
display=False, pause=True, transparent=True, canvas_width=None,
canvas_height=None, camera_dist=50., image_plane=None,
camera_type='orthographic', point_lights=[],
area_light=[(2., 3., 40.), 'White', .7, .7, 3, 3],
background='White', textures=None, transmittances=None,
depth_cueing=False, cue_density=5e-3,
celllinewidth=0.05, bondlinewidth=0.10, bondatoms=[],
exportconstraints=False

But!!! Usually, one can use classmethod POVRAY.from_PlottingVariables to init.
If one want change textures from ase3 to jmol, one need set list textures: povray_settings["textures"] = ['jmol', 'jmol'...]
If one want change trans from 0.0 to 0.9 or higher, one need set list transmittances: povray_settings["transmittances"] = [0.9, 0.9, ...]
If one want change colors from jmol color to self-define, one need set numpy.ndarray colors: POVRAY.colors = [[0.58, 0.  , 0.58], [1., 1., 1.]...]
'''
# For function write_pov, param generic_projection_settings is param for Class PlottingVariables.
'''**generic_projection_settings:
atoms, rotation='', show_unit_cell=2, radii=None, bbox=None, 
colors=None, scale=20, maxwidth=500, extra_offset=(0., 0.)

for function write_pov, atoms and scale shouldn't be given, colors are not suggest to given.
'''


#=======================Render execute=======================
# D:\Draw\POV-Ray\core\bin\pvengine64.exe +W2556 +H2000 -Iexample.pov -OC2H161Cu144IO83.png +P +X +A0.1 +FN +C +UA

png_path = pov_ini_put.render(povray_executable='D:\\Draw\\POV-Ray\\core\\bin\\pvengine64.exe', clean_up=True)
