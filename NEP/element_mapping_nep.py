# mapping.py
# change type according to nep.txt first line
import os
import argparse
import ast
import numpy as np
from ase.io import read
from ase.atoms import Atoms
from ase.calculators.lammps import Prism, convert
from ase.utils import writer


@writer
def write_lammps_data(
    fd,
    atoms: Atoms,
    *,
    specorder: list = None,
    force_skew: bool = False,
    prismobj: Prism = None,
    masses: bool = False,
    velocities: bool = False,
    units: str = "metal",
    atom_style: str = "atomic",
):
    """Write atomic structure data to a LAMMPS data file.

    Parameters
    ----------
    fd : file|str
        File to which the output will be written.
    atoms : Atoms
        Atoms to be written.
    specorder : list[str], optional
        Chemical symbols in the order of LAMMPS atom types, by default None
    force_skew : bool, optional
        Force to write the cell as a
        `triclinic <https://docs.lammps.org/Howto_triclinic.html>`__ box,
        by default False
    prismobj : Prism|None, optional
        Prism, by default None
    masses : bool, optional
        Whether the atomic masses are written or not, by default False
    velocities : bool, optional
        Whether the atomic velocities are written or not, by default False
    units : str, optional
        `LAMMPS units <https://docs.lammps.org/units.html>`__,
        by default "metal"
    atom_style : {"atomic", "charge", "full"}, optional
        `LAMMPS atom style <https://docs.lammps.org/atom_style.html>`__,
        by default "atomic"
    """

    # FIXME: We should add a check here that the encoding of the file object
    #        is actually ascii once the 'encoding' attribute of IOFormat objects
    #        starts functioning in implementation (currently it doesn't do
    #         anything).

    if isinstance(atoms, list):
        if len(atoms) > 1:
            raise ValueError(
                "Can only write one configuration to a lammps data file!"
            )
        atoms = atoms[0]

    fd.write("(written by ASE)\n\n")

    symbols = atoms.get_chemical_symbols()      # return a list ['O', 'H', 'H', 'O',......,'H'] len=n_atoms
    n_atoms = len(symbols)
    fd.write(f"{n_atoms} atoms\n")

    if specorder is None:
        # This way it is assured that LAMMPS atom types are always
        # assigned predictably according to the alphabetic order
        species = sorted(set(symbols))
    else:
        # To index elements in the LAMMPS data file
        # (indices must correspond to order in the potential file)
        species = specorder
    n_atom_types = len(species)
    fd.write(f"{n_atom_types} atom types\n\n")

    if prismobj is None:
        p = Prism(atoms.get_cell())
    else:
        p = prismobj

    # Get cell parameters and convert from ASE units to LAMMPS units
    xhi, yhi, zhi, xy, xz, yz = convert(p.get_lammps_prism(), "distance",
                                        "ASE", units)

    fd.write(f"0.0 {xhi:23.17g}  xlo xhi\n")
    fd.write(f"0.0 {yhi:23.17g}  ylo yhi\n")
    fd.write(f"0.0 {zhi:23.17g}  zlo zhi\n")

    if force_skew or p.is_skewed():
        fd.write(f"{xy:23.17g} {xz:23.17g} {yz:23.17g}  xy xz yz\n")
    fd.write("\n")

    if masses:
        _write_masses(fd, atoms, species, units)

    # Write (unwrapped) atomic positions.  If wrapping of atoms back into the
    # cell along periodic directions is desired, this should be done manually
    # on the Atoms object itself beforehand.
    fd.write(f"Atoms # {atom_style}\n\n")
    pos = p.vector_to_lammps(atoms.get_positions(), wrap=False)

    if atom_style == 'atomic':
        for i, r in enumerate(pos):
            # Convert position from ASE units to LAMMPS units
            r = convert(r, "distance", "ASE", units)
            s = species.index(symbols[i]) + 1       # if symbols[i] return 'O' and species is ['C', 'N', 'Fe', 'O', 'H'], then species.index(symbols[i]) return 3
            fd.write(
                "{0:>6} {1:>3} {2:23.17g} {3:23.17g} {4:23.17g}\n".format(
                    *(i + 1, s) + tuple(r)
                )
            )
    elif atom_style == 'charge':
        charges = atoms.get_initial_charges()
        for i, (q, r) in enumerate(zip(charges, pos)):
            # Convert position and charge from ASE units to LAMMPS units
            r = convert(r, "distance", "ASE", units)
            q = convert(q, "charge", "ASE", units)
            s = species.index(symbols[i]) + 1
            fd.write("{0:>6} {1:>3} {2:>5} {3:23.17g} {4:23.17g} {5:23.17g}\n"
                     .format(*(i + 1, s, q) + tuple(r)))
    elif atom_style == 'full':
        charges = atoms.get_initial_charges()
        # The label 'mol-id' has apparenlty been introduced in read earlier,
        # but so far not implemented here. Wouldn't a 'underscored' label
        # be better, i.e. 'mol_id' or 'molecule_id'?
        if atoms.has('mol-id'):
            molecules = atoms.get_array('mol-id')
            if not np.issubdtype(molecules.dtype, np.integer):
                raise TypeError((
                    "If 'atoms' object has 'mol-id' array, then"
                    " mol-id dtype must be subtype of np.integer, and"
                    " not {:s}.").format(str(molecules.dtype)))
            if (len(molecules) != len(atoms)) or (molecules.ndim != 1):
                raise TypeError((
                    "If 'atoms' object has 'mol-id' array, then"
                    " each atom must have exactly one mol-id."))
        else:
            # Assigning each atom to a distinct molecule id would seem
            # preferableabove assigning all atoms to a single molecule
            # id per default, as done within ase <= v 3.19.1. I.e.,
            # molecules = np.arange(start=1, stop=len(atoms)+1,
            # step=1, dtype=int) However, according to LAMMPS default
            # behavior,
            molecules = np.zeros(len(atoms), dtype=int)
            # which is what happens if one creates new atoms within LAMMPS
            # without explicitly taking care of the molecule id.
            # Quote from docs at https://lammps.sandia.gov/doc/read_data.html:
            #    The molecule ID is a 2nd identifier attached to an atom.
            #    Normally, it is a number from 1 to N, identifying which
            #    molecule the atom belongs to. It can be 0 if it is a
            #    non-bonded atom or if you don't care to keep track of molecule
            #    assignments.

        for i, (m, q, r) in enumerate(zip(molecules, charges, pos)):
            # Convert position and charge from ASE units to LAMMPS units
            r = convert(r, "distance", "ASE", units)
            q = convert(q, "charge", "ASE", units)
            s = species.index(symbols[i]) + 1
            fd.write("{0:>6} {1:>3} {2:>3} {3:>5} {4:23.17g} {5:23.17g} "
                     "{6:23.17g}\n".format(*(i + 1, m, s, q) + tuple(r)))
    else:
        raise NotImplementedError

    if velocities and atoms.get_velocities() is not None:
        fd.write("\n\nVelocities\n\n")
        vel = p.vector_to_lammps(atoms.get_velocities())
        for i, v in enumerate(vel):
            # Convert velocity from ASE units to LAMMPS units
            v = convert(v, "velocity", "ASE", units)
            fd.write(
                "{0:>6} {1:23.17g} {2:23.17g} {3:23.17g}\n".format(
                    *(i + 1,) + tuple(v)
                )
            )

    fd.flush()

def _write_masses(fd, atoms: Atoms, species: list, units: str):
    symbols_indices = atoms.symbols.indices()       # return a dict {'element_symbol_1': array([a,b,c,...]), 'element_symbol_2': array([x,y,z,...])}
    fd.write("Masses\n\n")
    for i, s in enumerate(species):
        # Find the first atom of the element `s` and extract its mass
        if s not in symbols_indices:
            # Still write mass if the system does not contain the element `s`.
            # Cover by `float` to make a new object for safety
            mass = float(Atoms(s)[0].mass)
            # Convert mass from ASE units to LAMMPS units
            mass = convert(mass, "mass", "ASE", units)
            atom_type = i + 1
            fd.write(f"{atom_type:>6} {mass:23.17g} # {s}\n")
        else:
            # Cover by `float` to make a new object for safety
            mass = float(atoms[symbols_indices[s][0]].mass)
            # Convert mass from ASE units to LAMMPS units
            mass = convert(mass, "mass", "ASE", units)
            atom_type = i + 1
            fd.write(f"{atom_type:>6} {mass:23.17g} # {s}\n")
    fd.write("\n")

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('path_nep_txt', help='NEP model path')
    parser.add_argument('path_old_data', help='old data file to be mapping')
    parser.add_argument('old_data_type', help='old data file type-element map, format: dict, exp: "{1:1,2:8}"')
    args = parser.parse_args()

    if args.path_old_data:
        path_data = args.path_old_data
    else:
        path_data = '/home/yuqinghan/Project/Fe-C-N4/lmp/frozen_model/data-structure/water/water_60.data'
    if args.path_nep_txt:
        path_nep_txt = args.path_nep_txt
    else:
        path_nep_txt = '/home/yuqinghan/Project/DAC/NEP-infer/nep.txt'

    if args.old_data_type and isinstance(ast.literal_eval(args.old_data_type), dict):
        old_data_type = ast.literal_eval(args.old_data_type)
    else:
        print(f'======Warning: use default old_data_type {old_data_type}!!!======')
        old_data_type = {1:1, 2:8}

    # get nep_type_mapping_list
    with open(path_nep_txt, 'r') as file:
        line_1 = file.readline()
    line_1_parts = line_1.split()
    nep_version = line_1_parts[0]; type_len = int(line_1_parts[1]); nep_type_mapping_list = line_1_parts[2:]
    if type_len != len(nep_type_mapping_list):
        raise Exception(f'''type length({type_len}) not equal to length of types({nep_type_mapping_list})''')


    atoms =  read(path_data, format='lammps-data', Z_of_type=old_data_type, style='atomic', sort_by_id=True, units='metal')

    pre_dir = os.path.dirname(path_data)
    old_filename = path_data.split('/')[-1].split('.')[-2]
    new_filename = old_filename + '_nep_fix.data'
    new_file_path = os.path.join(pre_dir, new_filename)

    write_lammps_data(new_file_path, atoms, masses=True, units='metal', atom_style='atomic',  specorder=nep_type_mapping_list)


if __name__ == "__main__":
    main()
