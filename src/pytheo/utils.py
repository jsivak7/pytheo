# general utilities

from ase import Atoms


def read(file_path: str):
    """Read in structure file to ASE Atoms Object

    Args:
        file_path (str): relative path to structure file

    Returns:
        Atoms: modular object to perform other operations
    """
    from ase import io

    s = io.read(f"{file_path}")
    return s


def write(struc: Atoms, file_path: str):
    """Write ASE Atoms object to structure file
    NOTE that always writes "direct" coordinates

    Args:
        struc (Atoms): structure to be written
        file_path (str): relative path to write structure file, include suffix for desired file type (*.vasp, *.cif, etc.)

    Raises:
        FileExistsError: if given path for output file already exists, ensures that previously generated SQSs are not overwritten
    """
    import os
    from ase import io

    if os.path.exists(file_path):
        raise FileExistsError(file_path)

    io.write(f"{file_path}", struc, direct=True)


def rattle_atoms(struc: Atoms, stddev=0.02):
    """Rattles atoms of a given input structure

    Args:
        struc (Atoms): structure to be rattled
        stddev (float, optional): standard deviation for amount of rattling to perform in Angstroms. Defaults to 0.02.

    Returns:
        Atoms: ASE Atoms object for rattled structure
    """
    import random

    struc.rattle(stddev, seed=int(random.uniform(0, 2000)))  # random seed

    return struc
