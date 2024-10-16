# may in the future change this to be a class
# only have functions allows for the input to be retained however since the object itself is not operated on...

from ase import Atoms


def read(file_path: str):
    """Read in structure file to ASE Atoms Object

    Parameters
    ----------
    file_path : str
        relative path to structure file

    Returns
    -------
    ASE Atoms object
        modular format to perform other operations on
    """
    from ase import io

    s = io.read(f"{file_path}")
    return s


def write(struc: Atoms, file_path: str):
    """Write ASE Atoms object to structure file
    NOTE that always writes "direct" coordinates

    Parameters
    ----------
    struc: Atoms
        structure to be written
    output_path : str
        relative path to write structure file, include suffix for desired file type (*.vasp, *.cif, etc.)

    Raises
    ------
    FileExistsError
        if given path for output file already exists
        ensures that previously generated SQSs are not overwritten
    """
    import os
    from ase import io

    if os.path.exists(file_path):
        raise FileExistsError(file_path)

    io.write(f"{file_path}", struc, direct=True)


def make_supercell(struc: Atoms, dimensions: tuple):
    """Make supercell from ASE Atoms object

    Parameters
    ----------
    struc : Atoms
        structure to be made into a supercell, usually a unit cell
    dimensions : tuple
        (x, y, z) multipliers for supercell generation

    Returns
    -------
    Atoms
        ASE Atoms object for supercell
    """
    return struc.repeat(dimensions)


def rattle_atoms(struc: Atoms, stddev=0.02):
    """Rattles atoms of a given input structure

    Parameters
    ----------
    struc : Atoms
        structure to be rattled
    stddev : float, optional
        standard deviation for amount of rattling to perform in Angstroms, by default 0.02

    Returns
    -------
    Atoms
        ASE Atoms object for rattled structure
    """
    import random

    struc.rattle(stddev, seed=int(random.uniform(0, 2000)))  # random seed

    return struc


def generate_sqs(
    struc: Atoms,
    dimensions: tuple,
    chemical_symbols: list,
    cutoffs: list,
    concentrations: dict,
    num_steps: int,
):
    """Generate special quasirandom structure (SQS) for an arbitrary input structure using ICET (https://icet.materialsmodeling.org/)

    Parameters
    ----------
    struc : Atoms
        structure to be made into an SQS, usually a unit cell
    dimensions : tuple
        (x, y, z) multipliers for supercell generation
    chemical_symbols : list
        list of lists for allowed elements following same order as unitcell structure
        example: [["Sr"], ["Ti", "Cr"], ["O"], ["O"], ["O"]] for perovskite
    cutoffs : list
        cutoffs in order of multiplicity (pair, triplet, quadruplet, etc.)
    concentrations : dict
        fractions of different elements for each lattice site
        only need to specify those that are not 1.0
    num_steps : int
        number of Monte Carlo steps to run the SQS generation

    Returns
    -------
    Atoms
        ASE Atoms object for SQS
    """

    import os
    from icet import ClusterSpace
    from icet.tools.structure_generation import (
        generate_sqs_from_supercells,
        _get_sqs_cluster_vector,
    )
    from icet.input_output.logging_tools import set_log_config

    set_log_config(level="INFO")

    supercell = [struc.repeat(dimensions)]

    cs = ClusterSpace(struc, cutoffs, chemical_symbols)
    print(cs)

    sqs = generate_sqs_from_supercells(
        cluster_space=cs,
        supercells=supercell,
        target_concentrations=concentrations,
        n_steps=num_steps,
    )

    trial_cluster_vector = cs.get_cluster_vector(sqs)
    perfectly_random_cluster_vector = _get_sqs_cluster_vector(
        cluster_space=cs, target_concentrations=concentrations
    )

    print(f"\nTrial Cluster Vector ->\n{trial_cluster_vector}")
    print(f"\nPerfectly Random Cluster Vector ->\n{perfectly_random_cluster_vector}")

    return sqs


def decorate_randomly(
    struc: Atoms,
    dimensions: tuple,
    chemical_symbols: list,
    concentrations: dict,
):
    """Randomly decorate an arbitrary input structure with ICET (https://icet.materialsmodeling.org/)

    Parameters
    ----------
    struc : Atoms
        structure to be made into an randomly decorated structure, usually a unit cell
    dimensions : tuple
        (x, y, z) multipliers for supercell generation
    chemical_symbols : list
        list of lists for allowed elements following same order as unitcell structure
        example: [["Sr"], ["Ti", "Cr"], ["O"], ["O"], ["O"]] for perovskite
    concentrations : dict
        fractions of different elements for each lattice site
        only need to specify those that are not 1.0

    Returns
    -------
    Atoms
        ASE Atoms object for randomly decorated structure
    """

    import os
    from icet import ClusterSpace
    from icet.tools.structure_generation import occupy_structure_randomly

    cs = ClusterSpace(
        struc, [0], chemical_symbols
    )  # cutoffs do not matter for random decoration, but needed for cluster space

    supercell = struc.repeat(dimensions)
    occupy_structure_randomly(supercell, cs, concentrations)
    random = supercell

    return random
