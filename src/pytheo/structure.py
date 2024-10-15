# might should make all of the inputs and outputs only as ASE Atoms objects and/or pymatgen structures objects for better clarity?? Then I won't need to worry about the output names and all of this...


def write_to_cif(input_struc_path: str, output_struc_path: str):
    """Given an input structure file, outputs a *.cif file format

    Parameters
    ----------
    input_struc_path : str
        relative path to input structure
    output_struc_path : str
        relative path to output structure as *.cif

    Returns
    -------
    _type_
        _description_
    """
    from pymatgen.core.structure import Structure

    print(f"========== writing .cif file to {output_struc_path} ==========")

    s = Structure.from_file(struc_path)
    s.to_file(output_struc_path)
    return None


def rattle_atoms(file_path: str, stddev=0.02, write_pristine=True):
    """Rattles atoms of a given input structure

    Parameters
    ----------
    file_path : str
        relative file path to structure file
    stddev : float, optional
        standard deviation for amount of rattling to perform in Angstroms, by default 0.02
    write_pristine : bool, optional
        retain pristine (unrattled) structure by writing *_pristine, by default True

    Returns
    -------
    ASE Atoms object
        rattled structure that was written
    """

    from ase.io import read
    from ase.io.vasp import write_vasp
    import random

    s = read(f"{file_path}")
    write_vasp(f"{file_path}_pristine", s, direct=True)
    s.rattle(
        stddev, seed=int(random.uniform(0, 2000))
    )  # random seed is used to ensure different distortions are applied
    write_vasp(f"{file_path}", s, direct=True)

    return s


def make_sqs(
    unitcell_path: str,
    output_path: str,
    output_name: str,
    chemical_symbols: list,
    cutoffs: list,
    supercell_dimensions: tuple,
    target_concentrations: dict,
    num_mc_steps: int,
):
    """Generate SQS for an arbitrary input structure using ICET (https://icet.materialsmodeling.org/)

    Parameters
    ----------
    unitcell_path : str
        relative path to unit cell with which SQS will be generated from, can be any ASE file format
    output_path : str
        relative path to write generated SQS as a *.vasp
    output_name : str
        name of file generated SQS will be written to, only a *.vasp file currently
    chemical_symbols : list
        list of lists for allowed elements following same order as unitcell structure
        example: [["Sr"], ["Ti", "Cr"], ["O"], ["O"], ["O"]]
    cutoffs : list
        cutoffs in order of multiciplity (pair, triplet, quadruplet, etc.)
    supercell_dimensions : tuple
        x, y, and z multipliers for supercell
    target_concentrations : dict
        fractions of different elements for each lattice site
        only need to specify those that are not 1.0
    num_mc_steps : int
        number of Monte Carlo steps to run the SQS generation

    Returns
    -------
    ASE Atoms Object
        generated SQS that was written to a *.vasp file

    Raises
    ------
    FileExistsError
        if given path for output file already exists
        ensures that previously generated SQSs are not overwritten
    """

    import os
    from ase import Atom
    from ase.build import bulk
    from ase.io import write, read
    from icet import ClusterSpace
    from icet.tools.structure_generation import (
        generate_sqs,
        generate_sqs_from_supercells,
        _get_sqs_cluster_vector,
    )
    from icet.input_output.logging_tools import set_log_config
    import os

    set_log_config(level="INFO")

    if os.path.exists(f"{output_path}/{output_name}.vasp"):
        raise FileExistsError(f"{output_path}/{output_name}.vasp")

    unitcell = read(unitcell_path)

    cs = ClusterSpace(unitcell, cutoffs, chemical_symbols)
    print(f"\n{cs}\n")

    supercell = [unitcell.repeat(supercell_dimensions)]

    sqs = generate_sqs_from_supercells(
        cluster_space=cs,
        supercells=supercell,
        target_concentrations=target_concentrations,
        n_steps=num_mc_steps,
    )

    write(f"{output_path}/{output_name}.vasp", sqs)

    trial_cluster_vector = cs.get_cluster_vector(sqs)
    perfectly_random_cluster_vector = _get_sqs_cluster_vector(
        cluster_space=cs, target_concentrations=target_concentrations
    )

    print("\nTrial Cluster Vector -> \n" + str(trial_cluster_vector))
    print(
        "\n\nPerfectly Random Cluster Vector -> \n"
        + str(perfectly_random_cluster_vector)
    )

    return sqs


def make_random(
    unitcell_path: str,
    output_path: str,
    output_name: str,
    chemical_symbols: list,
    supercell_dimensions: tuple,
    concentrations: dict,
):
    """Generate randomly decorated structure for an arbitrary input structure using ICET (https://icet.materialsmodeling.org/)

    Parameters
    ----------
    unitcell_path : str
        relative path to unit cell with which SQS will be generated from, can be any ASE file format
    output_path : str
        relative path to write generated SQS as a *.vasp
    output_name : str
        name of file generated SQS will be written to, only a *.vasp file currently
    chemical_symbols : list
        list of lists for allowed elements following same order as unitcell structure
        example: [["Sr"], ["Ti", "Cr"], ["O"], ["O"], ["O"]]
    supercell_dimensions : tuple
        x, y, and z multipliers for supercell
    concentrations : dict
        fractions of different elements for each lattice site

    Returns
    -------
    ASE Atoms Object
        generated randomly decorated structure that was written to a *.vasp file

    Raises
    ------
    FileExistsError
        if given path for output file already exists
        ensures that previously generated structures are not overwritten
    """

    import os
    from ase import Atom
    from ase.build import bulk
    from ase.io import write, read
    from icet import ClusterSpace
    from icet.tools.structure_generation import occupy_structure_randomly
    from icet.input_output.logging_tools import set_log_config

    set_log_config(level="INFO")

    if os.path.exists(f"{output_path}/{output_name}.vasp"):
        raise FileExistsError(
            f"This file already exists ({output_path}/{output_name}.vasp)"
        )

    unitcell = read(unitcell_path)

    cs = ClusterSpace(unitcell, [0], chemical_symbols)

    supercell = unitcell.repeat(supercell_dimensions)

    occupy_structure_randomly(supercell, cs, concentrations)

    write(f"{output_path}/{output_name}.vasp", supercell)

    return supercell
