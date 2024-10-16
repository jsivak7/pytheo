# general utilities


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
