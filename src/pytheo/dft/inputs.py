# for making DFT calculation input files

from ase import Atoms


def make_relax(
    struc: Atoms,
    output_path: str,
    kpt_mesh: tuple,
    user_incar_changes=None,
    functional="r2scan",
):
    """Function for making VASP input files (INCAR, KPOINTS, POSCAR, POTCAR) using Pymatgen for a relaxation calculation

    Args:
        struc (Atoms): structure to be relaxed
        output_path (str): relative path to output generated files
        kpt_mesh (tuple): k-point mesh to use for relaxation, gamma-centered mesh only
        user_incar_changes (dict, optional): changes that deviate from default INCAR parameters given in *.yaml files. Defaults to None.
        functional (str, optional): XC functional to be used. Defaults to "r2scan".

    Raises:
        FileExistsError: if output_path already exists
        ValueError: if an invalid functional is given
    """

    import os
    import yaml
    from pymatgen.core import Structure
    from pymatgen.io.vasp.sets import MPRelaxSet, MPScanRelaxSet
    from pymatgen.io.vasp.inputs import Kpoints

    if os.path.exists(output_path):
        raise FileExistsError(output_path)
    else:
        os.mkdir(output_path)

    module_dir = os.path.dirname(__file__)
    s = Structure.from_ase_atoms(struc)

    with open(f"{module_dir}/incar_defaults/{functional}.yaml", "r") as f:
        incar_settings = yaml.load(f, Loader=yaml.SafeLoader)

    if (
        bool(user_incar_changes) != False
    ):  # only update INCAR settings if user asks for changes
        incar_settings.update(user_incar_changes)

    kpoints = Kpoints.from_dict(
        {
            "comment": "KPOINTS",
            "generation_style": "Gamma",
            "kpoints": [kpt_mesh],
        }
    )

    if functional == "pbe":
        calc = MPRelaxSet(
            structure=s,
            user_incar_settings=incar_settings,
            user_kpoints_settings=kpoints,
            user_potcar_functional="PBE",
            sort_structure=True,
        )
    elif functional == "r2scan":
        calc = MPScanRelaxSet(
            structure=s,
            user_incar_settings=incar_settings,
            user_kpoints_settings=kpoints,
            user_potcar_functional="PBE_54",
            sort_structure=True,
        )
    else:
        raise ValueError(f"{functional} is not available.")

    calc.write_input(f"{output_path}")


# it assumed that you are using custodian here, which I highly recommend for error handling
# I feel that this is overly complicated, but will leave for now and adjust in the future...
def write_psu_roar_collab_submission(
    name: str,
    output_path: str,
    type: str,
    nodes=1,
    cpu=64,
    mem_per_cpu="3500MB",
    hours=120,
    partition="sla-prio",
    account="ixd4_n_bc",
):
    """Slurm submission script writer for Roar Collab VASP jobs

    Args:
        name (str): job name
        output_path (_type_): relative path to write submission file
        type (str): type of calculation for writing correct slurm submission ["basic", "standard", "open"]
        nodes (int, optional): number of nodes. Defaults to 1.
        cpu (int, optional): number of cpu per node. Defaults to 64.
        mem_per_cpu (str, optional): amount of memory per cpu. Defaults to "3500MB".
        hours (int, optional): number of hours. Defaults to 120.
        partition (str, optional): SLURM submission partition. Defaults to "sla-prio".
        account (str, optional): Paid allocations. Defaults to "ixd4_n_bc".
    """
    if type == "basic":

        submitvasp = f"""#!/bin/bash
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={cpu}
#SBATCH --mem-per-cpu={mem_per_cpu}
#SBATCH --time={hours}:00:00
#SBATCH --partition={partition}
#SBATCH --qos=burst4x
#SBATCH --account={account}
#SBATCH --output=vasp.out
#SBATCH --error=vasp.err
#SBATCH --job-name={name}

export UCX_TLS=all

cd $SLURM_SUBMIT_DIR
module purge
module use /storage/icds/RISE/sw8/modules/
module load vasp/vasp-6.4.1v

eval "$(conda shell.bash hook)"
conda activate research

python cstdn.py"""

    elif type == "standard":

        submitvasp = f"""#!/bin/bash
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={cpu}
#SBATCH --mem-per-cpu={mem_per_cpu}
#SBATCH --time={hours}:00:00
#SBATCH --partition={partition}
#SBATCH --qos=burst4x
#SBATCH --account={account}
#SBATCH --output=vasp.out
#SBATCH --error=vasp.err
#SBATCH --job-name={name}
#SBATCH --constraint=icelake

export UCX_TLS=all

cd $SLURM_SUBMIT_DIR
module purge
module use /storage/icds/RISE/sw8/modules/
module load vasp/vasp-6.4.1v

eval "$(conda shell.bash hook)"
conda activate research

python cstdn.py"""

    elif type == "open":

        submitvasp = f"""#!/bin/bash
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={cpu}
#SBATCH --mem-per-cpu={mem_per_cpu}
#SBATCH --time={hours}:00:00
#SBATCH --partition={partition}
#SBATCH --output=vasp.out
#SBATCH --error=vasp.err
#SBATCH --job-name={name}
#SBATCH --constraint=icelake

export UCX_TLS=all

cd $SLURM_SUBMIT_DIR
module purge
module use /storage/icds/RISE/sw8/modules/
module load vasp/vasp-6.4.1v

eval "$(conda shell.bash hook)"
conda activate research

python cstdn.py"""

    with open(f"{output_path}/submitvasp", "w") as write_runvasp:
        write_runvasp.write(submitvasp)
