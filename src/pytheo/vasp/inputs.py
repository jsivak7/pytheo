"""
For making files to run DFT calculations

I have opted to use the `KSPACING` INCAR flag rather than the `KPOINTS` file as this is better for high-throughput calculations of structures with different shapes and sizes

It is assumed that you are using the CUSTODIAN job management and error handler package for VASP calculations
    - this greatly assists in minimizing human input time for calculations, especially for high-throughput
    - https://github.com/materialsproject/custodian


"""

from ase import Atoms


def make_relax(
    struc: Atoms,
    output_path: str,
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

    if functional == "r2scan":
        calc = MPScanRelaxSet(
            structure=s,
            user_incar_settings=incar_settings,
            user_potcar_functional="PBE_54",
            sort_structure=True,
        )
    elif functional == "pbe":
        calc = MPRelaxSet(
            structure=s,
            user_incar_settings=incar_settings,
            user_potcar_functional="PBE",
            sort_structure=True,
        )
    else:
        raise ValueError(f"{functional} is not available.")

    calc.write_input(f"{output_path}")


# I feel that this is overly complicated, but will leave for now and adjust in the future...
def write_psu_submission(
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
        output_path (str): relative path to write submission file
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
#SBATCH --partition=open
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


def write_custodian_relax(output_path: str, kspacing=0.25, half_kmesh_first_relax=True):
    """Write a generic double relaxation script for calculation workflow and error handling using Custodian (https://github.com/materialsproject/custodian).

    Args:
        output_path (str): relative path to write submission file
        kspacing (float, optional): K-point mesh spacing with VASP KSPACING tag (https://www.vasp.at/wiki/index.php/KSPACING). Defaults to 0.25.
        half_kmesh_first_relax (bool, optional): Use more sparse k-mesh for initial relax. Defaults to True.
    """
    if half_kmesh_first_relax == True:
        cstdn_script = f"""kspacing_initial = {kspacing*2}\nkspacing = {kspacing}\n\n"""
    else:
        cstdn_script = f"""kspacing_initial = {kspacing}\n,kspacing = {kspacing}\n\n"""

    cstdn_script += """import os
from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler
from custodian.vasp.jobs import VaspJob

subset = list(VaspErrorHandler.error_msgs.keys())
handlers = [VaspErrorHandler(errors_subset_to_catch=subset)]

step1 = VaspJob(
    vasp_cmd=["srun", "vasp_std"],
    final=False,
    suffix=".1",
    settings_override=[
        {
            "dict": "INCAR",
            "action": {
                "_set": {
                    "KSPACING": kspacing_initial
                }
            },
        },
    ],
)

step2 = VaspJob(
    vasp_cmd=["srun", "vasp_std"],
    final=True,
    settings_override=[
        {
            "dict": "INCAR",
            "action": {
                "_set": {
                    "KSPACING": kspacing
                }
            },
        },
        {"file": "CONTCAR", "action": {"_file_copy": {"dest": "POSCAR"}}},
    ],
)


jobs = [step1, step2]
c = Custodian(handlers, jobs, max_errors=3)
c.run()"""
    with open(f"{output_path}/cstdn.py", "w+") as f:
        f.writelines(cstdn_script)


def set_up_dos_calc():
    """
    Call in the same location as a relaxation
    Moves all relaxation files to a new directory called '01_relax' and gets everything ready for a density of states calculation.
    """

    import os
    from pymatgen.io.vasp.outputs import Eigenval, Vasprun

    os.mkdir("01_relax")
    os.system("mv * 01_relax")
    os.system(
        "cp 01_relax/CONTCAR 01_relax/INCAR 01_relax/WAVECAR 01_relax/CHGCAR 01_relax/POTCAR 01_relax/submitvasp ."
    )

    os.system("mv CONTCAR POSCAR")
    os.system(
        "perl -pi -e 's/python cstdn.py/srun vasp_std/g' submitvasp"
    )  # since not using custodian here

    e = Eigenval("01_relax/EIGENVAL")
    num_elec = e.nelect
    new_num_elec = num_elec - 1 + 0.999999
    os.system(f"echo NELECT = {new_num_elec} >> INCAR")
    print(f"number of electrons: {num_elec} -> {new_num_elec}")

    v = Vasprun("01_relax/vasprun.xml")
    efermi = v.efermi
    print(f"fermi energy = {efermi} eV")
    emin = efermi - 8  # for a reasonable energy window
    emax = efermi + 6  # for a reasonable energy window
    os.system(f"echo EMIN = {emin} >> INCAR")
    os.system(f"echo EMAX = {emax} >> INCAR")

    os.system("perl -pi -e 's/ISMEAR = 0/ISMEAR = -5/g' INCAR")
    os.system("perl -pi -e 's/NSW = 250/NSW = 0/g' INCAR")
    os.system("perl -pi -e 's/LCHARG = True/LCHARG = False/g' INCAR")
    os.system("perl -pi -e 's/LWAVE = True/LWAVE = False/g' INCAR")

    os.system(f"echo NEDOS = 501 >> INCAR")


# TODO just copied in for now from the perovskite HEO project - still needs to be fixed
def set_up_optical_calc(sigma=0.1):
    """
    For making a new directory for optical calculations (should be called in directory where "02_static" exists).
    - NBANDS is doubled from the default value of 02_static calculation

    Args:
        sigma (float, optional): Smearing parameter than is usually increased for optical calculations compared to relaxation. Defaults to 0.2.
    """
    import os
    from pymatgen.io.vasp.outputs import Eigenval, Vasprun

    print(f"-> Making optical directory...")

    os.mkdir("04_optical")
    os.system(
        f"cp 02_static/INCAR 02_static/POSCAR 02_static/KPOINTS 02_static/POTCAR 02_static/CHGCAR 02_static/WAVECAR 02_static/submitvasp ./04_optical"
    )
    os.chdir("04_optical")

    e = Eigenval("../02_static/EIGENVAL")

    nbands = e.nbands
    new_nbands = nbands * 2  # doubling NBANDS
    os.system(f"echo NBANDS = {new_nbands} >> INCAR")
    print(f"number of bands: {nbands} -> {new_nbands}")

    os.system("perl -pi -e 's/LWAVE = True/LWAVE = False/g' INCAR")
    os.system("perl -pi -e 's/LCHARG = True/LCHARG = False/g' INCAR")
    os.system(f"echo NEDOS = 20000 >> INCAR")  # for adequate sampling
    os.system(f"perl -pi -e 's/SIGMA = 0.05/SIGMA = {sigma}/g' INCAR")
    os.system("perl -pi -e 's/ALGO = All/ALGO = Normal/g' INCAR")
    os.system("perl -pi -e 's/EDIFF = 1e-06/EDIFF = 1e-08/g' INCAR")
    os.system("perl -pi -e 's/LREAL = Auto/LREAL = False/g' INCAR")
    os.system(f"echo LOPTICS = True >> INCAR")
    os.system(f"echo CSHIFT = 0.01 >> INCAR")


# TODO just copied in for now from the perovskite HEO project - still needs to be fixed
def set_up_bandstructure_calc():
    """For making a new directory for density of states calculations (should be called in directory where "02_static" exists)
    - NBANDS is 1.5x from the default value of 02_static calculation

    Additional, manual steps are required for this function still due to band structure process:
        NOTE that this is for the 'regular' band structure calculations (i.e., not unfolding)
        1. run 01_pbe calculation
        2. get high-symmetry k-path with sumo (> sumo-kgen --hybrid --symprec 0.1 --pymatgen)
        3. copy KPOINTS_band to KPOINTS
        4. copy WAVECAR from finished 01_pbe calculation to the 'bandstructure' directory
        5. now can run actual metaGGA band structure calculation
    """

    import os
    from pymatgen.io.vasp.outputs import Eigenval

    print(f"-> Making bandstructure directory...")

    os.mkdir("05_bandstructure")
    os.system(
        f"cp 02_static/INCAR 02_static/IBZKPT 02_static/POSCAR 02_static/KPOINTS 02_static/POTCAR 02_static/CHGCAR 02_static/WAVECAR 02_static/submitvasp ./05_bandstructure"
    )
    os.chdir("05_bandstructure")

    e = Eigenval("../02_static/EIGENVAL")

    nbands = e.nbands
    new_nbands = int(nbands * 1.5)  # increasing NBANDS by 1.5x
    os.system(f"echo NBANDS = {new_nbands} >> INCAR")
    print(f"number of bands: {nbands} -> {new_nbands}")

    os.system("perl -pi -e 's/LCHARG = True/LCHARG = False/g' INCAR")
    os.system(f"echo NEDOS = 5001 >> INCAR")  # for adequate sampling
    os.system("perl -pi -e 's/LREAL = Auto/LREAL = False/g' INCAR")

    os.mkdir("01_pbe")
    os.system("cp * 01_pbe")
    os.system("rm CHGCAR WAVECAR")
    os.chdir("01_pbe")

    os.system("perl -pi -e 's/METAGGA = R2scan/GGA = PE/g' INCAR")
    os.system(f"echo ICHARG = 11 >> INCAR")
