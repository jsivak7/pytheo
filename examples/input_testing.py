from pytheo import utils
from pytheo.dft import inputs

s = utils.read("MgO.poscar")
inputs.make_relax(s, "MgO_unitcell")
inputs.write_psu_roar_collab_submission("MgO", "MgO_unitcell", type="open")
inputs.write_custodian_relax("MgO_unitcell")
