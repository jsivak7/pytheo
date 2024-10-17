from pytheo.dft.inputs import make_relax, write_psu_roar_collab_submission
from pytheo.utils import read

s = read("MgO.poscar")
print(s)
make_relax(s, "relax", (4, 4, 4), functional="pbe")
write_psu_roar_collab_submission("MgO", "relax", type="basic")
