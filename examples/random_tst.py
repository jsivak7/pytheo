from pytheo.structure import read, write, decorate_randomly, make_supercell

unitcell = read("MgO.poscar")

random = decorate_randomly(
    unitcell,
    (10, 5, 5),
    chemical_symbols=[
        ["Mg", "Ni", "Co", "Zn"],
        ["Mg", "Ni", "Co", "Zn"],
        ["Mg", "Ni", "Co", "Zn"],
        ["Mg", "Ni", "Co", "Zn"],
        ["O"],
        ["O"],
        ["O"],
        ["O"],
    ],
    concentrations={"Mg": 0.25, "Ni": 0.25, "Co": 0.25, "Zn": 0.25},
)

write(random, "random.vasp")
