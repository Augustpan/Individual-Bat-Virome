from pymol import cmd

with open("../../../modified_scripts/interface.pdb") as f:
    lines = f.readlines()

sele = []
for line in lines:
    line = line.strip().split()
    chain = line[4]
    resi = line[5]
    sele.append(f"(c. {chain} & resi {resi})")
sele = "|".join(sele)

cmd.load("conf_slt.pdb")

cmd.select("c. A")
cmd.set_name("sele", "rmsd-")


if cmd.select("c. A") > cmd.select("c. B"):
    A = "ACE"
    B = "RBD"
else:
    A = "RBD"
    B = "ACE"

cmd.select("c. A & bb.")
cmd.set_name("sele", "rmsd-{}".format(A))
cmd.select("c. B & bb.")
cmd.set_name("sele", "rmsd-{}".format(B))

cmd.select(sele)
cmd.set_name("sele", "rmsd-INTER-all")
cmd.select("rmsd-INTER-all & bb.")
cmd.set_name("sele", "rmsd-INTER")

cmd.save("PLUMED/rmsd-ACE.pdb", "rmsd-ACE", -1, "pdb")
cmd.save("PLUMED/rmsd-RBD.pdb", "rmsd-RBD", -1, "pdb")
cmd.save("PLUMED/rmsd-INTER-raw.pdb", "rmsd-INTER", -1, "pdb")