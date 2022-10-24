# use Spack AMD optimzed build of gromacs@2022.2
GMX=gmx_mpi
MPI=mpirun

# in total 16 * 8 = 128 threads
NP=4   # num of mpi process
NT=8    # num of thread per mpi process

# 1) energy minimization
$GMX grompp -f ../../../DATA/gromacs-mdps-charmm/0-em-steep.mdp -c ../0-TOPO/step3_input_water_ions.gro -p ../0-TOPO/topol.top
$MPI -np $NP $GMX mdrun -ntomp $NT -c conf_emin.gro

# 2) 1-ns-long NPT run with positional restraints 
$GMX grompp -f ../../../DATA/gromacs-mdps-charmm/1-npt-posres.mdp -c conf_emin.gro -p ../0-TOPO/topol.top -r conf_emin.gro -n ../0-TOPO/index.ndx -maxwarn 2
$MPI -np $NP $GMX mdrun -ntomp $NT -c conf_npt.gro -nsteps 500000

# 3) 1-ns-long NVT run with positional restraints
$GMX grompp -f ../../../DATA/gromacs-mdps-charmm/2-nvt-posres.mdp -c conf_npt.gro -p ../0-TOPO/topol.top -r conf_emin.gro -n ../0-TOPO/index.ndx -maxwarn 1
$MPI -np $NP $GMX mdrun -ntomp $NT -c conf_nvt.gro -nsteps 500000
