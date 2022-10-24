# Set name of GROMACS executable (gmx/gmx_mpi)
GMX=gmx_mpi
# Set number of threads - adjust if needed
# if using SLURM this can be set automatically using the variable SLURM_CPUS_PER_TASK
NT=10
# Set MPI launcher (srun/mpirun) - you might need to specify the number of MPI tasks
MPI=srun

# 1) energy minimization
$GMX grompp -f ../../../DATA/gromacs-mdps-charmm/0-em-steep.mdp -c ../0-TOPO/step3_input_water_ions.gro -p ../0-TOPO/topol.top
$MPI $GMX mdrun -ntomp $NT -c conf_emin.gro

# 2) 1-ns-long NPT run with positional restraints 
$GMX grompp -f ../../../DATA/gromacs-mdps-charmm/1-npt-posres.mdp -c conf_emin.gro -p ../0-TOPO/topol.top -r conf_emin.gro -n ../0-TOPO/index.ndx -maxwarn 1
$MPI $GMX mdrun -ntomp $NT -c conf_npt.gro -nsteps 500000


# 3) 1-ns-long NVT run with positional restraints
$GMX grompp -f ../../../DATA/gromacs-mdps-charmm/2-nvt-posres.mdp -c conf_npt.gro -p ../0-TOPO/topol.top -r conf_emin.gro -n ../0-TOPO/index.ndx -maxwarn 1
$MPI $GMX mdrun -ntomp $NT -c conf_nvt.gro -nsteps 500000
