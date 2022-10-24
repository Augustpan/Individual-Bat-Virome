# Set name of GROMACS executable (gmx/gmx_mpi)
GMX=gmx_mpi
# Set number of threads - adjust if needed
# if using SLURM this can be set automatically using the variable SLURM_CPUS_PER_TASK
NT=10
# Set MPI launcher (srun/mpirun) - you might need to specify the number of MPI tasks
MPI=srun

# 1) creation of tpr file from equilibrated conformation - run once before first production run
$GMX grompp -f ../../../DATA/gromacs-mdps-charmm/3-nvt-production.mdp -c ../1-EQUIL/conf_nvt.gro -n ../0-TOPO/index.ndx -p ../0-TOPO/topol.top

# 2) production in chunks of 24h (adjust based on available resources) - resubmit the line below to restart the simulation from the state (.cpt) file
$MPI $GMX mdrun -ntomp $NT -maxh 24.0 -cpi -noappend
