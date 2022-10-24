# use Spack AMD optimzed build of gromacs@2022.2
GMX=/home1/panyuanfei/spack/opt/spack/linux-ubuntu18.04-zen2/aocc-3.2.0/gromacs-2022.2-s4i42f44xw6id5ovhs7dyfrtw3rl6ff2/bin/gmx_mpi
MPI=/home1/panyuanfei/spack/opt/spack/linux-ubuntu18.04-zen2/aocc-3.2.0/openmpi-4.1.1-mgk2k5fnbudpsnzxqf5jxj6dxkeajpuq/bin/mpirun

# in total 16 * 8 = 128 threads
NP=16   # num of mpi process
NT=8    # num of thread per mpi process

# 1) creation of tpr file from equilibrated conformation - run once before first production run
$GMX grompp -f ../../../DATA/gromacs-mdps-charmm/3-nvt-production.mdp -c ../1-EQUIL/conf_nvt.gro -n ../0-TOPO/index.ndx -p ../0-TOPO/topol.top

# 2) production in chunks of 24h (adjust based on available resources) - resubmit the line below to restart the simulation from the state (.cpt) file
$MPI -np $NP $GMX mdrun -ntomp $NT -maxh 12.0 -cpi -noappend
