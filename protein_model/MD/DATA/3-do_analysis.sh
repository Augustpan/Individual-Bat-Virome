# Set name of GROMACS executable (gmx/gmx_mpi)
GMX=gmx_mpi
GMX="/usr/local/gromacs/bin/gmx"
# Set MPI launcher (srun/mpirun) - you might need to specify the number of MPI tasks
MPI=srun
# Set FoldX executable
FOLDX=~/bin/FoldX/foldx
FOLDX="/usr/local/bin/foldx_20221231"
# Set ROSETTA executable
ROSETTA=~/rosetta/3.11/rosetta_bin_linux_2019.35.60890_bundle/main/source/bin/InterfaceAnalyzer.static.linuxgccrelease

###
### preparing trajectory
###
# 1) first we need to cat together all the production trajectories
$GMX trjcat -f ../2-PRODUCTION/traj_comp.part*.xtc -o traj_all.xtc -cat
# 2) now we fix PBCs - print every 10 steps (only protein atoms). This will result in 10001 frames
echo 1 | $GMX trjconv -f traj_all.xtc -o traj_all_PBC.xtc -s ../1-EQUIL/topol.tpr -pbc nojump -skip 10 -n ../0-TOPO/index.ndx -e 1000000
# 3) remove temporary trajectory
rm traj_all.xtc

###
### now we do the analysis on the full trajectory
###
# 4) PLUMED post-processing: RMSD calculation and salt bridge analysis
cd PLUMED
plumed driver --plumed plumed_driver.dat --mf_xtc ../traj_all_PBC.xtc
# postprocess COLVAR_SB file to calculate frequency of salt bridge formation
bash ../../../../DATA/calculate_sb.sh
# done
cd ../

# 5) Hbonds analysis
cd HBONDS
# execute python script: requires python3 and MDAnalysis
python ../../../../DATA/count-HB.py topol_slt.tpr ../traj_all_PBC.xtc  
# done
cd ../

# 6) Create individual PDBs from traj_all_PBC.xtc
# you might want to parallelize this using a job-array
mkdir PDBs; cd PDBs
# execute python script: requires python3 and MDTraj
python ../../../../DATA/get_PDB_frames.py ../conf_slt.pdb ../traj_all_PBC.xtc 0 10001
# done
cd ../

# 7) FoldX scoring
mkdir FOLDX; cd FOLDX
# you might want to parallelize this loop using a job-array
for j in `seq 0 10000`
do
 # padding for PDB filename
 jj=$(printf "%05d" $j)
 # run FoldX
 $FOLDX --command=AnalyseComplex --pdb=frame_${jj}.pdb --analyseComplexChains=A,B --complexWithDNA=false --output-file=out_$jj clean-mode --pdb-dir="../PDBs/" --output-dir="./"
done
# collect "Interaction Energy" of all frames in FoldX.dat
grep -A2 "Interaction Energy" Summary_out_* | grep PDB | awk '{print $6}' > FoldX.dat
# done
cd ../

# 8) ROSETTA scoring
#mkdir ROSETTA; cd ROSETTA
# you might want to parallelize this using a job-array
#$ROSETTA -s ../PDBs/frame_*.pdb -use_jobname 1 -tracer_data_print 1 -out:file:score_only score.dat -add_regular_scores_to_scorefile 1
# collect dG_separated (binding energy) of all frames in ROSETTA.dat
#tail -n +3 score.dat | awk '{print $6}' > ROSETTA.dat
# done
#cd ../
