PYTHON="/home1/panyuanfei/.conda/envs/md/bin/python"
GMX="/home1/panyuanfei/bin/gmx"
PYMOL="/home1/panyuanfei/pymol/bin/pymol"
FOLDX="/home1/panyuanfei/bin/foldx_20221231"
ROSETTA="/home1/panyuanfei/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin/InterfaceAnalyzer.static.linuxgccrelease"

# do clean up
rm -rf *
mkdir HBONDS
mkdir PLUMED
mkdir PDBs
mkdir FOLDX
mkdir ROSETTA
cp ../../../modified_scripts/plumed_driver_SARS_CoV_2_like.dat PLUMED/plumed_driver.dat

###
### preparing trajectory
###
# 1) first we need to cat together all the production trajectories
$GMX trjcat \
    -f ../2-PRODUCTION/traj_comp.part*.xtc \
    -o traj_all.xtc \
    -cat > /dev/null 2>&1

# 2) now we fix PBCs - print every 10 steps (only protein atoms). This will result in 10001 frames
echo 1 \
    | $GMX trjconv \
        -f traj_all.xtc \
        -o traj_all_PBC.xtc \
        -s ../1-EQUIL/topol.tpr \
        -pbc nojump \
        -skip 10 \
        -n ../0-TOPO/index.ndx \
        -e 1000000 \
    > /dev/null 2>&1

#echo 1 \
#    | $GMX trjconv \
#        -f traj_all.xtc \
#        -o traj_all_PBC.pdb \
#        -s ../1-EQUIL/topol.tpr \
#        -pbc nojump \
#        -skip 10 \
#        -n ../0-TOPO/index.ndx \
#        -e 1000000 \
#    > /dev/null 2>&1

# capture first frame, base RMSD caculation on first frame
echo 1 \
    | $GMX trjconv \
        -f traj_all.xtc \
        -o conf_slt.pdb \
        -s ../1-EQUIL/topol.tpr \
        -pbc nojump \
        -skip 10 \
        -n ../0-TOPO/index.ndx \
        -e 1 \
    > /dev/null 2>&1

$PYMOL -c ../../../modified_scripts/generate_RMSD_templates.py
sed '/^TER   /d' PLUMED/rmsd-INTER-raw.pdb > PLUMED/rmsd-INTER.pdb
rm PLUMED/rmsd-INTER-raw.pdb

# 3) remove temporary trajectory
rm traj_all.xtc

###
### now we do the analysis on the full trajectory
###
# 4) PLUMED post-processing: RMSD calculation and salt bridge analysis
cd PLUMED
plumed driver \
    --plumed plumed_driver.dat \
    --mf_xtc ../traj_all_PBC.xtc > plumed.log
# postprocess COLVAR_SB file to calculate frequency of salt bridge formation
bash ../../../../modified_scripts/calculate_sb.sh
# done
cd ../ 

# 5) Hbonds analysis
cd HBONDS

# preparing the "topol_slt.tpr" file
cp -r ../../0-TOPO/toppar .
cp ../../../../modified_scripts/topol.top .
$GMX grompp \
    -f ../../../../DATA/gromacs-mdps-charmm/0-em-steep.mdp \
    -c ../conf_slt.pdb \
    -p topol.top \
    -o topol_slt.tpr \
    -maxwarn 1 \
    > /dev/null 2>&1
rm mdout.mdp
rm topol.top
rm -r toppar

# execute python script: requires python3 and MDAnalysis
nohup $PYTHON ../../../../modified_scripts/count-HB.py topol_slt.tpr ../traj_all_PBC.xtc  > hb.log 2>&1 &

cd ..

# 6) Create individual PDBs from traj_all_PBC.xtc
# you might want to parallelize this using a job-array
cd PDBs
# execute python script: requires python3 and MDTraj
$PYTHON ../../../../DATA/get_PDB_frames.py ../conf_slt.pdb ../traj_all_PBC.xtc 0 10001
#python ../../../../modified_scripts/get_PDB_frames.py ../conf_slt.pdb ../traj_all_PBC.xtc 0 10001
# done
cd ../

# 7) FoldX scoring
cd FOLDX
seq 0 10000 | parallel -j32 sh ../../../../modified_scripts/do_foldx_remote.sh

# collect "Interaction Energy" of all frames in FoldX.dat
grep -A2 "Interaction Energy" Summary_out_* | grep PDB | awk '{print $6}' > FoldX.dat
# done

# clean up
tar -cf fx_log.tar.gz -I pigz *.log
tar -cf fx_summary_out.tar.gz -I pigz Summary_out_*
tar -cf fx_interface_residues_out.tar.gz -I pigz Interface_Residues_out_*
tar -cf fx_interaction_out.tar.gz -I pigz Interaction_out_*
tar -cf fx_indiv_energies_out.tar.gz -I pigz Indiv_energies_out_*

rm *.log
rm Summary_out_*
rm Interface_Residues_out_*
rm Interaction_out_*
rm Indiv_energies_out_*

cd ../

# 8) ROSETTA scoring
cd ROSETTA
seq 0 10000 | parallel -j32 sh ../../../../modified_scripts/do_rosetta_remote.sh
tail -n +3 score.dat | awk '{print $6}' > ROSETTA.dat
tar -cf - *.log | pigz > rosetta_log.tar.gz
rm *.log
cd ../

# clean up
cd PDBs
tar -cf pdb_frames.tar.gz -I pigz *.pdb
rm *.pdb
cd ..

#Rscript ../../../modified_scripts/visualize.R