ROSETTA="/Users/aug/resetta/3.13/rosetta_bin_mac_2021.16.61629_bundle/main/source/bin/InterfaceAnalyzer.static.macosclangrelease"

# padding for PDB filename
jj=$(printf "%05d" $1)
# run FoldX
if [ -f "../PDBs/frame_${jj}.pdb" ]; then
    $ROSETTA -s ../PDBs/frame_${jj}.pdb -use_jobname 1 -tracer_data_print 1 -out:file:score_only score.dat -add_regular_scores_to_scorefile 1 > rosetta_frame_${jj}.log 2>&1
fi