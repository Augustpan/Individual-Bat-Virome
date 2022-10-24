FOLDX=foldx_20221231

# padding for PDB filename
jj=$(printf "%05d" $1)
# run FoldX
if [ -f "../PDBs/frame_${jj}.pdb" ]; then
    $FOLDX --command=AnalyseComplex --pdb=frame_${jj}.pdb --analyseComplexChains=A,B --complexWithDNA=false --output-file=out_$jj clean-mode --pdb-dir="../PDBs/" --output-dir="./" > foldx_frame_${jj}.log
fi