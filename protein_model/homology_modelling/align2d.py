from modeller import *

env = Environ()
env.io.hetatm = True
aln = Alignment(env)
mdl = Model(env, file='6m0j', model_segment=('FIRST:A','LAST:F'))
aln.append_model(mdl, align_codes='6m0j', atom_files='6m0j.pdb')
aln.append(file='SARS2_like_RBD_with_hACE2.ali', align_codes='SARS2_like_RBD_with_hACE2')
aln.align2d(max_gap_length=50)
aln.write(file='SARS2_like_RBD_with_hACE2-6m0j.ali', alignment_format='PIR')
aln.write(file='SARS2_like_RBD_with_hACE2-6m0j.pap', alignment_format='PAP')