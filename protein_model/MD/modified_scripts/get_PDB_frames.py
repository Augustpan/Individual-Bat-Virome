import mdtraj as md
import sys

# PDB file
PDB_=sys.argv[1]
# TRAJECTORY in GROMACS xtc format
XTC_=sys.argv[2]
# first (included) and last (excluded) frame to be converted
f0=int(sys.argv[3])
f1=int(sys.argv[4])

# load trajectory
t = md.load(XTC_, top=PDB_)

# save frame in PDB format
#for i in range(f0,f1):
#    # add leading zeros
#    l='%05d' % i 
#    # save
#    t[i].save('frame_'+l+'.pdb')

def save_pdb_frame(i):
    l='%05d' % i 
    # save
    t[i].save('frame_'+l+'.pdb')

from multiprocessing import Pool

pool = Pool(2)
pool.map(save_pdb_frame, range(f0,f1))