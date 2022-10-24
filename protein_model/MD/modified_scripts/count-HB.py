import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import sys

# get info
def get_info(at):
    # index
    index = at.id + 1
    # name
    name = at.name
    # get chain
    ch = at.segid[-1]
    # get resid
    resid = at.resid
    # correct resid id
    #if(ch=="A"): resid += 19
    #if(ch=="B"): resid -= 264 
    # residue name
    res = at.resname
    # return stuff
    return index,name,res,resid,ch

# load universe
u = MDAnalysis.Universe(sys.argv[1], sys.argv[2])
# number of frames
nframes = float(len(u.trajectory))
# get atoms
atoms = u.atoms
# calculate hydrogen bonds
hbonds = HBA(universe=u)
results = hbonds.run()
# organize results by frequency of formation
hbs = results.count_by_ids()

# open log file
log=open('hbonds.dat','w')
# print stats
for hb in hbs:
    # get info on donor and acceptor
    d_id,d_at,d_resn,d_resid,d_ch = get_info(atoms[hb[0]])
    a_id,a_at,a_resn,a_resid,a_ch = get_info(atoms[hb[2]])
    # print only inter-chains
    if(d_ch!=a_ch):
       log.write("DONOR %6d %4s %4s %4d %1s "      % (d_id, d_at, d_resn, d_resid, d_ch))
       log.write("* ACCEPTOR %6d %4s %4s %4d %1s " % (a_id, a_at, a_resn, a_resid, a_ch))
       # calculate population
       pop = float(hb[3])/nframes
       log.write("* POP %6.3lf\n" % pop)
