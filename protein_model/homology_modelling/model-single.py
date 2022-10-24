from modeller import *
from modeller.automodel import *
#from modeller import soap_protein_od

env = Environ()
env.io.hetatm = True
a = AutoModel(env, alnfile='SARS2_like_RBD_with_hACE2-6m0j.ali',
              knowns='6m0j',                                        # Template 
              sequence='SARS2_like_RBD_with_hACE2',                 # Target
              assess_methods=(assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341))
a.starting_model = 1
a.ending_model = 100
a.make()