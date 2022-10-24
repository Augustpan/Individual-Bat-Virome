# Repository of scripts and GROMACS input (mdp) files
This directory contains the following items:

* `gromacs-mdps-charmm`: GROMACS input files for energy minimization, equilibration, and production;
* `1-do_emin_equil.sh`: protocol for energy minimization and equilibration;
* `2-do_production.sh`: protocol for production simulations;
* `3-do_analysis.sh`: protocol for the analysis of the production simulations;
* `calculate_sb.sh`: bash script to compute the frequency of formation of inter-subunits salt bridges from PLUMED output;
* `count-HB.py`: python script to compute the frequency of formation of inter-subunits hydrogen bonds;
* `get_PDB_frames.py`: convert a GROMACS xtc trajectory to individual PDB files, one for each frame.
