# Individual Bat Virome

This repository contains codes and data of the manuscript "Wang et al. (2022) The viromes of individual Chiroptera reveal the co-infection, spillover and emergence risk of potential zoonotic viruses".

## Run

1. To reproduce the statistical analyses in our study, pls run the R codes in `R` folder. The script `virus_table_cleaning.R` should be run first.

2. To conduct homology modelling, run `protein_model/homology_modelling/align2d.py` then `protein_model/homology_modelling/model-single.py`.

3. To conduct molecular dynamics analysis, run `1-do_emin_equil.sh`, `2-do_production.sh` and `3-do_analysis.sh` in the `protein_model/MD/modified_scripts` folder sequentially. Please note that our MD simulation pipeline is modified from https://github.com/maxbonomi/bat-MD (Temmam et al., Discovery of bat coronaviruses close to SARS-CoV-2 and infectious for human cells. Nature 604 (2022) 330.).

## Note

1. Genome sequences of viruses, as well as bat COI gene sequences will be stored in this repository until NCBI accessions become available. See the `sequences` folder.

2. For any technical questions, please write to yfpan21 at m.fudan.edu.cn