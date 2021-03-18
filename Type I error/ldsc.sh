#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=0-00:15           # time (DD-HH:MM)

conda activate ldsc

# ldsc regression
python ~/ldsc.py --h2 ${file} --ref-ld-chr ~/baselineLD. --w-ld-chr ~/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --overlap-annot --frqfile-chr ~/1000G_Phase3_frq/1000G.EUR.QC. --out ${file} --print-coefficients

# FINDOR
python FINDOR.py --ref-ld-chr ~/baselineLD. --gwas-data ${file} --regression-results ${file}.results --out ${file}.reweighted
