#!/bin/bash
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=40
#SBATCH --time=0-10:00           # time (DD-HH:MM)

module load gnu-parallel
export PATH="/home/<user>/anaconda3/bin:$PATH"
source activate ldsc
scontrol show hostname ${SLURM_JOB_NODELIST} > ./node_list_${SLURM_JOB_ID}


cd /scratch/l/leisun/jianhuig/FINDOR/ldsc/ldsc/

ls /scratch/l/leisun/jianhuig/all_phonotype/heritability_sumstats/ | parallel --sshloginfile ./node_list_${SLURM_JOB_ID} -j 20 --env PATH /home/l/leisun/jianhuig/miniconda3/envs/ldsc/bin/python /scratch/l/leisun/jianhuig/FINDOR/ldsc/ldsc/ldsc.py 
--h2 scratch/l/leisun/jianhuig/all_phonotype/heritability_sumstats/{} 
--ref-ld-chr /scratch/l/leisun/jianhuig/FINDOR/ldsc/ldsc/baseline_v1.2/baseline. 
--w-ld-chr /scratch/l/leisun/jianhuig/FINDOR/ldsc/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. 
--print-coefficients 
--out /scratch/l/leisun/jianhuig/all_phonotype/heritability_sumstats/{} 
--overlap-annot 
--frqfile-chr /scratch/l/leisun/jianhuig/FINDOR/ldsc/ldsc/1000G_Phase3_frq/1000G.EUR.QC.