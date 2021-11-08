#!/bin/bash
#SBATCH --account=nn9244k
#SBATCH --time=01:00:00
#SBATCH --job-name=msprime
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

#set -o errexit


module --quiet purge
module load Anaconda3/2019.03
module load PLINK/1.9b_6.13-x86_64
module load Perl/5.32.0-GCCcore-10.2.0
module load ADMIXTURE/1.3.0
export PS1=\$

source ${EBROOTANACONDA3}/etc/profile.d/conda.sh

conda deactivate &>/dev/null
conda activate /cluster/projects/nn9244k/for_lane/msprime-env


#for chr {1..22}; do for size in $(less sample_size.list); \
#do for dem in $(less dems.list); do \
#sbatch run_msprime.sh $chr $dem $size; done; done; done

chrom=$1
dem_option=$2
sample_pop1=$3

python msprime_demographies.py $chrom $dem_option $sample_pop1

module load PLINK/1.9b_6.13-x86_64

plink --vcf $chrom.$dem_option.$sample_pop1.vcf --double-id --make-bed --out $chrom.$dem_option.$sample_pop1

python Dependencies/bim_fix.py $chrom.$dem_option.$sample_pop1 $chrom.$dem_option.$sample_pop1.fixed.bim $chrom

mv $chrom.$dem_option.$sample_pop1.fixed.bim $chrom.$dem_option.$sample_pop1.bim

plink --bfile $chrom.$dem_option.$sample_pop1 --recode vcf-iid --out $chrom.$dem_option.$sample_pop1
