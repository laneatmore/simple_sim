# simple_sim

Use msprime 1.0 to generate simple demographies

Can be used much like model_admix.py, but only requires three arguments:

arg1 - chrom \
arg2 - demography \
arg3 - population size

Where chromosomes are from 1-22 (chromosome mapping can be edited by hand if desired). The run_msprime.sh script will loop through demographies and population sizes and output VCFs by chromosome. 
