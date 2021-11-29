#!/usr/bin/env python3

import numpy as np
import sys
import msprime

chrom = int(sys.argv[1]) #chromosome
dem_option = sys.argv[2] #which model?
sample_pop1 = int(sys.argv[3]) #specify sample size

def model_constant(chrom, dem_option, sample_pop1):
	print('begin constant population model with args: ' + str(chrom) + ' ' + str(dem_option) + ' ' + str(sample_pop1))
	#first map out the chromosomes
	chroms = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','18',
	'19','20','21','22']
	
	length = ['248956422','242193529','198295559','190214555','181538259','170805979',
	'159345973','145138636','138394717','133797422','135086622','133275309','114364328',
	'107043718','101991189','90338345','83257441','80373285','58617616','64444167',
	'46709983','50818468']
	
	recomb_rate = ['1.14856e-08','1.10543e-08', '1.12796e-08','1.12312e-08','1.12809e-08',
	'1.12229e-08','1.17646e-08','1.14785e-08','1.17807e-08','1.33651e-08','1.17193e-08',
	'1.30502e-08','1.09149e-08','1.11973e-08','1.38358e-08','1.48346e-08','1.58249e-08',
	'1.5076e-08','1.82201e-08','1.71783e-08','1.30452e-08','1.4445e-08']
	#chromosome lengths and recombination rates from stdpopsim catalogue
	
	chrom_tuple = list(zip(chroms, length, recomb_rate))
	chrom_map = pd.DataFrame(chrom_tuple, columns = ['chroms', 'length', 'recomb_rate'])
	
	chrom_map_pos = (1-int(chrom))
	
	#now identify which length and recombination rate we need to use
	chrom_recomb_rate = chrom_map.iloc[chrom_map_pos]['recomb_rate']
	chrom_length = chrom_map.iloc[chrom_map_pos]['length']
	
	print('chrom length and recomb rate set', flush = True)
		#generate the demography
	model = msprime.Demography()
	print('demography set', flush = True)
	#we just have one population 
	model.add_population(name='pop1', initial_size=10000)
	print('population 1 added', flush = True)

	#make sure it writes out in nucleotides
	#model = msprime.InfiniteSites(msprime.NUCLEOTIDES)
	#generate the ancestry
	sim = msprime.sim_ancestry(samples={'pop1' : sample_pop1}, demography=model, random_seed = 13486, recombination_rate=float(chrom_recomb_rate), sequence_length=float(chrom_length), ploidy = 2)
	
	print('model simulated', flush = True)

	steps=[10,20,50,100,200, 400, 600]
	out_file = open('msprime_' + str(chrom) + '_' + str(dem_option) + '_' + str(sample_pop1) + '.log', 'w')
	sys.stdout = out_file
	
	print(msprime.DemographyDebugger(demography=model).population_size_trajectory(steps), flush = True)
	print(model.debug(), flush = True)
	print(model, flush = True)	
	sys.stdout = sys.__stdout__
	out_file.close
	#generate the SNPs
	#using human mutation rate, our ancestry model, random seed
	sim = msprime.sim_mutations(sim, rate = 1.29e-8, random_seed = 145697)

	print('SNPs generated', flush = True)
	#write to VCF
	with open(str(chrom) + '.' + str(dem_option) + '.' + str(sample_pop1) + '.vcf', "w") as vcf_file: 
		sim.write_vcf(vcf_file)
		
	print('VCF created', flush = True)

def model_expansion(chrom, dem_option, sample_pop1):
	print('begin expansion population model with args: ' + str(chrom) + ' ' + str(dem_option) + ' ' + str(sample_pop1))
	#first map out the chromosomes
	chroms = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','18',
	'19','20','21','22']
	
	length = ['248956422','242193529','198295559','190214555','181538259','170805979',
	'159345973','145138636','138394717','133797422','135086622','133275309','114364328',
	'107043718','101991189','90338345','83257441','80373285','58617616','64444167',
	'46709983','50818468']
	
	recomb_rate = ['1.14856e-08','1.10543e-08', '1.12796e-08','1.12312e-08','1.12809e-08',
	'1.12229e-08','1.17646e-08','1.14785e-08','1.17807e-08','1.33651e-08','1.17193e-08',
	'1.30502e-08','1.09149e-08','1.11973e-08','1.38358e-08','1.48346e-08','1.58249e-08',
	'1.5076e-08','1.82201e-08','1.71783e-08','1.30452e-08','1.4445e-08']
	#chromosome lengths and recombination rates from stdpopsim catalogue
	
	chrom_tuple = list(zip(chroms, length, recomb_rate))
	chrom_map = pd.DataFrame(chrom_tuple, columns = ['chroms', 'length', 'recomb_rate'])
	
	chrom_map_pos = (1-int(chrom))
	
	#now identify which length and recombination rate we need to use
	chrom_recomb_rate = chrom_map.iloc[chrom_map_pos]['recomb_rate']
	chrom_length = chrom_map.iloc[chrom_map_pos]['length']
	
	print('chromosome mapped', flush = True)
	
	#generate the demography
	model = msprime.Demography()
	print('demography set', flush = True)
	#we just have one population 
	model.add_population(name='pop1', initial_size=1000000, growth_rate = 0.03)
	print('population 1 added', flush = True)
	#add the population collapse
	#Change population size in the past
	#add instantaneous bottleneck with growth rate of zero afterwards?
	
	model.add_population_parameters_change(time = 150, population = 'pop1', growth_rate = 0)
	print('growth rate added', flush = True)

	#make sure it writes out in nucleotides
	#model = msprime.InfiniteSites(msprime.NUCLEOTIDES)
	#generate the ancestry
	sim = msprime.sim_ancestry(samples={'pop1' : sample_pop1}, demography=model, random_seed = 13486, recombination_rate=float(chrom_recomb_rate), sequence_length=float(chrom_length), ploidy = 2)

	print('model simulated', flush = True)

	#generate the SNPs
	#using human mutation rate, our ancestry model, random seed
	sim = msprime.sim_mutations(sim, rate = 1.29e-8, random_seed = 145697)
	steps=[10,20,50,100,200, 400, 600]
	
	out_file = open('msprime_' + str(chrom) + '_' + str(dem_option) + '_' + str(sample_pop1) + '.log', 'w')
	sys.stdout = out_file
	
	print(msprime.DemographyDebugger(demography=model).population_size_trajectory(steps), flush = True)
	print(model.debug(), flush = True)
	print(model, flush = True)	
	sys.stdout = sys.__stdout__
	out_file.close
		
	print('SNPs generated', flush = True)
	#write to VCF
	with open(str(chrom) + '.' + str(dem_option) + '.' + str(sample_pop1) + '.vcf', "w") as vcf_file: 
		sim.write_vcf(vcf_file)
		
	print('VCF created', flush = True)
	
def model_collapse(chrom, dem_option, sample_pop1):
	print('begin collapse population model with args: ' + str(chrom) + ' ' + str(dem_option) + ' ' + str(sample_pop1))
	#first map out the chromosomes
	chroms = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','18',
	'19','20','21','22']
	
	length = ['248956422','242193529','198295559','190214555','181538259','170805979',
	'159345973','145138636','138394717','133797422','135086622','133275309','114364328',
	'107043718','101991189','90338345','83257441','80373285','58617616','64444167',
	'46709983','50818468']
	
	recomb_rate = ['1.14856e-08','1.10543e-08', '1.12796e-08','1.12312e-08','1.12809e-08',
	'1.12229e-08','1.17646e-08','1.14785e-08','1.17807e-08','1.33651e-08','1.17193e-08',
	'1.30502e-08','1.09149e-08','1.11973e-08','1.38358e-08','1.48346e-08','1.58249e-08',
	'1.5076e-08','1.82201e-08','1.71783e-08','1.30452e-08','1.4445e-08']
	#chromosome lengths and recombination rates from stdpopsim catalogue
	
	chrom_tuple = list(zip(chroms, length, recomb_rate))
	chrom_map = pd.DataFrame(chrom_tuple, columns = ['chroms', 'length', 'recomb_rate'])
	
	chrom_map_pos = (1-int(chrom))
	
	#now identify which length and recombination rate we need to use
	chrom_recomb_rate = chrom_map.iloc[chrom_map_pos]['recomb_rate']
	chrom_length = chrom_map.iloc[chrom_map_pos]['length']
	
	print('chromosome mapped', flush = True)
	
	#generate the demography
	model = msprime.Demography()
	print('demography set', flush = True)
	#we just have one population 
	model.add_population(name='pop1', initial_size=10000, growth_rate = -0.03)
	print('population 1 added', flush = True)
	#add the population expansion
	model.add_population_parameters_change(time=150, population = 'pop1', growth_rate = 0)
	#model.add_population_parameters_change(time = 100, initial_size = 500, population = 'pop1')
	print('growth rate added', flush = True)

	#make sure it writes out in nucleotides
	#model = msprime.InfiniteSites(msprime.NUCLEOTIDES)
	#generate the ancestry
	sim = msprime.sim_ancestry(samples={'pop1' : sample_pop1}, demography=model, random_seed = 13486, recombination_rate=float(chrom_recomb_rate), sequence_length=float(chrom_length), ploidy = 2)

	print('model simulated', flush = True)

	steps=[10,20,50,100,200, 400, 600]
	#generate the SNPs
	#using human mutation rate, our ancestry model, random seed
	sim = msprime.sim_mutations(sim, rate = 1.29e-8, random_seed = 145697)

	print('SNPs generated', flush = True)
	#write to VCF
	
	with open(str(chrom) + '.' + str(dem_option) + '.' + str(sample_pop1) + '.vcf', "w") as vcf_file: 
		sim.write_vcf(vcf_file)
		
	print('VCF created', flush = True)
	
	out_file = open('msprime_' + str(chrom) + '_' + str(dem_option) + '_' + str(sample_pop1) + '.log', 'w')
	sys.stdout = out_file
	
	print(msprime.DemographyDebugger(demography=model).population_size_trajectory(steps), flush = True)
	print(model.debug(), flush = True)
	print(model, flush = True)	
	sys.stdout = sys.__stdout__
	out_file.close
	
def main(chrom, dem_option, sample_pop1):
	if sys.argv[2] == 'constant':
		model_constant(chrom, dem_option, sample_pop1)
	elif sys.argv[2] == 'collapse':
		model_collapse(chrom, dem_option, sample_pop1)
	elif sys.argv[2] == 'expansion':
		model_expansion(chrom, dem_option, sample_pop1)
	else:
		sys.exit('Did you specify constant, collapse, or expansion models?', flush = True)

if __name__ == '__main__':
	main(chrom, dem_option, sample_pop1)
		
