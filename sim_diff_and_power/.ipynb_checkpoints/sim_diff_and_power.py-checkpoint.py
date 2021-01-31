################################################################
### sim_diff_and_power.py
# - 20201028 Natalie Lowell
#
# PURPOSE: this script runs a simulation model to evalute under which
# migration rates, population sizes, and number of generations of drift
# different Fst values are achieved and how easily they are detected 
# using a chi square test.
#
# How to use this file, minmig_pow_sim_NL.py
# (1) Make sure you're operaing python 3 and have all the necessary modules installed
# (2) Prepare a file of your global allele frequencies, where each line
# includes the number of alleles, followed by the frequency of each allele, delimited by white 
# space. These must add to 1 for each locus. For example, the first rows should look like
#      2    0.935714  0.064286
#      2    0.75  0.25
#      2    0.819643  0.180357
# (3) Prepare a parameters file (see example)
# (4) Call this script using:
#         python minmig_pow_sim_NL.py -a allele_freq_file -p params_file -o output_name

################################################################

# import modules
import argparse
import numpy as np
import matplotlib.pyplot as plt
import simuPOP as sim
import time
import datetime
import scipy.stats as stats

################################################################

# print start time to screen and save for log
startTime = time.time()
print("startTime", datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))

# read in arguments
parser = argparse.ArgumentParser(description="Run simulation model to evalute under which migration rates and population sizes different Fst values are achieved and how much power one has to detect differentiation using a chi square test.")
parser.add_argument("-a", "--afs", help="path to allele frequency file", type=str, required=True)
parser.add_argument("-p", "--params", help="path to parameter file", type=str, required=True)
parser.add_argument("-o", "--output", help="name for output files", type=str, required=True)
args = parser.parse_args()

# read in parameters from parameter file and store in dictionary
params_file = open(args.params, "r")
params = ["reps", "popsizes", "migrates", "genss"]
params_dict = {}
param_index = 0
for line in params_file:
    if line.startswith("#") == False:
        params_dict[params[param_index]] = [float(x) for x in line.split("#")[0].strip().split()]
        param_index += 1
params_file.close()

# write log 
log = open(args.output + "_log.txt", "w")
log.write("startTime " + str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))+"\n")
log.write("\n-----Parameters----\n")
log.write("replicates " + str(params_dict["reps"][0]) + "\n")
log.write("popsizes")
for popsize in params_dict["popsizes"]:
    log.write(" " + str(popsize))
log.write("\n")
log.write("migrationRates ")
for mig_rate in params_dict["migrates"]:
    log.write(" " + str(mig_rate))
log.write("\n")
log.write("generationsDrift")
for gens in params_dict["genss"]:
    log.write(" " + str(gens))
log.write("\n\n")
log.write("-----Simulation-----\n")


# get empirical afs
all_afs_file = open(args.afs, "r") # while working in NB!
afs = {}
locus_index = 0 
for line in all_afs_file:
    afs[locus_index] = []
    linelist = line.strip().split()
    for freq in linelist[1:]:
        afs[locus_index].append(float(freq))
    locus_index += 1
all_afs_file.close()
num_loci = len(afs.keys())

# init results storage
fst = {}
global_chi = {}
dropped_chi = {}
for rep in range(int((params_dict["reps"][0]))):
    fst[rep] = {}
    global_chi[rep] = {}
    dropped_chi[rep] = {}
    for popsize in params_dict["popsizes"]:
        fst[rep][popsize] = {}
        global_chi[rep][popsize] = {}
        dropped_chi[rep][popsize] = {}
        for mig_rate in params_dict["migrates"]:
            fst[rep][popsize][mig_rate] = {}
            global_chi[rep][popsize][mig_rate] = {}
            dropped_chi[rep][popsize][mig_rate] = {}
            for gens in params_dict["genss"]:
                dropped_chi[rep][popsize][mig_rate][gens] = 0
                
# run simulations
for rep in range(int(params_dict["reps"][0])):
    print("rep", rep, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
    log.write("rep " + str(rep) + str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))+"\n")
    for popsize in params_dict["popsizes"]:
        print(" popsize", popsize)
        log.write(" popsize " + str(popsize) + "\n")
        for mig_rate in params_dict["migrates"]:
            print("  migrate", mig_rate)
            log.write("  migrate " + str(mig_rate) + "\n")
            for gens in params_dict["genss"]:
                print("   gens", gens)
                log.write("   gens " + str(gens) + "\n")

                # initialize pop
                pop = sim.Population(size=[popsize]*2, # two subpops
                                        loci=num_loci, # empirical number of loci
                                        subPopNames = ["sp1", "sp2"],
                                    infoFields='migrate_to')
                for locus_index in range(num_loci): # empirical allele frequencies
                    sim.initGenotype(pop, loci=locus_index, freq=afs[locus_index], subPops=[0,1]) 

                # evolve pop across gens
                for generation in range(int(gens)):
                    pop.evolve(
                        initOps=sim.InitSex(), # initialize sex randomly
                        preOps=sim.Migrator(rate=[[0, mig_rate], # by probability (and not proportion)
                                                  [mig_rate, 0]]), 
                        matingScheme=sim.RandomMating(), # random mating
                        postOps=[
                            sim.Stat(popSize=True),
                            sim.Stat(structure=range(num_loci), vars = ['F_st']), # measure Fst
                            sim.Stat(alleleFreq=range(num_loci), vars=['alleleNum','alleleNum_sp'])
                        ],
                        gen = 1
                    )  

                # store final FST in results dictionary
                fst[rep][popsize][mig_rate][gens] = pop.vars()['F_st']

                # run chi square test across loci, store global pvalue
                df = 0
                chi_stats = {}
                chi_pvals = {}
                for locus_index in range(num_loci):
                    allele_counts = []
                    for subpop_index in [0,1]:
                        subpop_row = []
                        this_sp_loc_acs = pop.vars()['subPop'][subpop_index]['alleleNum'][locus_index]
                        alleles = list(range(len(afs[locus_index])))
                        for allele in alleles:
                            if allele in this_sp_loc_acs:
                                subpop_row.append(this_sp_loc_acs[allele])
                            else:
                                subpop_row.append(0)
                        allele_counts.append(subpop_row)
                    
                    # check for expected freq = 0, i.e., no alleles lost from both pops
                    test = True
                    for allele in alleles:
                        count_across_sp_freq_0 = 0
                        for subpop_index in [0,1]:
                            if allele_counts[subpop_index][allele] == 0:
                                count_across_sp_freq_0 += 1
                        if count_across_sp_freq_0 == 2:
                            test=False 
                    if test == False: # could have multiple alleles lost for triallelic loci, so count per locus
                        dropped_chi[rep][popsize][mig_rate][gens] += 1
                            
                    # test as long as no alleles were lost 
                    if test == True:
                        chi_stats[locus_index] = stats.chi2_contingency(observed = allele_counts)[0]
                        chi_pvals[locus_index] = stats.chi2_contingency(observed = allele_counts)[1]
                        df += stats.chi2_contingency(observed = allele_counts)[2]

                global_p = stats.chi2.pdf(x=sum(list(chi_stats.values())), df=df)
                global_chi[rep][popsize][mig_rate][gens] = global_p
                dropped_chi[rep][popsize][mig_rate][gens] = dropped_chi[rep][popsize][mig_rate][gens] / num_loci
                
# write results to file
sp_for_cg = open(args.output + "_globFstChi2_results.txt", "w")
sp_for_cg.write("Rep\tPopsize\tMigRate\tGensDrift\tFst\tChi2Pval\tPropTestsDropped\n")
for rep in range(int(params_dict["reps"][0])):
    for popsize in params_dict["popsizes"]:
        for mig_rate in params_dict["migrates"]:
            for gens in params_dict["genss"]:
                
                # any NaNs? can't be having any NaNs!
                if global_chi[rep][popsize][mig_rate][gens] == np.NaN:
                    print("ruh roh, NaN alert! will cause downstream issues")
                
                sp_for_cg.write(str(rep) + "\t" + str(popsize) + "\t" + str(mig_rate) + "\t" + str(gens) + "\t")
                sp_for_cg.write(str(fst[rep][popsize][mig_rate][gens]) + "\t")
                sp_for_cg.write(str(global_chi[rep][popsize][mig_rate][gens]) + "\t")
                sp_for_cg.write(str(dropped_chi[rep][popsize][mig_rate][gens]) + "\n")
sp_for_cg.close()

log.write("endTime " + str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))+"\n")
log.close()

endTime = time.time()
print("endTime", datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
print("runTime", str(round(endTime-startTime, 2)), "seconds")
