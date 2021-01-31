## sim_diff_and_power.py

This script runs a simulation model to evalute under which conditions of migration rates, population sizes, and generations of drift different Fst values are achieved and how easily they are detected using a chi square test. It can be used as a power analysis, as well as to see under what conditions migration and population sizes the empirically observed global Fst can be expected.

The script simulate two populations, connected by some amount of migration. Loci and allele frequencies are parameterized using empirical data, supplied as an input argument. Then, a range of parameter values for migration rates, population sizes, and generations of drift are supplied in a parameter file, as well as a number of replicate times to run the model with each set of input parameters.



#### How to use this file, minmig_pow_sim_NL.py
(1) Make sure you're operaing python 3 and have all the necessary modules installed
(2) Prepare a file of your global allele frequencies, where each line includes the number of alleles, followed by the frequency of each allele, delimited by white space. These must add to 1 for each locus. For example, the first rows should look like

```
      2    0.935714  0.064286
      2    0.75  0.25
      2    0.819643  0.180357
      ```


(3) Prepare a parameters file (see example)
(4) Call this script using:
   ``      python minmig_pow_sim_NL.py -a allele_freq_file -p params_file -o output_name``