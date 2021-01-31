## sim_diff_and_power.py

This script runs a simulation model to evalute under which conditions of migration rates, population sizes, and generations of drift different Fst values are achieved and how easily they are detected using a chi square test. It can be used as a power analysis, as well as to see under what conditions the empirically observed global Fst can be expected.

The script simulates two populations, connected by some amount of migration. Loci and allele frequencies are parameterized using empirical data, supplied as an input argument (``-a``). Then, a range of parameter values for migration rates, population sizes, and generations of drift are supplied in a parameter file (``-p``), as well as a number of replicate times to run the model with each set of input parameters. 

For each replicate of parameters (migration rate, population size, and generations of drift), it estimates Fst between the two populations, a global (locus) chi-squared statistic, and an associated p-value. To do so, it runs a chi-squared test for each locus and reports this global test and p-value by summing the degrees of freedom and chi-squared statitics across loci. If any locus loses an allele in at least one of the two populations, the test will produce a ``NaN``. For this reason, I added some code to remove these ``NaNs`` and calculate the global test statistic and p-value without them. The model output also reports the proportion of loci that produced ``NaNs`` (and thus had dropped chi-squared tests for some loci) so that the user can interpret results in light of this.

#### How to use this file, minmig_pow_sim_NL.py
1. Make sure you're operaing python 3 and have all the necessary modules installed (i.e.,``argparse``, ``numpy``, ``matplotlib``, ``simuPOP``, ``time``, ``datetime``, ``scipy``)
2. Prepare a file of your global allele frequencies (see example file [here](https://github.com/nclowell/side_projects/blob/main/sim_diff_and_power/all_AFs_for_sp.txt)), where each line includes the number of alleles, followed by the frequency of each allele, delimited by white space. Allele frequencies must add to 1 for each locus. For example, the first rows should look something like this below, if you have biallelic loci. 

```
      2    0.935714  0.064286
      2    0.75  0.25
      2    0.819643  0.180357
```

3. Prepare a parameters file (see [example](https://github.com/nclowell/side_projects/blob/main/sim_diff_and_power/params.txt))
4. Call this script at the command line using:
   ``      python minmig_pow_sim_NL.py -a allele_freq_file -p params_file -o output_name``
   
#### Output

The output of the model is a file called ``*_globFstChi2_results`` (example [here](https://github.com/nclowell/side_projects/blob/main/sim_diff_and_power/test_globFstChi2_results.txt))and it looks like this:

![sim_out](https://github.com/nclowell/side_projects/blob/main/sim_diff_and_power/sim_out.PNG?raw=true)

