# two-alleles-waves-regime-
Two interfering alleles, both are waves 
## waves.cpp:
### Input:
trails: the number of trials the simulation will run, distance: The distance between the sweeping allele and the focal allele,
foc: initial location of the focal allele, wrap: whether periodic bc or not, mig: rate of migration, N: number of inds in a deme,
s2: the fitness of the focal allele, r: recombination rate, infile: name of the file for the data of the profile of the focal allele,
outfile: name of outputfile. 
### Output:

The output is file gives information on every trial. For each trial the first row corresponds to the time the focal allele was inserted
to the system. It gives (for a specific time. The sweeping allele is on the left while the focal allele is on the right): (1) The last deme 
fully occupied by allele 1, (2) The last deme with allele 1, (3) The first deme fully occupied with allele 2, (4) The first deme with 
allele 2, the last deme with allele 2, (5) The last deme fully occupied by allele 2, (6) The number of demes fully occupied by allele 1,
(7) The number of demes fully occupied by allele 2, (8) The total number of allele 2.
The second row for each trial gives 1 or 0 depending whether allele 2 is fixed or extinct. 
### Comments:
Needs gsl to run (uses the random functions). Loads a precalculated density profile of the focal allele from the output of eq.c
