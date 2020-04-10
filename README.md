# Parbores–MW

Parbores–MW is a 'Master–Worker' parallel implementation of Arbores; a Markov chain Monte Carlo (MCMC) algorithm for simulating ancestral
recombination graphs (ARG's) from a Bayesian posterior distribution for given
DNA polymorphism data. The original Arbores can be found here:
https://github.com/heinekmp/Arbores

## Requirements

1. Ensure that you have a **mpicc** compiler installed. Any type *mpicc* is fine. However, the one we use is *MVAPICH* for the *C* programming language.
*(Note that the MVAPICH compiler allows negative tags in the sending of packages)*
2. Ensure that you have **14 cores** on the same node; connected by *Infiniband*. (Best to remotely access high performance computing center)
3. A terminal that accepts linux commands.

## Instructions for compiling and execution.

1. Change directory into the *src*.
2. Compile the program by running the following command: 

**mpicc -std=c99 -I../include main.c ARGdecomposer.c MCMCutils.c backtracking.c commandlineoutput.c constants.c data.c debugging.c exhaustiveSearch.c fileaccess.c freeTimes.c graph2tikz.c initialTimeProposal.c initialisation.c jittering.c likelihood.c mt19937ar.c pathutils.c randomness.c recombinationTimes.c results.c shrub.c smcPrior.c sorting.c timeadjustment.c treeutils.c utils.c -lm -o Parbores_MW**

The above line compiles the executable with the name *Parbores_MW*, located in the src folder. The toy dataset we work with is arb_data4.txt.

3. Run the program with the following arguments. The arguments are the same as the original Arbores. For detailed explanation, look into the manual found in the *doc* file. We add the *mpicc* wrapper and specify the number of processes to execute the program. The following shows a sample line of execution.

**mpirun -np 14 ./Parbores_MW arb_data4.txt 20000 1e-8 0.001 1 output_file**

The above line executed Parbores_MW with 14 cores, 7 for the masters and 7 for the workers. It is run for 20,000 iterations, with mutation and recombination rates of 1e-8 and 0.001 respectively. The random seed is 1, and the results are stored in *output_file*. 

The experiments in the written report are executed for 200,000 iterations instead.
There is a sample job submission script to the NUS HPC.
