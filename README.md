# Sequential

*Sequential* is the sequential version of Parbores. The original code can be found here:
https://github.com/heinekmp/Arbores

## Requirements

1. Ensure that you have a **gcc** compiler installed.

## Instructions for compiling and execution.

1. Change directory into the *src*.
2. Compile the program by running the following command: 

**gcc -std=c99 -I../include main.c ARGdecomposer.c MCMCutils.c backtracking.c commandlineoutput.c constants.c data.c debugging.c exhaustiveSearch.c fileaccess.c freeTimes.c graph2tikz.c initialTimeProposal.c initialisation.c jittering.c likelihood.c mt19937ar.c pathutils.c randomness.c recombinationTimes.c results.c shrub.c smcPrior.c sorting.c timeadjustment.c treeutils.c utils.c -lm -o Sequential**

The above line compiles the executable with the name *Sequential*, located in the src folder. The toy dataset we work with is arb_data4.txt.

3. Run the program with the following arguments. The arguments are the same as the original Arbores. For detailed explanation, look into the manual found in the *doc* file. The following shows a sample line of execution.

**./Sequential arb_data4.txt 20000 1e-8 0.001 1 output_file**

The above line executed *Sequential* with 1 core. It is run for 20,000 iterations, with mutation and recombination rates of 1e-8 and 0.001 respectively. The random seed is 1, and the results are stored in *output_file*. 

The experiments in the written report are executed for 200,000 iterations instead.

There is a sample job submission script to the NUS HPC.
