

/*
 * main.c
 *
 *  Created on: 5.10.2016
 *      Author: heine
 *  Written by Kari Heine.
 *
 * Copyright (C) 2016, Kari Heine, Maria De Iorio, Alex Beskos, Ajay Jasra,
 * David Balding. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. The names of its contributors may not be used to endorse or promote
 *      products derived from this software without specific prior written
 *      permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <unistd.h>
#include <stdbool.h>
#include <sys/stat.h>
#include "ARGdecomposer.h"
#include "commandlineoutput.h"
#include "constants.h"
#include "data.h"
#include "datastructures.h"
#include "debugging.h"
#include "fileaccess.h"
#include "graph2tikz.h"
#include "initialisation.h"
#include "jittering.h"
#include "likelihood.h"
#include "MCMCutils.h"
#include "pathutils.h"
#include "randomness.h"
#include "recombinationTimes.h"
#include "results.h"
#include "shrub.h"
#include "smcPrior.h"
#include "sorting.h"
#include "timeadjustment.h"
#include "treeutils.h"
#include "utils.h"


int *intArray;
char *result_folder;

void deallocateChain(struct MCMCSummary *chain, long iter);
void deallocateDiagnostics(struct MCMCDiagnostics diag);
void deleteBridgePoints(struct BridgePoints bp);
void deleteData(struct Data d);

int main(int argc, const char * argv[]) {
    MPI_Init(NULL, NULL);
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);
//    intArray = malloc(sizeof(int) * 5);
//    intArray[0] = 1;
//
//    printf("Array entry for rank %d is %d\n", world_rank, intArray[0]);


    char *filename;
    struct Data data;
    struct Smc path;
    struct BridgePoints bp;
    struct MCMCSummary *chain;
    struct Parameters parm;
    struct LikelihoodData like_data;
    struct SmcPriorData prior_data;
    struct MRCA mrca;
    long N, iter = 0, full_scan_count = 0;
    short segment_sampler_on = 1;

    clock_t total_program_begin;
    clock_t parallel_begin;
    clock_t parallel_end;
    double total_parallel_time;
    if (world_rank == 0){
        total_program_begin = clock();
        printf("Arbores algorithm for simulating ancestral recombination graphs ");
        printf("conditional of observed DNA polymorphism data.\n");
        printf("Copyright (c) 2016, Kari Heine, Maria De Iorio, Alex Beskos, Ajay Jasra, David Balding\n\n");

        if (argc < 6) {
            printf("Incorrect number of arguments. See instructions for help.\n");
            return 0;
        }
    }

    filename = (char *) argv[1];
    N = (long) atol(argv[2]);
    parm.mu = (double) atof(argv[3]); //1.3e-6;
    parm.rho = (double) atof(argv[4]); //3.5e-7;
    parm.n_eff = 1; //atoi(argv[5]); // 10000;
    init_genrand((unsigned long) atoi(argv[5]));
    result_folder = (char*) argv[6];

    parm.verb = 0;

    if (world_rank ==0 ){
        mkdir(result_folder, S_IRWXU);
        createResultFullPahts(result_folder);
    }

    /* Read a data file */
    data = readData(filename);

    if(data.n_sites >= 30) {
        printf("WARNING: More than 30 recombinations within the data may cause the algorithm be inaccurate or unstable.\n");
        printf("Press any key to continue (or Crtl+C to quit).\n");
        getchar();
    }
    if (parm.verb > 0 && data.M != NULL)
        printData(&data);

    path = initialisation(data, parm);

    /* If initialization introduces recombinations at sites that are not
     * segregating, include non-segregating sites as segregating
     * and augment the data accordingly. */
    data = augmentWithNonSegregatingSites(data, path);

    if (parm.verb > 0)
        printData(&data);

//	bp = createBridgePoints(data, BR_LEN);
    struct BridgePoints bp_one = createPhaseOneBridgePoints(data);
    struct BridgePoints bp_two = createPhaseTwoBridgePoints(data);
    struct BridgePoints bp_three = createPhaseThreeBridgePoints(data);
    struct BridgePoints all_bp[3] = {bp_one, bp_two, bp_three};

    if (world_rank ==0){
        printf("\n%i segments\n\n", (bp_one.length + bp_two.length + bp_three.length));
    }
//	printf("Segregating sites\n");
//	printIntArray(data.segregating_sites, 1, data.n_sites,bp_one 1);

    /* Read initial path from file, if file is provided */
    if (argc >= 8) {
        createInitFilePath((char *) argv[7]);
        deallocatePath(path);
        path = readInitialisationFileRowFormat(data);
        printf("Initialisation read from a file.\n");
    }

    assert(checkTreePathCompletely(path) == 1);
    assert(checkCompatibility(path, data) == 1);

    if (world_rank ==0){
        chain = malloc(sizeof(struct MCMCSummary) * N);
        removeMRCAFile();
        removeChainFile();
        chain[iter].path = path;
        chain[iter].full_scan = 0;
        chain[iter].data.accept_indicator = 1;
        chain[iter].data.alpha = 1;
        chain[iter].data.cardinality_ratio = 1;
        chain[iter].data.current_free_time_density = -1;
        chain[iter].data.current_log_likelihood = -1;
        chain[iter].data.current_log_prior = -1;
        chain[iter].data.current_number_of_free_times = -1;
        chain[iter].data.current_recombination_density = -1;
        chain[iter].data.irreducibility = 0;
        chain[iter].data.jitter_step = 0;
        like_data = likelihood(path, data, parm);
        chain[iter].data.log_likelihood = like_data.log_likelihood;
        prior_data = smcprior(path, parm, data);
        chain[iter].data.log_prior = prior_data.density;
        chain[iter].data.log_posterior = like_data.log_likelihood + prior_data.density;
        deallocateLikelihood(like_data);
        deallocatePriorData(prior_data);
        chain[iter].data.proposed_free_time_density = -1;
        chain[iter].data.proposed_log_likelihood = -1;
        chain[iter].data.proposed_log_prior = -1;
        chain[iter].data.proposed_number_of_free_times = -1;
        chain[iter].data.proposed_number_of_recombinations = countRecombinations(path);
        chain[iter].data.proposed_recombination_density = -1;
    }
    iter++;

//	writeStateToFile(chain,NULL,1);
    if (world_rank ==0){
        writePathToChainFile(path);
    }

    struct Smc old_path;
    struct Smc path_p;
    struct Smc combined_path;
    struct MCMCSummary out;
    struct MCMCDiagnostics dgn;
    struct ShortVector selector;
    struct LikelihoodData like;
    struct SmcPriorData prior;
    int chain_length = iter;
    int old_chain_length;
    struct arraySegmentOutput temp_summary[7];
    struct Smc_array_version old_path_array_form;
    struct arraySegmentOutput to_send;
    int test = 1;

//	initChainTikzFile(chain[iter - 1].path, data);
    MPI_Barrier(MPI_COMM_WORLD);

    if (parm.verb == 0)
        printf("                     ");
    while (true){ //while (true) {

        if (world_rank ==0){
            old_chain_length = chain_length;
            chain_length = jittering(chain, chain_length, N, parm, data);
            iter += (chain_length - old_chain_length);
            if (parm.verb == 1)
                printf("ITERATION %ld\n", iter);

        }
//        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&iter, 1, MPI_INT, 0, MPI_COMM_WORLD);
//        printf("Broadcasted Iter, iter = %d\n", iter);

        if (iter >= N)
            break;

        if (segment_sampler_on == 1) {

            for (int j = 0; j < 3; j++) {
                bp = all_bp[j];

                if (world_rank == 0){
                    convertPathToArray(chain[chain_length-1].path, &old_path_array_form);
                    parallel_begin = clock();
                }


//                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Bcast(&old_path_array_form, sizeof(old_path_array_form), MPI_BYTE, 0, MPI_COMM_WORLD);
                struct Smc old_path_original = convertPathArrayToSmc(old_path_array_form, data);

                if (world_rank < bp.length){
                    struct segment_output output = truncatedSegmentSampler(old_path_original, world_rank, bp, data, parm);
                    convertSegmentOutputToArray(output, &to_send);
                    deallocatePath(output.new_segment);
                }
                deallocatePath(old_path_original);

//                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Gather(&to_send, sizeof(to_send), MPI_BYTE, &temp_summary, sizeof(to_send), MPI_BYTE, 0, MPI_COMM_WORLD);

                if (world_rank == 0){
                    //Combining of segments
                    struct MCMCDiagnostics dgn;
                    struct ShortVector selector;
                    struct MCMCSummary out;

                    struct Smc new_path = combineSegments(&(temp_summary[0]), bp.length, data);

                    for (int i =0; i < bp.length; i++){
                        if (temp_summary[i].accept_indicator == 1)
                            dgn.indicators[i] ='1';
                        else if (temp_summary[i].accept_indicator == 0)
                            dgn.indicators[i] ='0';
                    }
                    dgn.indicators[bp.length] = '\0';
                    parallel_end = clock();
                    total_parallel_time += ((double)(parallel_end - parallel_begin) / CLOCKS_PER_SEC);


                    assert(checkCompatibility(new_path, data) == 1);

                    assert(checkTreePathCompletely(new_path) == 1);
                    assert(checkOperations(new_path) == 1);
                    new_path = removeNoOps(new_path);
                    assert(checkOperations(new_path) == 1);

                    if (new_path.tree_selector != NULL)
                        free(new_path.tree_selector);
                    selector = createTreeSelector(new_path, data);
                    new_path.selector_length = (int) selector.length;
                    new_path.tree_selector = selector.v;

                    dgn = getDiagnostics(data, parm, new_path, dgn);

                    out.data = dgn;
                    out.path = createPathCopy(new_path);

                    //write diagnostics
                    writeDiagnosticsFile(dgn);

                    chain[chain_length] = out;

                    chain_length++;
                    iter += bp.length;

                    deallocatePath(new_path);

                    if (parm.verb == 1)
                        printf("ITERATION %ld\n", iter);
                    else {
                        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%10ld/%-10ld", iter, N);
                        fflush(stdout);
                    }
                    if (j==2) {
                        path = chain[chain_length - 1].path;
                        chain[chain_length - 1].full_scan = 1;
                        mrca = timesToMRCA(path);
                        writeMrcaToFile(mrca);
                        deallocateMRCA(mrca);
                        writePathToChainFile(path);
                        full_scan_count++;
                        map(chain, chain_length, full_scan_count, data);
                    }
                    else
                        chain[chain_length - 1].full_scan = 0;
                }
//                if (iter >= N)
//                    break;
            }
            //End of loop j
//            test += 1;
//            break;
        }
//        if (iter >= N)
//            break;
    }
    if (world_rank ==0){
        deallocateChain(chain, chain_length);
        clock_t total_program_end = clock();
        double total_program_time = (double)(total_program_end - total_program_begin) / CLOCKS_PER_SEC;
        printf("\n Total time for the entire program is %f seconds\n",total_program_time);
        printf("\n Total parallel time is %f seconds\n", total_parallel_time);
    }
    deleteData(data);
    free(bp_one.points);
    free(bp_two.points);
    free(bp_three.points);
    MPI_Finalize();
    return 0;
}

void deallocateChain(struct MCMCSummary *chain, long iter) {
    for (int i = 0; i < iter; i++) {
        deallocatePath(chain[i].path);
        deallocateDiagnostics(chain[i].data);
    }
    free(chain);
}

void deallocateDiagnostics(struct MCMCDiagnostics diag) {
    return;
}

void deleteBridgePoints(struct BridgePoints bp) {
    free(bp.points);
}

void deleteData(struct Data d) {
    free(d.M);
    free(d.name);
    free(d.segregating_sites);
}