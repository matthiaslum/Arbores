

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
//#include "/usr/local/include/omp.h"
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
#include <time.h>

int *intArray;
char *result_folder;
//#pragma omp threadprivate(intArray)

void deallocateChain(struct MCMCSummary *chain, long iter);
void deallocateDiagnostics(struct MCMCDiagnostics diag);
void deleteBridgePoints(struct BridgePoints bp);
void deleteData(struct Data d);

int main(int argc, const char * argv[]) {

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

    total_program_begin = clock();

    printf("Arbores algorithm for simulating ancestral recombination graphs ");
    printf("conditional of observed DNA polymorphism data.\n");
    printf("Copyright (c) 2016, Kari Heine, Maria De Iorio, Alex Beskos, Ajay Jasra, David Balding\n\n");

    if (argc < 6) {
        printf("Incorrect number of arguments. See instructions for help.\n");
        return 0;
    }

    filename = (char *) argv[1];
    N = (long) atol(argv[2]);
    parm.mu = (double) atof(argv[3]); //1.3e-6;
    parm.rho = (double) atof(argv[4]); //3.5e-7;
    parm.n_eff = 1; //atoi(argv[5]); // 10000;
    init_genrand((unsigned long) atoi(argv[5]));
    result_folder = (char*) argv[6];

    parm.verb = 0;

    mkdir(result_folder, S_IRWXU);
    createResultFullPahts(result_folder);

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

    printf("\n%i segments\n\n", (bp_one.length + bp_two.length + bp_three.length));
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

    /* main loop */
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
    iter++;

//	writeStateToFile(chain,NULL,1);
    writePathToChainFile(path);

    struct Smc old_path;
    struct MCMCSummary out;
    struct MCMCDiagnostics dgn;
    struct ShortVector selector;
    struct LikelihoodData like;
    struct SmcPriorData prior;
    int chain_length = iter;
    int old_chain_length;
    struct segment_output temp_summary[7];

//	initChainTikzFile(chain[iter - 1].path, data);
    if (parm.verb == 0)
        printf("                     ");
    while (true) {

//        iter = jittering(chain, iter, N, parm, data);
        old_chain_length = chain_length;
        chain_length = jittering(chain, chain_length, N, parm, data);
        iter += (chain_length - old_chain_length);

        if (parm.verb == 1)
            printf("ITERATION %ld\n", iter);
//        path = chain[iter - 1].path;

        if (iter >= N)
            break;
        if (segment_sampler_on == 1) {

            for (int j = 0; j <3; j++) {

                bp = all_bp[j];

//                old_path = chain[iter-1].path;
                old_path = chain[chain_length-1].path;
                parallel_begin = clock();
                for (int i = 0; i < bp.length; i++) {
                    temp_summary[i] = truncatedSegmentSampler(old_path, i, bp, data, parm);
                }

                //Combining of segments
                struct MCMCDiagnostics dgn;
                struct ShortVector selector;
                struct MCMCSummary out;

                struct Smc new_path = combineSegments(temp_summary, bp.length, data);

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

                dgn = getDiagnostics(data, parm, old_path, new_path, dgn);

                out.data = dgn;
                out.path = createPathCopy(new_path);

                //write diagnostics
                writeDiagnosticsFile(dgn);

//                chain[iter] = out;
                chain[chain_length] = out;

                chain_length++;
                iter += bp.length;

//                free(temp_summary);
                deallocatePath(new_path);
//                path = chain[iter - 1].path;
                path = chain[chain_length - 1].path;
//                end = clock();
//                cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
//                printf("Time taken for combining is %f\n", cpu_time_used);

                if (parm.verb == 1)
                    printf("ITERATION %ld\n", iter);
                else {
                    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%10ld/%-10ld", iter, N);
                    fflush(stdout);
                }
                //				appendChainTikzFile(chain[iter - 1].path, data, 0);
                //				writeTikzTexFileForTreePath(NULL, path.tree_path, path.opers, NULL, NULL, NULL,
                //						path.sites, path.rec_times, path.path_len);

                if (j==2) {
                    chain[chain_length -1].full_scan = 1;
//                    chain[iter - 1].full_scan = 1;
                    mrca = timesToMRCA(path);
                    writeMrcaToFile(mrca);
                    deallocateMRCA(mrca);
                    writePathToChainFile(path);
                    full_scan_count++;
//                    map(chain, iter, full_scan_count, data);
                    map(chain, chain_length, full_scan_count, data);
                }
                else {
//                    chain[iter - 1].full_scan = 0;
                    chain[chain_length - 1].full_scan = 0;
                }

                if (iter >= N)
                    break;
                //				getchar();
            } //for loop for j
        } //closing brace for segment sampler
        if (iter >= N)
            break;
    } //closing brace for while (true)

    //	closeChainTikzFile();
//    deallocateChain(chain, iter);
    deallocateChain(chain, chain_length);
    clock_t total_program_end = clock();
    double total_program_time = (double)(total_program_end - total_program_begin) / CLOCKS_PER_SEC;
    printf("\n Total time for the entire program is %f seconds\n",total_program_time);
    printf("\n Total parallel time is %f seconds\n", total_parallel_time);

    //	deleteBridgePoints(bp);
    deleteData(data);
//    free(temp_summary);
    free(bp_one.points);
    free(bp_two.points);
    free(bp_three.points);
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