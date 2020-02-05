

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

    printf("%i segments\n\n", (bp_one.length + bp_two.length + bp_three.length));
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
    struct Smc path_p;
    struct Smc combined_path;
    struct MCMCSummary_modified temp_summary;
    struct MCMCSummary out;
    struct MCMCDiagnostics dgn;
    struct ShortVector selector;
    struct LikelihoodData like;
    struct SmcPriorData prior;

//	initChainTikzFile(chain[iter - 1].path, data);
    if (parm.verb == 0)
        printf("                     ");
    while (true) {

        iter = jittering(chain, iter, N, parm, data);
        if (parm.verb == 1)
            printf("ITERATION %ld\n", iter);
//        path = chain[iter - 1].path;

        if (iter >= N)
            break;

        //variables needed, old_path, combined_path, new_segment, MCMCsummary, diagnostics;
        if (segment_sampler_on == 1) {

            for (int j = 0; j <3; j++) {

                bp = all_bp[j];

//                deallocatePath(old_path);
//                deallocatePath(combined_path);

                old_path = chain[iter-1].path;
                combined_path = createPathCopy(old_path);

                for (int i = 0; i < bp.length; i++) {

//                    chain[iter++] = segmentSampler(path, i, bp, data, parm); //return segment and conditioning that's important
                    temp_summary = truncatedSegmentSampler(old_path, i, bp, data, parm);

                    //new_segment, oldsegment, conditioning
                    //generate complete proposal
                    //generate alpha
                    //generate diagnostics
                    //UPDATE PATH WITH IT

                    path_p = createCompleteProposal(temp_summary.new_segment, combined_path, parm, temp_summary.condition, data);
                    assert(checkCompatibility(path_p, data) == 1);

                    dgn = calculateAlpha(temp_summary.new_segment, path_p, temp_summary.old_segment, combined_path, data, parm, temp_summary.condition,
                                         temp_summary.card_ratio);

                    dgn.u = genrand_real3();
                    dgn.accept_indicator = dgn.u < dgn.alpha ? 1 : 0;
                    dgn.proposed_number_of_recombinations = countRecombinations(path_p);
                    dgn.current_number_of_recombinations = countRecombinations(combined_path);

                    if (dgn.accept_indicator == 1) // accept
                        combined_path = createPathCopy2(path_p);
                    else
                        combined_path = createPathCopy2(combined_path);

                    assert(checkTreePathCompletely(combined_path) == 1);
                    assert(checkOperations(combined_path) == 1);
                    combined_path = removeNoOps(combined_path);
                    assert(checkOperations(combined_path) == 1);

                    if (combined_path.tree_selector != NULL)
                        free(combined_path.tree_selector);
                    selector = createTreeSelector(combined_path, data);
                    combined_path.selector_length = (int) selector.length;
                    combined_path.tree_selector = selector.v;

                    dgn.jitter_step = 0;

                    like = likelihood(combined_path, data, parm);
                    dgn.log_likelihood = like.log_likelihood;
                    prior = smcprior(combined_path, parm, data);
                    dgn.log_prior = prior.density;
                    deallocatePriorData(prior);
                    dgn.log_posterior = dgn.log_likelihood + dgn.log_prior;
                    deallocateLikelihood(like);
                    dgn.irreducibility = temp_summary.irreducibility;


                    free(temp_summary.condition.M);
                    free(temp_summary.condition.sites);
                    deallocatePath(temp_summary.old_segment);
                    deallocatePath(temp_summary.new_segment);
                    deallocatePath(path_p);

                    out.data = dgn;
                    out.path = createPathCopy(combined_path);
                    chain[iter] = out;
                    iter++;

                    path = chain[iter - 1].path;
                    if (parm.verb == 1)
                        printf("ITERATION %ld\n", iter);
                    else {
                        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%10ld/%-10ld", iter, N);
                        fflush(stdout);
                    }
                    //				appendChainTikzFile(chain[iter - 1].path, data, 0);
                    //				writeTikzTexFileForTreePath(NULL, path.tree_path, path.opers, NULL, NULL, NULL,
                    //						path.sites, path.rec_times, path.path_len);

                    if (j==2 && i == bp.length - 1) {
                        chain[iter - 1].full_scan = 1;
                        mrca = timesToMRCA(path);
                        writeMrcaToFile(mrca);
                        deallocateMRCA(mrca);
                        writePathToChainFile(path);
                        full_scan_count++;
                        map(chain, iter, full_scan_count, data);
                    } else {
                        chain[iter - 1].full_scan = 0;
                    }

                    if (iter >= N)
                        break;
                    //				getchar();
                }
            }
            if (iter >= N)
                break;
        }
    }

//	closeChainTikzFile();

    deallocateChain(chain, iter);
//	deleteBridgePoints(bp);
    deleteData(data);
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