/*
 * MCMCutils.c
 *
 *  Created on: 14.11.2016
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

#include <math.h>
#include "backtracking.h"
#include "commandlineoutput.h"
#include "constants.h"
#include "data.h"
#include "datastructures.h"
#include "debugging.h"
#include "exhaustiveSearch.h"
#include "fileaccess.h"
#include "freeTimes.h"
#include "initialTimeProposal.h"
#include "likelihood.h"
#include "MCMCutils.h"
#include "pathutils.h"
#include "randomness.h"
#include "recombinationTimes.h"
#include "smcPrior.h"
#include "sorting.h"
#include "timeadjustment.h"
#include "treeutils.h"
#include "utils.h"

//static short synchronising = 1;
static short verbose = 0;
//static const short set_chooser[21] = {0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,0};
//static const short acc_chooser[21] = {1,0,0,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,1,0};
//short choose_count = 0;
//short a_choose_count = 0;

struct MCMCSummary segmentSampler(struct Smc path_c, int seg, struct BridgePoints bps,
		struct Data data, struct Parameters parm) {

	struct MCMCSummary out;
	struct MCMCDiagnostics dgn;
	struct Conditioning condition;
	struct SamplingSet set;
	struct SamplingSetData sampling_set_data;
	struct Smc path_p, sgmnt_p, sgmnt_c, sampled_path;
	struct ShortVector selector;
	struct LikelihoodData like;
	struct SmcPriorData prior;
	long pick;
	double card_ratio;
	short irreducibility = 0;

	verbose = parm.verb;

	assert(checkOperations(path_c) == 1);

	condition = prepareConditioning(seg, bps, data, path_c);

	if (verbose == 1)
		printSegmentPreamble(condition, data);

	/* We extract a path segment from the current full path and reverse it if we are
	 * in the first segment. This makes the segment comparable with the proposed
	 * segment, which we generate in the reversed direction. */
	sgmnt_c = extractAndReverseSegment(path_c, condition);
	assert(checkOperations(sgmnt_c) == 1);

	sampling_set_data = generateSamplingSet(sgmnt_c, condition, data);
	set = sampling_set_data.set;

	if (verbose == 1)
		printf("Sampling set cardinalities: %ld %ld\n", set.set_cardinalities[0],
				set.set_cardinalities[1]);

	/* in standard setting the sampling set is not empty */
	if (*set.n_paths > 0) {

		pick = drawPathIndex(set);
		card_ratio = calculateCardinalityRatio(sampling_set_data, pick);
		sampled_path = set.paths[pick];

	} else {
		// empty sampling set, the segment is that of the input
//		printf("Empty sampling set\n");
		sampled_path = sgmnt_c;
		card_ratio = (double) 1;
		irreducibility = 1;
	}
	assert(checkTreePathCompletely(sampled_path) == 1);

	if (condition.backward == 1) {
		// turn the right way
		sgmnt_p = reversePathSegment(createPathCopy(sampled_path));
		assert(checkOperations(sgmnt_p) == 1);
	} else {
		sgmnt_p = createPathCopy(sampled_path);
	}
	assert(checkTreePathCompletely(sgmnt_p) == 1);
	path_p = createCompleteProposal(sgmnt_p, path_c, parm, condition, data);
	assert(checkCompatibility(path_p, data) == 1);

	dgn = calculateAlpha(sgmnt_p, path_p, sgmnt_c, path_c, data, parm, condition,
			card_ratio);

	dgn.u = genrand_real3();
	dgn.accept_indicator = dgn.u < dgn.alpha ? 1 : 0;

	if (dgn.accept_indicator == 1) // accept
		out.path = createPathCopy2(path_p);
	else
		out.path = createPathCopy2(path_c);

	assert(checkTreePathCompletely(out.path) == 1);
	assert(checkOperations(out.path) == 1);
	out.path = removeNoOps(out.path);
	assert(checkOperations(out.path) == 1);

	if (out.path.tree_selector != NULL)
		free(out.path.tree_selector);
	selector = createTreeSelector(out.path, data);
	out.path.selector_length = (int) selector.length;
	out.path.tree_selector = selector.v;

	dgn.proposed_number_of_recombinations = countRecombinations(path_p);
	dgn.current_number_of_recombinations = countRecombinations(path_c);
	dgn.jitter_step = 0;

	like = likelihood(out.path, data, parm);
	dgn.log_likelihood = like.log_likelihood;
	prior = smcprior(out.path, parm, data);
	dgn.log_prior = prior.density;
	deallocatePriorData(prior);
	dgn.log_posterior = dgn.log_likelihood + dgn.log_prior;
	deallocateLikelihood(like);
	dgn.irreducibility = irreducibility;

	writeDiagnosticsFile(dgn);

	out.data = dgn;

	free(condition.M);
	free(condition.sites);
	deallocateSamplingSet(set);
	deallocatePath(sgmnt_c);
	deallocatePath(sgmnt_p);
	deallocatePath(path_p);

	return out;
}

long drawPathIndex(struct SamplingSet set) {

	long pick, u;

	if (genrand_res53() < .5 && set.set_cardinalities[0] > 0) {
//	if (set_chooser[choose_count++]==0) {
		// non forced set
//		printf("Choose NON FORCED set\n");
		u = (long) (floor(genrand_res53() * set.set_cardinalities[0]));
		pick = 0;
		for (long i = 0; i < *set.n_paths; i++) {
			if (set.set_selector[i] == 0)
				pick++;
			if (pick - 1 == u) {
				pick = i;
				break;
			}
		}
	} else {
		// forced set
//		printf("Choose FORCED set\n");
		u = (long) (floor(genrand_res53() * set.set_cardinalities[1]));
		pick = 0;
		for (long i = 0; i < *set.n_paths; i++) {
			if (set.set_selector[i] == 1)
				pick++;
			if (pick - 1 == u) {
				pick = i;
				break;
			}
		}
	}

	return pick;
}

double calculateCardinalityRatio(struct SamplingSetData ss_data, long pick) {

	if (ss_data.i_curr >= 0) {
//		printf("Card: %ld %ld\n",
//				ss_data.set.set_cardinalities[ss_data.set.set_selector[pick]],
//				ss_data.set.set_cardinalities[ss_data.set.set_selector[ss_data.i_curr]]);
		return (double) ss_data.set.set_cardinalities[ss_data.set.set_selector[pick]]
				/ (double) ss_data.set.set_cardinalities[ss_data.set.set_selector[ss_data.i_curr]];
	} else {
		return (double) 1;
	}
}

struct Smc extractAndReverseSegment(struct Smc path, struct Conditioning cond) {

	struct Smc out;
	struct LongVector sites;

	assert(path.selector_length>0 && path.tree_selector != NULL);

	sites = getSites(cond);
//	printf("Extracting to sites\n");
//	printLongArray(sites.values, 1, (int) sites.length, 1);

	out = extractToSites(path, sites);

	fillDoubleArray(out.rec_times, -1, 1, out.path_len - 1, 1);
	free(sites.values);

	if (cond.backward == 1)
		out = reversePathSegment(out);

	assignFreeTimeIndicators(out, cond);

	return out;
}

struct Smc extendToFullSequence(struct Smc segment, struct Smc path,
		struct Conditioning cond, struct Data data) {

	struct Smc out;
	long *segment_edges;
	short n_global = 0;
	long i_start, i_end, last, new_path_length, i_tail;
	struct ShortVector selector;

	assert(segment.rec_times!=NULL);

	segment_edges = segmentEndIndices(path, cond);
	i_start = segment_edges[0];
	i_end = segment_edges[1];

	new_path_length = i_start - 1 + segment.path_len + path.path_len - i_end;
	out = createPath(new_path_length);
	out.path_len = new_path_length;
	out.selector_length = 0;

	/* before the segment */
	for (long i = 0; i < i_start; i++) {
		out.tree_path[i] = createCopy(path.tree_path[i]);
		out.opers[i] = createOperationCopy(path.opers[i]);
		out.rec_times[i] = path.rec_times[i];
		out.sites[i] = path.sites[i];
	}

	/* segment
	 * NB: To keep consistent operations, we need to copy the first entry of the segment
	 * (operation 0 --> 1) but the last entry need not be copied. */
	for (long i = 0; i < segment.path_len - 1; i++) {
		out.tree_path[i + i_start] = createCopy(segment.tree_path[i]);
		out.opers[i + i_start] = createOperationCopy(segment.opers[i]);
		out.rec_times[i + i_start] = segment.rec_times[i];
		out.sites[i + i_start] = segment.sites[i];
	}
	out.sites[i_start] = path.sites[i_start];

	/* after the segment */
	for (long i = i_end; i < path.path_len - 1; i++) {
		i_tail = i - i_end + i_start + segment.path_len - 1;
		out.tree_path[i_tail] = createCopy(path.tree_path[i]);
		out.opers[i_tail] = createOperationCopy(path.opers[i]);
		out.rec_times[i_tail] = path.rec_times[i];
		out.sites[i_tail] = path.sites[i];
	}
	/* last tree and site */
	last = path.path_len - 1;
	if (cond.twosided == 0 && cond.backward == 0) {
		out.tree_path[last - i_end + i_start + segment.path_len - 1] = createCopy(
				segment.tree_path[segment.path_len - 1]);
		out.sites[last - i_end + i_start + segment.path_len - 1] =
				segment.sites[segment.path_len - 1];
	} else {
		out.tree_path[last - i_end + i_start + segment.path_len - 1] = createCopy(
				path.tree_path[last]);
		out.sites[last - i_end + i_start + segment.path_len - 1] = path.sites[last];
	}

	out.sites[i_start + segment.path_len - 1] = segment.sites[segment.path_len - 1];

	/* recalculate the global indexing */
	for (int i = 0; i < out.path_len; i++)
		assignGlobalIndices(out.global_index, &n_global, i, out.tree_path);

	selector = createTreeSelector(out, data);
	out.tree_selector = selector.v;
	out.selector_length = (int) selector.length;

	free(segment_edges);
	return out;
}

void deallocateSamplingSet(struct SamplingSet s) {

	for (int i = 0; i < *s.n_paths; i++) {

		for (int j = 0; j < s.paths[i].path_len - 1; j++) {
			deleteTree(s.paths[i].tree_path[j]);
			free(s.paths[i].opers[j]);
		}
		deleteTree(s.paths[i].tree_path[s.paths[i].path_len - 1]);
		free(s.paths[i].tree_path);
		free(s.paths[i].opers);

		if (s.paths[i].global_index != NULL) {
			for (int j = 0; j < s.paths[i].path_len; j++)
				free(s.paths[i].global_index[j]);
			free(s.paths[i].global_index);
		}

		if (s.paths[i].is_free != NULL)
			free(s.paths[i].is_free);
		assert(s.paths[i].rec_times != NULL);
		free(s.paths[i].rec_times);
		free(s.paths[i].sites);
		if (s.paths[i].tree_selector != NULL)
			free(s.paths[i].tree_selector);
	}
	free(s.n_paths);
	free(s.paths);
	free(s.set_cardinalities);
	free(s.set_selector);
}

void printSegmentPreamble(struct Conditioning condition, struct Data data) {
	printf("Segment data:\n");
	printShortArray(condition.M, data.n_seq, condition.length, data.n_seq);
	printf("Sites within the segment:\n");
	for (int i = 0; i < condition.length; i++)
		printf("%6i ", condition.sites[i]);
	printf("\n\n");
}

struct Smc reversePathSegment(struct Smc path) {

	struct Smc reversed;
	short n_global = 0;

	reversed = createPathCopy(path);

// reverse the tree sequence
	for (int i = 0; i < path.path_len; i++) {
		copyTree(&reversed.tree_path[i], path.tree_path[path.path_len - 1 - i]);
		reversed.sites[i] = path.sites[path.path_len - 1 - i];
	}

// update operations
	for (int i = 0; i < path.path_len - 1; i++) {
		free(reversed.opers[i]);
		reversed.opers[i] = calculateOperation(reversed.tree_path[i],
				reversed.tree_path[i + 1]);
	}

// global indices
	if (path.global_index == NULL) {
		reversed.global_index = malloc(sizeof(short*) * path.path_len);
		for (int i = 0; i < path.path_len; i++)
			assignGlobalIndices(reversed.global_index, &n_global, i, reversed.tree_path);
	} else {
		for (int i = 0; i < path.path_len; i++) {
			free(reversed.global_index[i]);
			assignGlobalIndices(reversed.global_index, &n_global, i, reversed.tree_path);
		}
	}

	if (reversed.is_free != NULL)
		free(reversed.is_free);
	reversed.is_free = malloc(sizeof(short) * reversed.path_len);

	if (reversed.rec_times != NULL)
		assert(reversed.rec_times[0] < 0);

	if (reversed.tree_selector != NULL)
		free(reversed.tree_selector);
	reversed.tree_selector = NULL;
	reversed.selector_length = 0;

	deallocatePath(path);
	return reversed;
}

struct Smc createCompleteProposal(struct Smc segment, struct Smc path,
		struct Parameters parm, struct Conditioning condition, struct Data data) {

	struct Smc complete_path;

	/* free time generation */
	if (condition.backward == 1)
		initialTimeProposal(segment, parm, 1);
	else
		freetimes(segment, parm, condition, 1);

	recombinationTimes(segment, 1);
	assert(checkRecombinationTimes(segment) == 1);

	complete_path = extendToFullSequence(segment, path, condition, data);

	return complete_path;
}

struct MCMCDiagnostics calculateAlpha(struct Smc segment_p, struct Smc path_p,
		struct Smc segment_c, struct Smc path_c, struct Data data, struct Parameters parm,
		struct Conditioning condition, double extra) {

	struct MCMCDiagnostics dgn;
	struct LikelihoodData like_p, like_c;
	struct SmcPriorData prior_p, prior_c;
	struct Recombinations rec_c, rec_p;
	double like, prior, f_p, f_c, free_time, rec_time, alpha;

	/* free time generation */
	if (condition.backward == 1) {

		f_p = initialTimeProposal(segment_p, parm, 0);
		f_c = initialTimeProposal(segment_c, parm, 0);
		dgn.proposed_number_of_free_times = -1;
		dgn.current_number_of_free_times = -1;

	} else {

		f_p = freetimes(segment_p, parm, condition, 0);
		f_c = freetimes(segment_c, parm, condition, 0);

		dgn.proposed_number_of_free_times = shortSum(segment_p.is_free, segment_p.path_len);
		dgn.current_number_of_free_times = shortSum(segment_c.is_free, segment_p.path_len);
	}
	free_time = f_c / f_p;

	dgn.proposed_free_time_density = f_p;
	dgn.current_free_time_density = f_c;

	/* recombination times */
	rec_p = recombinationTimes(path_p, 0);
	rec_c = recombinationTimes(path_c, 0);
	rec_time = exp(rec_c.log_p - rec_p.log_p);
	dgn.proposed_recombination_density = rec_p.log_p;
	dgn.current_recombination_density = rec_c.log_p;

	/* likelihood */
	like_p = likelihood(path_p, data, parm);
	like_c = likelihood(path_c, data, parm);
	like = exp(like_p.log_likelihood - like_c.log_likelihood);
	dgn.proposed_log_likelihood = like_p.log_likelihood;
	dgn.current_log_likelihood = like_c.log_likelihood;

	/* prior */
	prior_p = smcprior(path_p, parm, data);
	prior_c = smcprior(path_c, parm, data);
	prior = exp(prior_p.density - prior_c.density);
	dgn.proposed_log_prior = prior_p.density;
	dgn.current_log_prior = prior_c.density;

	/* acceptance probability */
	alpha = like * prior * free_time * rec_time * extra;
	alpha = alpha > 1 ? 1 : alpha;

	dgn.cardinality_ratio = extra;
	if (verbose == 1) {
		printf("Acceptance probability\n----------------------\n");
		printf("%20s: %15.5e %15.5e --> %1.5f\n", "Likelihood", like_p.log_likelihood,
				like_c.log_likelihood, like);
		printf("%20s: %15.5e %15.5e --> %1.5f\n", "Prior", prior_p.density, prior_c.density,
				prior);
		printf("%20s: %15.5e %15.5e --> %1.5f\n", "Free times", f_c, f_p, free_time);
		printf("%20s: %15.5e %15.5e --> %1.5f\n", "Recombination times", rec_c.log_p,
				rec_p.log_p, rec_time);
		printf("%20s: %15.5e %15.5f\n", "Extra", extra, extra);
		printf("%20s: %15.5e %15.5f\n", "Alpha", alpha, alpha);
	}
	deallocateLikelihood(like_p);
	deallocateLikelihood(like_c);
	deallocatePriorData(prior_p);
	deallocatePriorData(prior_c);

	dgn.alpha = alpha;
	return dgn;
}

struct SamplingSetData generateSamplingSet(struct Smc segment_c, struct Conditioning cond,
		struct Data d) {

	struct Smc new_path;
	struct SamplingSetData out;
	struct LinkedSetArray linkedSets;
	struct SamplingSet sampling_set, new_sampling_set, recreations;
	struct ForcingConfig force_config;
	long i_curr = -1, match_count;

	sampling_set.n_sets = 2;
	sampling_set.paths = malloc(sizeof(struct Smc) * MAX_SAMPLING_SET_SIZE);
	sampling_set.set_selector = malloc(sizeof(short) * MAX_SAMPLING_SET_SIZE);
	sampling_set.set_cardinalities = malloc(sizeof(long) * sampling_set.n_sets);
	sampling_set.n_paths = malloc(sizeof(long));
	*sampling_set.n_paths = 0;

	force_config = createForcingConfigurations(cond.length);
	if (verbose == 1)
		printForcingConfigurations(force_config);

	for (int forcing_scheme = 0; forcing_scheme < 2; forcing_scheme++) {

		/* iterate over all forcing configurations */
		for (int i = 0; i < force_config.n_configs; i++) {

			if (verbose == 1) {
				printf("Forcing configuration %3i/%-3i\n", i + 1, force_config.n_configs);
				printf("-----------------------------\n");
			}

			linkedSets = exhaustivePathFinder(cond, force_config.indicators[i], d);

			if (linkedSets.completed == 1) {
				if (cond.twosided == 1) {

					match_count = timeAdjustment(linkedSets, cond, d, sampling_set,
							force_config.indicators[i], verbose);

				} else {

					new_sampling_set = simpleBacktrack(-1, cond.tl, cond.sites[0], linkedSets, d);

					for (long j = 0; j < *new_sampling_set.n_paths; j++) {

						/* add freedom indicator */
						recreations = recreatePath(new_sampling_set.paths[j]);
						new_path = createPathCopy(recreations.paths[0]);
						deallocateSamplingSet(recreations);
						assignFreeTimeIndicators(new_path, cond);

						assert(checkTreePathCompletely(new_path) == 1);
						/* insert to the sampling set */
						sampling_set.paths[*sampling_set.n_paths] = new_path;
						sampling_set.set_selector[*sampling_set.n_paths] = hasExtraRecombinations(
								new_path, force_config.indicators[i], cond);
						sampling_set.n_paths[0]++;

					}

					if (verbose == 1)
						printf("%ld new paths, ", *new_sampling_set.n_paths);

					deallocateSamplingSet(new_sampling_set);
				}
			} else {
				if (verbose == 1)
					printf("Incompleted path construction\n");
			}

			deleteLinkedSetArray(linkedSets);
			if (verbose == 1)
				printf("sampling set size %ld\n\n", *sampling_set.n_paths);
		}

		// if current state is found in the sampling set, no need to enter the emergency
		// forcing iteration of the outer loop
		if ((i_curr = findCurrentSegment(segment_c, sampling_set)) != -1)
			break;

		if (verbose == 1)
			printf("Entering emergency forcing configuration\n");
		deallocateForcingConfiguration(force_config);
		force_config = createEmergencyForcingConfiguration(cond.length);
	}

	if (verbose == 1) {
		printf("\n");
		printf("Path sampling set created with %ld entries\n\n", *sampling_set.n_paths);
	}

	calculateSetSizes(sampling_set);

	deallocateForcingConfiguration(force_config);

	out.set = sampling_set;
	out.i_curr = i_curr;

	return out;
}

long findCurrentSegment(const struct Smc path, struct SamplingSet set) {

	short nb = (path.tree_path[0].n_nodes + 1) / 2 - 1;
	int j, k;

// iterate over paths
	for (long i = 0; i < *set.n_paths; i++) {
		// iterate over trees in the path
		for (j = 0; j < path.path_len; j++) {
			//iterate over the graph matrix entries
			for (k = 0; k < 2 * nb; k++)
				if (path.tree_path[j].C[k] != set.paths[i].tree_path[j].C[k])
					break;
			if (k < 2 * nb) // match
				break;
		}
		if (j == path.path_len)
			return i;
	}

	return -1; // not found
}

void calculateSetSizes(struct SamplingSet set) {
	set.set_cardinalities[0] = 0;
	set.set_cardinalities[1] = 0;
	for (int i = 0; i < *set.n_paths; i++)
		set.set_cardinalities[set.set_selector[i]]++;
}

void deleteLinkedSetArray(struct LinkedSetArray ls) {

	for (int i = 0; i < ls.length; i++) {

		for (int j = 0; j < ls.sets[i].n_trees; j++) {
			deleteTree(ls.sets[i].trees[j]);
			deleteTree(ls.sets[i].trees[j + ls.sets[i].n_trees]);
			free(ls.sets[i].oper[j]);
			free(ls.sets[i].links[j]);
		}
		free(ls.sets[i].trees);
		free(ls.sets[i].links);
		free(ls.sets[i].oper);

		free(ls.sets[i].sites);
		free(ls.sets[i].set_size);
	}

	free(ls.sets);
}

void printForcingConfigurations(struct ForcingConfig fc) {
	printf("Forcing configurations:\n\n");
	for (int i = 0; i < fc.n_configs; i++) {
		for (int j = 0; j < fc.length; j++)
			printf("%hd ", fc.indicators[i][j]);
		printf("\n");
	}
}

struct ForcingConfig createForcingConfigurations(int length) {

	struct ForcingConfig config;

	config.indicators = malloc(sizeof(short *) * (length - 1));
	for (int i = 0; i < length - 1; i++) {
		config.indicators[i] = malloc(sizeof(short) * length);
		fillShortArray(config.indicators[i], 0, 1, length, 1);
		config.indicators[i][i] = 1;
	}

	config.n_configs = length - 1;
	config.length = length;

	return config;
}

struct ForcingConfig createEmergencyForcingConfiguration(int length) {

	struct ForcingConfig config;

	config.indicators = malloc(sizeof(short *));
	config.indicators[0] = malloc(sizeof(short) * length);
	fillShortArray(config.indicators[0], 1, 1, length, 1);
	config.n_configs = 1;
	config.length = length;

	return config;
}

void deallocateForcingConfiguration(struct ForcingConfig config) {
	for (int i = 0; i < config.n_configs; i++)
		free(config.indicators[i]);
	free(config.indicators);
}

struct LinkedSetArray exhaustivePathFinder(struct Conditioning cond, short *force_conf,
		struct Data d) {

	long domain_size, new_domain_size;
	struct Tree *domain;
	short n_n, n_l, n_linked_sets, i_last;
	struct AdjacencySets aSets;
	struct LinkedSetArray pathSet;
	struct LinkedSet *new_linked_sets;

	n_n = cond.tl.n_nodes;
	n_l = (n_n + 1) / 2;

// initial domain
	domain_size = 1;
	domain = malloc(sizeof(struct Tree) * domain_size);
	domain[0] = createCopy(cond.tl);
	if (checkTree(domain[0]) != 1)
		assert(false);

	pathSet.sets = malloc(sizeof(struct LinkedSet) * (cond.length - 1) * 2);
	pathSet.length = 0;

	for (int j = 0; j < cond.length - 1; j++) {

		if (verbose == 1) {
			printf("\tTRANSITION %3i -> %3i, last site %3i,", j, j + 1, cond.length - 1);
			if (force_conf[j] == 1)
				printf(" ALLOW extra recombination\n");
			else
				printf(" NO extra recombinations\n");
			printf("\t---------------------------------------------------------------\n");

			printf("\tSite data: [");
			for (int i = 0; i < d.n_seq; i++)
				printf("%hd ", cond.M[i + j * d.n_seq]);
			printf("] --> [");
			for (int i = 0; i < d.n_seq; i++)
				printf("%hd ", cond.M[i + (j + 1) * d.n_seq]);
			printf("]\n");
			printf("\tDomain size %ld\n", domain_size);
		}

		// Do not allow extra recombinations, i.e. try with no-operation
		if (force_conf[j] == 0)
			aSets = createNoopAdjacencySets(cond, domain, domain_size, j, n_l);
		else
			aSets = createEmptyAdjacencySets(domain_size);

		// the next site could not be reached with no-op
		if (aSets.n_active == 0) {
			if (verbose == 1) {
				if (force_conf[j] == 0)
					printf("\tNo compatible trees with no-op --> recursive search\n\n");
				else
					printf("\tForcing recombination --> recursive search\n\n");
			}
			deleteAdjacencySets(aSets);
			aSets = createAdjacencySetsRecursively(domain, domain_size, cond, j, force_conf);
		} else if (verbose == 1) {
			printf("\t%ld compatible trees with no-operation\n", aSets.n_active);
		}

		if (aSets.n_active == 0) {
			deallocateDomain(domain, domain_size);
			pathSet.completed = 0;
			return pathSet;
		}

		/* convert adjacency sets into linked set structures*/
		assert(checkAdjacencySets(aSets, cond, j) == 1);
		new_linked_sets = linkSets(aSets);
		n_linked_sets = aSets.length - 1;
		assert(checkLastSet(aSets, new_linked_sets, d, j) == 1);
		deleteAdjacencySets(aSets);

		/* merge newly created linked sets into the path set */
		for (int i = 0; i < n_linked_sets; i++)
			pathSet.sets[pathSet.length++] = new_linked_sets[i];
		free(new_linked_sets);

		/* prepare for the next iteration */
		i_last = pathSet.length - 1;
		new_domain_size = pathSet.sets[i_last].n_trees;
		if (verbose == 1)
			printf("\tNew domain size %ld\n", new_domain_size);

		// deallocate the old domain
		for (long i = 0; i < domain_size; i++)
			deleteTree(domain[i]);
		domain = realloc(domain, sizeof(struct Tree) * new_domain_size);

		// new domain consists of the latter trees in the last linked set
		for (long i = 0; i < new_domain_size; i++)
			domain[i] = createCopy(pathSet.sets[i_last].trees[i + new_domain_size]);
		domain_size = new_domain_size;
		if (verbose == 1)
			printf("\n");
	}

	deallocateDomain(domain, domain_size);
	pathSet.completed = 1;
	return pathSet;
}

void deallocateDomain(struct Tree *domain, long domain_size) {
	for (int i = 0; i < domain_size; i++)
		deleteTree(domain[i]);
	free(domain);
}

void printAdjSets(struct AdjacencySets a) {
	for (int i = 0; i < a.n_sets; i++) {
		printf("%i: Set size %ld:", i, a.set_size[i]);
		for (int j = 0; j < a.set_size[i]; j++)
			for (int k = 0; k < a.length; k++)
				printf("%3hd", a.tree_adj_set[i][j + k * a.set_size[i]].n_nodes);
	}
	printf("\n");
}

void deleteAdjacencySets(struct AdjacencySets a) {

	for (int i = 0; i < a.n_sets; i++) {
		if (a.set_size[i] > 0) {
			for (int j = 0; j < a.set_size[i]; j++)
				for (int k = 0; k < a.length; k++)
					deleteTree(a.tree_adj_set[i][j + k * a.set_size[i]]);
			free(a.tree_adj_set[i]);
			free(a.oper_adj_set[i]);
		}
	}
	free(a.set_size);
	free(a.sites);
	free(a.tree_adj_set);
	free(a.oper_adj_set);

}

struct LinkedSet *linkSets(const struct AdjacencySets A) {

	struct AdjSetUnion u;
	struct LinkedSet *out;
	long cnt = 0, n = 0, set_size;
	short len = A.length;

	/* determine the resulting set size */
	for (int i = 0; i < A.n_sets; i++)
		n += A.set_size[i];

	u.n_paths = n;
	u.trees = malloc(sizeof(struct Tree) * n * len);
	u.opers = malloc(sizeof(short) * 2 * n * (len - 1));
	u.links = malloc(sizeof(long) * n * len);
	u.length = len;

// iterate over sets
	for (int i = 0; i < A.n_sets; i++) {

		if (A.set_size[i] > 0) {

			// iterate over the paths in a single set
			for (int j = 0; j < A.set_size[i]; j++) {

				set_size = A.set_size[i];
				for (int k = 0; k < len - 1; k++) { // iterate over the path elements
					// trees
					u.trees[cnt + k * n] = A.tree_adj_set[i][j + k * set_size];
					// operations
					u.opers[cnt + (2 * k + 0) * n] = A.oper_adj_set[i][j + (2 * k + 0) * set_size];
					u.opers[cnt + (2 * k + 1) * n] = A.oper_adj_set[i][j + (2 * k + 1) * set_size];
					// links
					u.links[cnt + k * n] = (k == 0) ? i : cnt;
				}
				// the last tree
				u.trees[cnt + (len - 1) * n] = A.tree_adj_set[i][j + (len - 1) * set_size];
				cnt++;
			}
		}
	}

	out = malloc(sizeof(struct LinkedSet) * (A.length - 1));
	if (A.length == 3) {
		out[0] = plainLinkedSet(u, A.sites);
		out[1] = linkedSetWithEquivalenceClasses(u, A.sites, 1);
	} else {
		out[0] = linkedSetWithEquivalenceClasses(u, A.sites, 0);
	}

	free(u.trees);
	free(u.opers);
	free(u.links);

	return out;
}

struct LinkedSet plainLinkedSet(struct AdjSetUnion u, long *sites) {

	struct LinkedSet out;

	out.trees = malloc(sizeof(struct Tree) * u.n_paths * 2);
	out.oper = malloc(sizeof(short *) * u.n_paths);
	out.links = malloc(sizeof(long *) * u.n_paths);
	out.set_size = malloc(sizeof(long) * u.n_paths);
	out.sites = malloc(sizeof(long) * 2);

	for (int i = 0; i < u.n_paths; i++) {

		out.trees[i] = createCopy(u.trees[i]);
		out.trees[i + u.n_paths] = createCopy(u.trees[i + u.n_paths]);

		out.oper[i] = malloc(sizeof(short) * 2);
		out.oper[i][0] = u.opers[i];
		out.oper[i][1] = u.opers[i + u.n_paths];

		out.links[i] = malloc(sizeof(long));
		out.links[i][0] = u.links[i];

		out.set_size[i] = 1;

	}

	out.length = 2;
	out.n_trees = u.n_paths;
	out.sites[0] = sites[0];
	out.sites[1] = sites[1];

	return out;
}

struct LinkedSet linkedSetWithEquivalenceClasses(struct AdjSetUnion u, long *sites,
		short offset) {

	struct LinkedSet out;
	struct LongVector eq_data;
	struct CanonicalMatrix uniq_rows;
	struct ShortRowSortedMtx sorted_cmtx;
	long end_index, n;
	long *set, *is;
	short *entry;

	n = u.n_paths;

	sorted_cmtx = sortedCanonMtxOfTerminalTrees(u);
	is = sorted_cmtx.i_sort;
	uniq_rows = uniqueRows(sorted_cmtx);
//	printf("Original size: %ld %ld\n",sorted_cmtx.m,sorted_cmtx.n);
//	printf("number of unique canonical forms %ld\n",uniq_rows.m);

	/* find indices of the unique entries in the original canonical matrix */
	out.trees = malloc(sizeof(struct Tree) * uniq_rows.m * 2);
	out.oper = malloc(sizeof(short *) * uniq_rows.m);
	out.links = malloc(sizeof(long *) * uniq_rows.m);
	out.set_size = malloc(sizeof(long) * uniq_rows.m);
	out.sites = malloc(sizeof(long) * 2);

	out.n_trees = uniq_rows.m;
	out.length = 2;
	out.sites[0] = sites[offset];
	out.sites[1] = sites[offset + 1];

	entry = malloc(sizeof(short) * uniq_rows.n);
	end_index = 0;

// iterate over unique terminal tree set
	for (long i = 0; i < uniq_rows.m; i++) {
		// pick an entry from the unique set
		for (int j = 0; j < uniq_rows.n; j++)
			entry[j] = uniq_rows.data[i + j * uniq_rows.m];

		// matching entries in the SORTED but NON-UNIQUE canonical matrix
		eq_data = equivalence(sorted_cmtx, entry, end_index);

		set = eq_data.values;
		end_index += eq_data.length;

		// keep only one tree (pair) in the equivalence class but all ops and links
		out.trees[i] = createCopy(u.trees[is[set[0]] + offset * n]);
		out.trees[i + uniq_rows.m] = createCopy(u.trees[is[set[0]] + (offset + 1) * n]);

		out.set_size[i] = eq_data.length;

		out.oper[i] = malloc(sizeof(short) * 2 * eq_data.length);
		out.links[i] = malloc(sizeof(long) * eq_data.length);

		for (long j = 0; j < eq_data.length; j++) {
			out.oper[i][j] = u.opers[is[set[j]] + 2 * offset * n];
			out.oper[i][j + eq_data.length] = u.opers[is[set[j]] + (2 * offset + 1) * n];
			out.links[i][j] = u.links[is[set[j]] + offset * n];
		}
		free(eq_data.values);
	}
	free(entry);
	free(uniq_rows.data);
	free(sorted_cmtx.data);
	free(sorted_cmtx.i_sort);
	return out;
}

struct ShortRowSortedMtx sortedCanonMtxOfTerminalTrees(struct AdjSetUnion u) {

	struct CanonicalMatrix canon_mtx;
	struct Tree *terminal_trees;

	terminal_trees = getTerminalTrees(u);
	canon_mtx = canoniseTreeSet(terminal_trees, u.n_paths);
	for (int i = 0; i < u.n_paths; i++) {
		free(terminal_trees[i].C);
		free(terminal_trees[i].times);
	}
	free(terminal_trees);

	return sortCanonicalMtx(canon_mtx);

}

struct ShortRowSortedMtx sortCanonicalMtx(struct CanonicalMatrix canon_mtx) {

	struct ShortRowSortedMtx out;
	long *is;
	/* sort the canonical matrix */
	is = malloc(sizeof(long) * canon_mtx.m);
	for (long i = 0; i < canon_mtx.m; i++)
		is[i] = i;
	shortLongSortRows(canon_mtx.data, is, 1, canon_mtx.m, canon_mtx.n, canon_mtx.m); // sorts canon

	out.data = canon_mtx.data;
	out.m = canon_mtx.m;
	out.n = canon_mtx.n;
	out.i_sort = is;
	return out;
}

struct Tree* getTerminalTrees(struct AdjSetUnion u) {

	struct Tree *terminal_trees;
	long n_path = u.n_paths;
	int len = u.length;

	/* collect the terminal trees of all the path segments */
	terminal_trees = malloc(sizeof(struct Tree) * n_path);
	for (long i = 0; i < n_path; i++)
		terminal_trees[i] = createCopy(u.trees[i + (len - 1) * n_path]);

	return terminal_trees;
}

/* Reduce the SORTED CANONICAL MATRIX to only unique entries */
struct CanonicalMatrix uniqueRows(struct ShortRowSortedMtx canon) {

	struct CanonicalMatrix out;
	long count = 0;
	short match;

	out.data = malloc(sizeof(short) * canon.m * canon.n);
	out.m = canon.m;
	out.n = canon.n;

// copy first row
	for (long j = 0; j < out.n; j++)
		out.data[j * out.m] = canon.data[j * canon.m];
	count++;

	for (long r = 1; r < canon.m; r++) {
		match = 1;
		for (int j = 0; j < out.n; j++)
			if (out.data[count - 1 + j * out.m] != canon.data[r + j * out.m]) {
				match = 0;
				break;
			}

		if (match == 0) { // copy row if not equal
			for (long j = 0; j < out.n; j++)
				out.data[count + j * out.m] = canon.data[r + j * canon.m];
			count++;
		}
	}

	/* reshape data to enable reallocation */
	for (long i = 0; i < count * out.n; i++)
		out.data[i] = out.data[i % count + i / count * out.m];

	/* free memory via reallocation */
	out.data = realloc(out.data, sizeof(short) * out.n * count);
	out.m = count;

	return out;
}

struct UniqueRowsWithCounts uniqueRowsWithCounts(struct ShortRowSortedMtx canon) {

	struct UniqueRowsWithCounts out;
	long count = 0;
	short match;

	out.data = malloc(sizeof(short) * canon.m * canon.n);
	out.counts = malloc(sizeof(long) * canon.m);
	out.selector = malloc(sizeof(long) * canon.m);
	out.inds = malloc(sizeof(long) * canon.m);
	for (int i = 0; i < canon.m; i++)
		out.counts[i] = 0;
	out.m = canon.m;
	out.n = canon.n;

// copy first row
	for (long j = 0; j < out.n; j++)
		out.data[j * out.m] = canon.data[j * canon.m];
	out.counts[count]++;
	out.inds[0] = 0;
	out.selector[0] = 0;
	count++;

	for (long r = 1; r < canon.m; r++) {
		match = 1;
		for (int j = 0; j < out.n; j++)
			if (out.data[count - 1 + j * out.m] != canon.data[r + j * out.m]) {
				match = 0;
				break;
			}

		if (match == 0) { // copy row if not equal
			for (long j = 0; j < out.n; j++)
				out.data[count + j * out.m] = canon.data[r + j * canon.m];
			out.counts[count]++;
			out.inds[r] = count;
			out.selector[count] = r;
			count++;
		} else {
			out.counts[count - 1]++;
			out.inds[r] = count - 1;
		}
	}

	/* reshape data to enable reallocation */
	for (long i = 0; i < count * out.n; i++)
		out.data[i] = out.data[i % count + i / count * out.m];

	/* free memory via reallocation */
	out.data = realloc(out.data, sizeof(short) * out.n * count);
	out.counts = realloc(out.counts, sizeof(long) * count);
	out.selector = realloc(out.selector, sizeof(long) * count);
	out.m = count;

	return out;
}

struct LongVector equivalence(struct ShortRowSortedMtx cmtx, short *entry, long start_row) {

	struct LongVector out;
	short match, bigger;
	long match_count = 0, n_rows, n_cols;
	long *tmp;

	n_rows = cmtx.m;
	n_cols = cmtx.n;
	tmp = malloc(sizeof(long) * n_rows);

	for (long i = start_row; i < n_rows; i++) { // iterate ROWS in whole space
		match = 1;
		bigger = 0;

		for (int j = 0; j < n_cols; j++) { // iterate COLUMNS in whole space
			match = cmtx.data[i + j * n_rows] == entry[j] ? 1 : 0;
			bigger = cmtx.data[i + j * n_rows] > entry[j] ? 1 : 0;
			if (match == 0)
				break;
		}
		if (bigger == 1)
			break;
		if (match == 1)
			tmp[match_count++] = i;
	}

	/* Copy found indices from the temporary array to output */
	out.values = malloc(sizeof(long) * match_count);
	out.length = match_count;
	for (long i = 0; i < match_count; i++)
		out.values[i] = tmp[i];
	free(tmp);

	return out;
}

struct AdjacencySets createAdjacencySetsRecursively(struct Tree *domain, long domain_size,
		struct Conditioning cond, int j, short *forced) {

//	struct Tree visu[2];
//	short **visuop;

	struct AdjacencySets aSets;
	struct AdjacencySet aSet, new_aSet;
	short length = 2;
	long new_site = 0;

//	visuop = malloc(sizeof(short*));
	aSets.tree_adj_set = malloc(sizeof(struct Tree *) * domain_size);
	aSets.oper_adj_set = malloc(sizeof(short *) * domain_size);
	aSets.set_size = malloc(sizeof(long) * domain_size);
	aSets.n_active = 0;
	aSets.n_sets = domain_size;

	while (aSets.n_active == 0 && length <= MAX_STEPS + 1) {

		if (verbose == 1)
			printf("\tSearching with %i operations\n\t                       ", length - 1);

		for (long i = 0; i < domain_size; i++) {

			if (verbose == 1) {
				printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%10ld/%-10ld", i + 1,
						domain_size);
				fflush(stdout);
			}

			aSet = exhaustiveSearch(domain[i], cond, (short) j + 1, forced[j], length);

			if (*aSet.n_entries > 0) {

				new_aSet = uniqueOperationSequences(aSet);

				for (int i = 0; i < *aSet.n_entries * aSet.length; i++)
					deleteTree(aSet.trees[i]);

				free(aSet.trees);
				free(aSet.n_entries);
				free(aSet.opers);

				aSets.tree_adj_set[i] = new_aSet.trees;
				aSets.oper_adj_set[i] = extract2DOperations(new_aSet.opers, *new_aSet.n_entries,
						new_aSet.length);

				free(new_aSet.opers);

			} else {
				new_aSet = aSet;
				aSets.tree_adj_set[i] = aSet.trees;
				aSets.oper_adj_set[i] = aSet.opers;
			}

			aSets.set_size[i] = *new_aSet.n_entries;
			aSets.n_active = aSets.set_size[i] > 0 ? aSets.n_active + 1 : aSets.n_active;

			// everything else in aSet is freed when aSets structure is freed
			free(new_aSet.n_entries);
		}

		new_site = length == 3 ? (cond.sites[j] + cond.sites[j + 1]) / 2 : -1;

		aSets.length = length++;
		if (verbose == 1)
			printf("\n");
	}

	if (verbose == 1) {
		printf("\n");
		printf("\tGenerated %ld active sets\n", aSets.n_active);
	}

	aSets.sites = malloc(sizeof(long) * aSets.length);
	if (aSets.length == 2) {
		aSets.sites[0] = cond.sites[j];
		aSets.sites[1] = cond.sites[j + 1];
	} else {
		aSets.sites[0] = cond.sites[j];
		aSets.sites[1] = new_site;
		aSets.sites[2] = cond.sites[j + 1];
	}
	return aSets;
}

short *extract2DOperations(short *op3D, long m, short length) {

	short * out;
	out = malloc(sizeof(short) * m * 2 * (length - 1));
	for (int i = 0; i < m; i++)
		for (int j = 0; j < length - 1; j++) {
			out[i + (2 * j + 0) * m] = op3D[i + (3 * j + 0) * m];
			out[i + (2 * j + 1) * m] = op3D[i + (3 * j + 1) * m];
		}
	return out;
}

struct AdjacencySet uniqueOperationSequences(struct AdjacencySet aSet) {

	struct AdjacencySet out;
	struct CanonicalMatrix cmtx;
	struct UniqueRowsWithCounts uni;
	struct ShortRowSortedMtx sorted;

	cmtx.data = aSet.opers;
	cmtx.m = *aSet.n_entries;
	cmtx.n = 3 * (aSet.length - 1);

	sorted = sortCanonicalMtx(cmtx);
	uni = uniqueRowsWithCounts(sorted);

	out.trees = malloc(sizeof(struct Tree) * uni.m * aSet.length);
	out.opers = malloc(sizeof(short) * uni.m * 3 * (aSet.length - 1));
	out.n_entries = malloc(sizeof(long));
	*out.n_entries = uni.m;
	out.length = aSet.length;

	for (int i = 0; i < uni.m; i++) {
		for (int j = 0; j < out.length; j++)
			out.trees[i + j * uni.m] = createCopy(
					aSet.trees[sorted.i_sort[uni.selector[i]] + j * cmtx.m]);
		for (int j = 0; j < uni.n; j++)
			out.opers[i + j * uni.m] = uni.data[i + j * uni.m];
	}

	free(sorted.i_sort);
	free(uni.counts);
	free(uni.data);
	free(uni.inds);
	free(uni.selector);
	return out;
}

struct AdjacencySets createNoopAdjacencySets(struct Conditioning cond,
		struct Tree *domain, long domain_size, int j, short n_l) {

	struct AdjacencySets adjSets;
	short *cs;

// data column sums
	cs = malloc(sizeof(short) * cond.length);
	shortColumnSum(cond.M, cs, n_l, cond.length, n_l);

	if (j + 1 < cond.length - 1 && (cs[j + 1] <= 1 || cs[j + 1] == n_l)) {
		// next site, which is not the last, has only one one OR all ones.
		// In this case all trees in the domain can be extended to the next
		// site with no-op.
		adjSets = createTrivialAdjSets(domain, domain_size, cond, j);

	} else if (j + 1 < cond.length - 1) {
		// the next site is not the last one, and it is also not a trivial
		// site either. We need to check if any of the trees in the domain
		// are compatible in which case, recombination need not be introduced.
		adjSets = createNonTrivialAdjSets(domain, domain_size, cond, j);
	} else {
		// stepping to the last site.

		// if there is no two sided conditioning, we proceed the same way as for
		// the other sites
		if (cond.tr.C == NULL) {
			// last site, one-sided conditioning
			if (cs[j + 1] <= 1 || cs[j + 1] == n_l) {
				adjSets = createTrivialAdjSets(domain, domain_size, cond, j);
			} else {
				adjSets = createNonTrivialAdjSets(domain, domain_size, cond, j);
			}
		} else {
			// last site with two-sided conditioning

			// Since the conditioning tree is guaranteed to be compatible with the
			// data at the last site, a tree that is structurally equal to the
			// conditioning tree must be compatible
			adjSets = createFinalAdjSets(domain, domain_size, cond, j);
		}
	}
	free(cs);
	return adjSets;
}

struct AdjacencySets createFinalAdjSets(struct Tree *domain, long domain_size,
		struct Conditioning cond, int j) {

	struct AdjacencySets aSets;
	struct CanonicalMatrix canon_R, canon_t;
	short success;

	aSets.tree_adj_set = malloc(sizeof(struct Tree *) * domain_size);
	aSets.oper_adj_set = malloc(sizeof(short *) * domain_size);
	aSets.set_size = malloc(sizeof(long) * domain_size);
	aSets.sites = malloc(sizeof(long) * 2);
	aSets.n_active = 0;
	aSets.n_sets = domain_size;
	aSets.length = 2;
	canon_R = canoniseTreeSet(&cond.tr, 1);
	for (int i = 0; i < domain_size; i++) {
		canon_t = canoniseTreeSet(&domain[i], 1);
		success = 1;
		for (int k = 0; k < canon_t.n * canon_t.m; k++)
			success *= (canon_R.data[k] == canon_t.data[k] ? 1 : 0);
		free(canon_t.data);
		if (success == 1) {
			aSets.tree_adj_set[i] = noopTreeExtension(domain[i]);
			aSets.oper_adj_set[i] = noopOperExtension();
			aSets.set_size[i] = 1;
			aSets.sites[0] = cond.sites[j];
			aSets.sites[1] = cond.sites[j + 1];
			aSets.n_active++;
		} else {
			aSets.set_size[i] = 0;
		}
	}
	free(canon_R.data);
	return aSets;
}

struct AdjacencySets createNonTrivialAdjSets(struct Tree *domain, long domain_size,
		struct Conditioning cond, int j) {

	struct AdjacencySets aSets;

	aSets.tree_adj_set = malloc(sizeof(struct Tree *) * domain_size);
	aSets.oper_adj_set = malloc(sizeof(short *) * domain_size);
	aSets.set_size = malloc(sizeof(long) * domain_size);
	aSets.sites = malloc(sizeof(long) * 2);
	aSets.n_active = 0;
	aSets.n_sets = 0;

	for (int i = 0; i < domain_size; i++)
		if (checkFeasibility(domain[i], cond.M, j + 1) == 1) {
			aSets.tree_adj_set[i] = noopTreeExtension(domain[i]);
			aSets.oper_adj_set[i] = noopOperExtension();
			aSets.set_size[i] = 1;
			aSets.sites[0] = cond.sites[j];
			aSets.sites[1] = cond.sites[j + 1];
			aSets.n_active++;
		} else {
			aSets.set_size[i] = 0;
		}
	aSets.n_sets = domain_size;
	aSets.length = 2;
	return aSets;
}

struct AdjacencySets createEmptyAdjacencySets(long domain_size) {

	struct AdjacencySets aSets;

	aSets.tree_adj_set = malloc(sizeof(struct Tree *) * domain_size);
	aSets.oper_adj_set = malloc(sizeof(short *) * domain_size);
	aSets.set_size = malloc(sizeof(long) * domain_size);
	aSets.sites = malloc(sizeof(long) * 2);
	aSets.n_active = 0;
	aSets.n_sets = 0;

	for (int i = 0; i < domain_size; i++)
		aSets.set_size[i] = 0;
	aSets.n_sets = domain_size;

	return aSets;
}

struct AdjacencySets createTrivialAdjSets(struct Tree *domain, long domain_size,
		struct Conditioning cond, int j) {

	struct AdjacencySets aSets;

	aSets.tree_adj_set = malloc(sizeof(struct Tree *) * domain_size);
	aSets.oper_adj_set = malloc(sizeof(short *) * domain_size);
	aSets.set_size = malloc(sizeof(long) * domain_size);
	aSets.sites = malloc(sizeof(long) * 2);
	aSets.n_active = 0;
	aSets.n_sets = 0;

	for (long i = 0; i < domain_size; i++) {
		aSets.tree_adj_set[i] = noopTreeExtension(domain[i]);
		aSets.oper_adj_set[i] = noopOperExtension();
		aSets.set_size[i] = 1;
		aSets.n_active++;
		aSets.n_sets++;
		aSets.sites[0] = cond.sites[j];
		aSets.sites[1] = cond.sites[j + 1];
	}
	aSets.length = 2;

	return aSets;
}

struct Tree *noopTreeExtension(struct Tree t) {

	struct Tree *set;
	// create the adjacency set for tree domain[i]...
	set = malloc(sizeof(struct Tree) * 2);
	set[0] = createCopy(t);
	set[1] = createCopy(t);
	return set;
}

short * noopOperExtension() {
	short *opers;
	opers = malloc(sizeof(short) * 4);
	opers[0] = -1;
	opers[1] = -1;
	opers[2] = -1;
	opers[3] = -1;
	return opers;
}

struct Conditioning prepareConditioning(int seg, struct BridgePoints bps, struct Data d,
		struct Smc path) {

	struct Conditioning out;
	struct FindResult sites;
	short *M;
	int initial_site, terminal_site;

	sites = findSegmentSites(d, bps.points[seg]);
	initial_site = d.segregating_sites[bps.points[seg].start];
	terminal_site = d.segregating_sites[bps.points[seg].end];

	if (seg == 0) {

		if (verbose == 1)
			printf("One-sided backward conditioning\n");

		/* for the first segment the data has to be flipped left to right */
		M = shortChooseColumns(d.M, sites.indices, sites.length, d.n_seq, d.n_sites);
		out.M = flipArray(M, d.n_seq, sites.length);
		free(M);

		/* for the first segment l.h.s. tree is the terminal tree*/
		out.tl = path.tree_path[path.tree_selector[terminal_site - 1]];

		/* no r.h.s conditioning */
		out.tr.C = NULL;
		out.tr.n_nodes = 0;
		out.tr.times = NULL;
		out.backward = 1;
		out.twosided = 0;

		// reverse the order of sites
		out.sites = malloc(sizeof(int) * sites.length);
		for (int i = 0; i < sites.length; i++)
			out.sites[sites.length - i - 1] = sites.values[i];

	} else {

		if (seg == bps.length - 1) {

			if (verbose == 1)
				printf("One-sided forward conditioning\n");

			out.M = shortChooseColumns(d.M, sites.indices, sites.length, d.n_seq, d.n_sites);
			out.tl = path.tree_path[path.tree_selector[initial_site - 1]];

			out.tr.C = NULL;
			out.tr.n_nodes = 0;
			out.tr.times = NULL;
			out.twosided = 0;

		} else {

			if (verbose == 1)
				printf("Two-sided conditioning\n");

			out.M = shortChooseColumns(d.M, sites.indices, sites.length, d.n_seq, d.n_sites);
			out.tl = path.tree_path[path.tree_selector[initial_site - 1]];
			out.tr = path.tree_path[path.tree_selector[terminal_site - 1]];
			out.twosided = 1;

		}

		out.backward = 0;
		out.sites = malloc(sizeof(int) * sites.length);
		for (int i = 0; i < sites.length; i++)
			out.sites[i] = sites.values[i];

	}

	out.length = sites.length;
	out.n_seq = d.n_seq;

	free(sites.indices);
	free(sites.values);

	assert(checkTree(out.tl) * (out.twosided == 1 ? checkTree(out.tr) : 1) == 1);

	return out;
}

struct LongVector getSites(struct Conditioning cond) {
	struct LongVector out;
	out.length = cond.length;
	out.values = malloc(sizeof(long) * out.length);
	if (cond.backward == 0)
		for (int i = 0; i < cond.length; i++)
			out.values[i] = cond.sites[i];
	else
		for (int i = 0; i < cond.length; i++)
			out.values[i] = cond.sites[cond.length - i - 1];

	return out;
}

struct ShortVector createTreeSelector(struct Smc path, struct Data d) {

	struct ShortVector selector;
	long max_site;
	short i_tree = 0;

	max_site = d.segregating_sites[d.n_sites - 1];
	selector.v = malloc(sizeof(short) * max_site);
	for (int i = 0; i < max_site; i++) {
		if (i_tree < path.path_len - 1 && i >= path.sites[i_tree + 1] - 1)
			i_tree++;
		selector.v[i] = i_tree;
	}
	selector.length = max_site;

	return selector;
}

struct FindResult findSegmentSites(struct Data data, struct BridgePoint bp) {

	struct FindResult out;
	short cnt = 0;
	int *sites, *inds;
	sites = malloc(sizeof(int) * data.n_sites);
	inds = malloc(sizeof(int) * data.n_sites);

	/* count sites in the segment */
	for (int i = 0; i < data.n_sites; i++)
		if (i >= bp.start && i <= bp.end) {
			sites[cnt] = data.segregating_sites[i];
			inds[cnt++] = i;
		}

	/* create output */
	out.values = malloc(sizeof(int) * cnt);
	out.indices = malloc(sizeof(int) * cnt);
	intArrayCopy(out.values, sites, 1, cnt, 1);
	intArrayCopy(out.indices, inds, 1, cnt, 1);
	out.length = cnt;

	free(sites);
	free(inds);
	return out;
}

struct Data augmentWithNonSegregatingSites(struct Data data, struct Smc path) {

	int extra_sites[MAX_RECOMBINATIONS];
	int n_extra = 0, k;
	short found, prev_sum = 1, new_sum = 0;
	short *zeros;

	/* Set difference 'extra_sites = path.sites - data.segregating_sites' */
	for (int i = 0; i < path.path_len - 1; i++) {
		found = 0;
		for (int j = 0; j < data.n_sites; j++)
			if (data.segregating_sites[j] == path.sites[i]) {
				found = 1;
				break;
			}
		if (found == 0)
			extra_sites[n_extra++] = path.sites[i];
	}

	if (n_extra == 0)
		return data;

	zeros = malloc(sizeof(short) * data.n_seq); // create non-segregating site
	fillShortArray(zeros, 0, data.n_seq, 1, data.n_seq);

	for (int i = 0; i < n_extra; i++) {
		for (k = 0; k < data.n_sites; k++)
			if (data.segregating_sites[k] > extra_sites[i])
				break;
		data.M = shortInsertColumn(data.M, k, zeros, data.n_seq, data.n_sites);
		data.segregating_sites = intInsertValue(data.segregating_sites, k, extra_sites[i],
				data.n_sites);
		data.n_sites++;
	}
	free(zeros);

	for (int i = 0; i < data.n_sites; i++) {
		new_sum = 0;
		for (int j = 0; j < data.n_seq; j++) {
			new_sum += data.M[j + i * data.n_seq];
		}
		if (new_sum == 0 && prev_sum == 0) {
			printf(
					"WARNING: Initialisation requires more than two recombinations between consecutive\n");
			printf("segregating sites. Algorithm may be very slow or unstable.\n");
			printf("Press any key to continue (or Ctrl+C to quit)\n");
			getchar();
		}
		prev_sum = new_sum;
	}
	return data;
}
struct BridgePoints createBridgePoints(struct Data d, const int len) {

	struct BridgePoints bps;
	struct BridgePoint *bp;
	int *starts;
	int cnt = 0, step = (int) ceil((double) len / (double) 2), start = step - len + 1;

	starts = malloc(sizeof(int) * d.n_sites);
	starts[cnt++] = start;
	while (starts[cnt - 1] + (step - 1) <= d.n_sites) {
		starts[cnt] = starts[cnt - 1] + (step - 1);
		cnt++;
	}
	cnt = starts[cnt - 1] == d.n_sites ? cnt - 1 : cnt;
	bp = malloc(sizeof(struct BridgePoint) * cnt);
	for (int i = 0; i < cnt; i++) {
		bp[i].start = (starts[i] > 1 ? starts[i] : 1) - 1;
		bp[i].end = (starts[i] + len - 1 < d.n_sites ? starts[i] + len - 1 : d.n_sites) - 1;
	}
	free(starts);
	bps.points = bp;
	bps.length = cnt;

	return bps;
}

struct CanonicalMatrix canoniseTreeSet(struct Tree *trees, long n_trees) {

	struct CanonicalMatrix out;
	int *order;
	short nb;
	int canonical_width;

	nb = (trees[0].n_nodes + 1) / 2 - 1;
	canonical_width = 3 * (int) nb;
	order = malloc(sizeof(int) * nb);

	out.data = malloc(sizeof(short) * n_trees * canonical_width);
	out.m = n_trees;
	out.n = canonical_width;

	for (int i = 0; i < n_trees; i++) {

		/* Sorting */
		for (int j = 0; j < nb; j++)
			order[j] = j;

		sortDoubleArray(order, trees[i].times, nb, 0);

		/* Collect to output */
		for (int j = 0; j < nb; j++) {
			out.data[i + (j + 0 * nb) * n_trees] = trees[i].C[j];
			out.data[i + (j + 1 * nb) * n_trees] = trees[i].C[j + nb];
			out.data[i + (j + 2 * nb) * n_trees] = (short) order[j];
		}
	}

	free(order);

	return out;
}

short *canonisePath(struct Smc path) {

	short nb = (path.tree_path[0].n_nodes + 1) / 2 - 1;
	short *out;

	out = malloc(sizeof(short) * path.path_len * (2 * nb + 2));

	for (int i = 0; i < path.path_len; i++) {
		for (int j = 0; j < nb; j++) {
			out[j + (2 * nb + 2) * i] = path.tree_path[i].C[j];
			out[j + (2 * nb + 2) * i + nb] = path.tree_path[i].C[j + nb];
		}
		if (i < path.path_len - 1) {
			out[(2 * nb + 2) * i + 2 * nb] = path.opers[i][0];
			out[(2 * nb + 2) * i + 2 * nb + 1] = path.opers[i][1];
		} else {
			out[(2 * nb + 2) * i + 2 * nb] = -1;
			out[(2 * nb + 2) * i + 2 * nb + 1] = -1;
		}

	}
	return out;
}

short checkAdjacencySets(struct AdjacencySets a, struct Conditioning cond, int j) {
	for (int i = 0; i < a.n_sets; i++)
		if (a.set_size[i] > 0)
			for (int k = 0; k < a.set_size[i]; k++) {
				if (checkFeasibility(a.tree_adj_set[i][k + (a.length - 1) * a.set_size[i]],
						cond.M, j + 1) != 1) {
					assert(checkTree(a.tree_adj_set[i][k + (a.length - 1) * a.set_size[i]]) == 1);
					printTree(a.tree_adj_set[i][k + (a.length - 1) * a.set_size[i]]);
					return 0;
				}
			}
	return 1;
}

short checkLastSet(struct AdjacencySets a, struct LinkedSet *l, struct Data d, int j) {
	long n_trees = l[a.length - 2].n_trees;
	enum Node_Colors **colors;
	short dataCol;
	struct LinkedSet last_set;
	colors = malloc(sizeof(enum Node_Colors *));
	last_set = l[a.length - 2];
	dataCol = getDataColumn(d, last_set.sites[1]);
	for (int i = 0; i < n_trees; i++) {
		if (checkFeasibility(last_set.trees[i + n_trees], d.M, dataCol) != 1) {
			printf("Failing tree index %i, at %ith conditioning site\n", i, j + 1);
			colors[0] = colorTree(last_set.trees[i + n_trees], d.M, dataCol);
			writeTikzTexFileForTreePath(NULL, &last_set.trees[i + n_trees], NULL, colors, NULL,
			NULL, NULL, NULL, 1);
			free(colors[0]);
			return 0;
		}
	}
	free(colors);
	return 1;
}

short checkCompatibility(struct Smc path, struct Data data) {
	for (int i = 0; i < data.n_sites; i++)
		if (checkFeasibility(
				path.tree_path[path.tree_selector[data.segregating_sites[i] - 1]], data.M, i)
				== 0) {
			for (int j = 0; j < data.n_seq; j++)
				printf("%hd ", data.M[j + i * data.n_seq]);
			printf("\n");
			printTree(path.tree_path[path.tree_selector[data.segregating_sites[i] - 1]]);
			printf("%i %i\n", i, data.segregating_sites[i]);
			return 0;
		}
	return 1;
}

struct SegmentSamplerSynchro synchroniseSegmentSampler(struct Conditioning cond) {

	FILE *file;
	char line[1500];
	const char sep[1] = " ";
	char *token;
	short nl = 9, path_len = 0;
	short **C, **opers;
	double **times;
	int *sites;
	struct SegmentSamplerSynchro out;

	file = fopen(
			"/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/treesegmentdebug.txt", "r");
	C = malloc(sizeof(short*) * 20);
	times = malloc(sizeof(double*) * 20);

	while (fgets(line, sizeof line, file) != NULL) {
		C[path_len] = malloc(sizeof(short) * 2 * (nl - 1));
		times[path_len] = malloc(sizeof(double) * (nl - 1));
		token = strtok(line, sep);
		for (int i = 0; i < nl - 1; i++) {
			C[path_len][i] = (short) atoi(token) - 1;
			token = strtok(NULL, sep);
		}
		for (int i = 0; i < nl - 1; i++) {
			C[path_len][i + nl - 1] = (short) atoi(token) - 1;
			token = strtok(NULL, sep);
		}
		for (int i = 0; i < nl - 1; i++) {
			times[path_len][i] = (double) atof(token);
			token = strtok(NULL, sep);
		}
		path_len++;
	}
	printf("Path length is %i", (int) path_len);

	fclose(file);

	C = realloc(C, sizeof(short *) * path_len);
	times = realloc(times, sizeof(double *) * path_len);

	file = fopen(
			"/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/opersegmentdebug.txt", "r");

	opers = malloc(sizeof(short *) * (path_len - 1));
	for (int i = 0; i < path_len - 1; i++) {
		fgets(line, sizeof line, file);
		opers[i] = malloc(sizeof(short) * 2);
		token = strtok(line, sep);
		opers[i][0] = (short) atoi(token);
		token = strtok(NULL, sep);
		opers[i][1] = (short) atoi(token);
	}
	fclose(file);

	file = fopen("/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/sitesdebug.txt",
			"r");
	sites = malloc(sizeof(int) * path_len);
	fgets(line, sizeof line, file);
	token = strtok(line, sep);
	if (cond.backward == 0) {
		for (int i = 0; i < path_len; i++) {
			sites[i] = (int) atoi(token);
			token = strtok(NULL, sep);
		}
	} else {
		for (int i = 0; i < path_len; i++) {
			sites[path_len - 1 - i] = (int) atoi(token);
			token = strtok(NULL, sep);
		}
	}

	out.C = C;
	out.times = times;
	out.nl = nl;
	out.length = path_len;
	out.opers = opers;
	out.sites = sites;

	return out;
}
