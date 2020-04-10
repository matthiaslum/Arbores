/*
 * timeadjustment.c
 *
 *  Created on: 16.11.2016
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

#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <math.h>
#include "timeadjustment.h"
#include "graph2tikz.h"
#include "backtracking.h"
#include "utils.h"
#include "constants.h"
#include "datastructures.h"
#include "exhaustiveSearch.h"
#include "pathutils.h"
#include "MCMCutils.h"

long timeAdjustment(struct LinkedSetArray array, struct Conditioning cond,
		struct Data data, struct SamplingSet sampling_set, short *extra, short verb) {

	struct Smc new_path;
	struct SamplingSet subset, recreations;
	struct ShortVector topoMatch;
	long n_terminal_trees, match_count = 0, new_paths = 0;
	int i_last;
	short adjust_success;

	i_last = array.length - 1;
	n_terminal_trees = array.sets[array.length - 1].n_trees;

	// iterate over all terminal trees
	for (long i = 0; i < n_terminal_trees; i++) {

		// check consistency with the right conditioning tree
		topoMatch = topologicalMatching(
				&array.sets[i_last].trees[i + array.sets[i_last].n_trees], 1, cond.tr);

		if (topoMatch.v[0] == 1) {

			match_count++;
			// paths leading to the current terminal tree from the linked sets
			subset = simpleBacktrack(i, cond.tl, cond.sites[0], array, data);
			// for each generated path, see if the times can be adjusted
			for (int j = 0; j < *subset.n_paths; j++) {

				/* recreate paths to obtain completely consistent paths with global indexing
				 * and consistent times thgourhgout the path */
				recreations = recreatePath(subset.paths[j]);

				for (int k = 0; k < *recreations.n_paths; k++) {

					adjust_success = adjustTimes(recreations.paths[k], cond.tr);

					// add successfully adjuste path to the sampling set and assign free time indicator
					if (adjust_success == 1) {

						new_path = createPathCopy(recreations.paths[k]);
						assert(checkTreePathCompletely(new_path) == 1);
						assignFreeTimeIndicators(new_path, cond);
						sampling_set.paths[*sampling_set.n_paths] = new_path;
						sampling_set.set_selector[*sampling_set.n_paths] = hasExtraRecombinations(
								new_path, extra, cond);
						sampling_set.n_paths[0]++;
						new_paths++;
						break;
					}
				}
				deallocateSamplingSet(recreations);
			}
			deallocateSamplingSet(subset);
		}
		free(topoMatch.v);
	}

	if (match_count != 1 && verb == 1)
		printf("%ld matches found should be (exactly 1)\n", match_count);

	if (verb == 1) {
		printf("%ld matches found (should be exactly 1)\n", match_count);
		printf("%ld new paths, ", new_paths);
	}
	return match_count;
}

void assignFreeTimeIndicators(struct Smc path, struct Conditioning cond) {

	short length = path.path_len;
	short nb = (path.tree_path[0].n_nodes + 1) / 2 - 1;
	short next_new_node_gi = nb;
	int k;

	assert(path.global_index!=NULL);
	assert(path.is_free!=NULL);

	fillShortArray(path.is_free, 0, 1, length, 1);

	if (cond.twosided == 0) {
		for (int i = 1; i < length; i++)
			for (int j = 0; j < nb; j++)
				// see if next new node is introduced at this site
				if (path.global_index[i][j] == next_new_node_gi) {
					path.is_free[i] = 1;
					next_new_node_gi++;
				}
	} else {
		for (int i = 1; i < length - 1; i++)
			for (int j = 0; j < nb; j++)
				// see if next new node is introduced at this site
				if (path.global_index[i][j] == next_new_node_gi) {
					// see if the new node also is found in the terminal tree
					for (k = 0; k < nb; k++)
						if (next_new_node_gi == path.global_index[length - 1][k])
							break;
					if (k == nb)
						path.is_free[i] = 1;
					next_new_node_gi++;
				}
	}
}

/* Check for each tree in struct Tree *trees, if the graph-matrix and time ordering are
 * equal to t.
 */
struct ShortVector topologicalMatching(struct Tree *trees, long n_trees, struct Tree t) {

	struct ShortVector out;
	short nb = (t.n_nodes + 1) / 2 - 1, ok;
	out.v = malloc(sizeof(long) * n_trees);
	out.length = n_trees;

	for (int i = 0; i < n_trees; i++) {
		ok = 1;
		for (int j = 0; j < nb; j++) {

			assert(t.C[j] < t.C[j + nb]);
			assert(trees[i].C[j] < trees[i].C[j + nb]);

			if ((trees[i].C[j] != t.C[j]) || (trees[i].C[j + nb] != t.C[j + nb])) {
				ok = 0;
				break;
			}
		}
		out.v[i] = ok;
	}

	return out;
}

/* This function builds the tree path again by applying the SPR operations. This is to
 * ensure that trees form a feasible realisation of an SMc process;
 */
struct SamplingSet recreatePath(struct Smc path) {

	struct SamplingSet set;
	struct SimplePath new_path;
	short len = path.path_len;

	set.paths = malloc(sizeof(struct Smc) * MAX_RECREATIONS);
	set.set_selector = malloc(sizeof(short) * MAX_SAMPLING_SET_SIZE);
	set.set_cardinalities = malloc(sizeof(long) * 2);
	set.n_paths = malloc(sizeof(long));
	*set.n_paths = 0;

	new_path = createNewPath(path);

	recreationRecursion(new_path, path, set);

	set.paths = realloc(set.paths, sizeof(struct Smc) * *set.n_paths);

	for (int i = 0; i < len; i++) {
		free(new_path.trees[i].C);
		free(new_path.trees[i].times);
	}
	for (int i = 0; i < len - 1; i++)
		free(new_path.opers[i]);
	free(new_path.trees);
	free(new_path.opers);
	free(new_path.sites);

	return set;
}

void recreationRecursion(struct SimplePath path, struct Smc orig_path,
		struct SamplingSet set) {

	struct Smc new_entry;
	struct SimplePath new_path;
	struct RegraftTimes regraft_times;
	struct Tree tmp_tree1, tmp_tree2;
	short i_last, i_new, new_node, nl, nb, k, n_global;
	short o[2];

	nl = (path.trees[0].n_nodes + 1) / 2;
	nb = nl - 1;

	if (path.length == orig_path.path_len) {
		// Terminate recursion

		/* create new complete Smc */
		new_entry.tree_path = malloc(sizeof(struct Tree) * path.length);
		for (int i = 0; i < path.length; i++)
			new_entry.tree_path[i] = createCopy(path.trees[i]);
		new_entry.opers = malloc(sizeof(short *) * (path.length - 1));
		for (int i = 0; i < path.length - 1; i++) {
			new_entry.opers[i] = malloc(sizeof(short) * 2);
			new_entry.opers[i][0] = path.opers[i][0];
			new_entry.opers[i][1] = path.opers[i][1];
		}
		new_entry.sites = malloc(sizeof(int) * path.length);
		for (int i = 0; i < path.length; i++)
			new_entry.sites[i] = orig_path.sites[i];

		new_entry.is_free = malloc(sizeof(short) * path.length);
		fillShortArray(new_entry.is_free, -1, 1, path.length, 1);
		new_entry.global_index = malloc(sizeof(short*) * path.length);
		new_entry.rec_times = malloc(sizeof(double) * (path.length - 1));
		fillDoubleArray(new_entry.rec_times, -1, 1, path.length - 1, 1);
		new_entry.path_len = orig_path.path_len;
		new_entry.selector_length = -1;
		new_entry.tree_selector = NULL;

		for (int i = 0; i < path.length; i++)
			assignGlobalIndices(new_entry.global_index, &n_global, i, new_entry.tree_path);

		assert(*set.n_paths < MAX_RECREATIONS);
		set.paths[*set.n_paths] = new_entry;
		set.n_paths[0]++;

	} else {
		// deepen the recursion

		tmp_tree1 = createTree(nl);
		tmp_tree2 = createTree(nl);

		i_last = path.length - 1; //
		i_new = i_last + 1;

		new_path.trees = path.trees;
		new_path.opers = path.opers;
		new_path.sites = path.sites;
		new_path.length = path.length + 1;

		// operation (from last tree to the next)
		o[0] = path.opers[i_last][0];
		o[1] = path.opers[i_last][1];

		if (!(o[0] == -1 && o[1] == -1)) { // non no-op

			regraft_times = getRegraftTimeRange(new_path.trees[i_last], o);

			/* apply the SPR */
			copyTree(&tmp_tree1, new_path.trees[i_last]);
			new_node = SPR(tmp_tree1, o);
			timeSorting(tmp_tree1, NULL);

			for (int j = 0; j < regraft_times.n_times; j++) {

				/* try out a regraft time in a temporary tree */
				copyTree(&tmp_tree2, tmp_tree1);
				tmp_tree2.times[new_node - nl] = regraft_times.times[j];
				timeSorting(tmp_tree2, NULL);

				/* check if the resulting tree is a match */
				for (k = 0; k < 2 * nb; k++)
					if (tmp_tree2.C[k] != orig_path.tree_path[i_new].C[k])
						break;

				/* use this regrafting time if it produces a match */
				if (k == 2 * nb) {
					copyTree(&new_path.trees[i_new], tmp_tree2);
					recreationRecursion(new_path, orig_path, set);
				}
			}
			free(regraft_times.times);
		} else {
			copyTree(&new_path.trees[i_new], new_path.trees[i_last]);
			recreationRecursion(new_path, orig_path, set);
		}
		free(tmp_tree1.C);
		free(tmp_tree1.times);
		free(tmp_tree2.C);
		free(tmp_tree2.times);
	}
}

void assignGlobalIndices(short **global_idx, short *n_global, short step,
		struct Tree *trees) {

	short nl = (trees[0].n_nodes + 1) / 2;
	short nb = nl - 1, matching_node;

	global_idx[step] = malloc(sizeof(short) * nb);

	if (step == 0) {
		for (int i = 0; i < nb; i++)
			global_idx[step][i] = i;
		*n_global = nb;
	} else {
		for (int i = 0; i < nb; i++) {
			matching_node = findMatchingNode(trees[step - 1], trees[step].times[i]);
			if (matching_node < 0) {
				global_idx[step][i] = *n_global;
				n_global[0]++;
			} else {
				global_idx[step][i] = global_idx[step - 1][matching_node];
			}
		}
	}
}

short adjustTimes(struct Smc path, struct Tree tr) {

	struct LimitSummary summary;
	short adjustment_success = 0;

	if (checkFixedTimeCondition(path, tr) == 0)
		return adjustment_success;

	summary = firstPass(path, tr);
	adjustment_success = secondPass(summary, path, tr);
	deallocateSummary(summary);

	return adjustment_success;
}

void deallocateSummary(struct LimitSummary summary) {
	free(summary.conditioning_index);
	free(summary.adjusted_times);
	free(summary.freedom_indicator);
	free(summary.is_new);
	free(summary.new_node);
	for (int i = 0; i < summary.length; i++) {
		if (summary.limits[i].glb_inds != NULL)
			free(summary.limits[i].glb_inds);
		if (summary.limits[i].inds != NULL)
			free(summary.limits[i].inds);
		if (summary.limits[i].is_free != NULL)
			free(summary.limits[i].is_free);
		if (summary.limits[i].times != NULL)
			free(summary.limits[i].times);
	}
	free(summary.limits);
	free(summary.glb_times.step);
	free(summary.glb_times.t);
}

short checkFixedTimeCondition(struct Smc path, struct Tree tr) {

	short nb = (tr.n_nodes + 1) / 2 - 1;
	short length = path.path_len;

	/* if a node in the first tree in the path has the same global index as
	 * a node in the right condition, then their times must be equal */
	for (int i = 0; i < nb; i++)
		for (int j = 0; j < nb; j++)
			if (path.global_index[0][i] == path.global_index[length - 1][j]) {
//				printf("%hd %hd %15.10f %15.10f\n", path.global_index[0][i],
//						path.global_index[length - 1][j], tr.times[j], path.tree_path[0].times[i]);
				if (fabs(tr.times[j] - path.tree_path[0].times[i])
						> (path.tree_path[0].times[i] + tr.times[j]) * TOL)
					return 0;
			}

	for (int i = 0; i < nb; i++)
		for (int j = 0; j < nb; j++)
			if (fabs(path.tree_path[0].times[i] - tr.times[j])
					< (path.tree_path[0].times[i] + tr.times[j]) * TOL) {
//				printf("%hd %hd %15.10f %15.10f\n", path.global_index[0][i],
//						path.global_index[length - 1][j], tr.times[j], path.tree_path[0].times[i]);
				if (path.global_index[0][i] != path.global_index[length - 1][j])
					return 0;
			}
	return 1;
}

struct LimitSummary firstPass(const struct Smc path, struct Tree tr) {

	struct Tree t_curr, t_next;
	short nn = path.tree_path[0].n_nodes;
	short nl = (nn + 1) / 2;
	short nb = nl - 1;
	short n_glb_times, new_node_gi, length = path.path_len;
	short o[2];
	double new_node_time;
	struct LimitSummary smry;

	/* set up global time book keeping */
	n_glb_times = length + nb - 1;
	smry.glb_times.t = malloc(sizeof(double) * n_glb_times);
	smry.glb_times.step = malloc(sizeof(short) * n_glb_times);
	for (int i = 0; i < nb; i++) {
		smry.glb_times.t[i] = path.tree_path[0].times[i];
		smry.glb_times.step[i] = 0;
	}
	smry.glb_times.n_times = nb;

	/* set up other arrays in the summary */
	smry.is_new = malloc(sizeof(short) * length);
	smry.is_new[0] = 0; // first tree does not contain a 'new node'
	smry.conditioning_index = malloc(sizeof(short) * length);
	smry.freedom_indicator = malloc(sizeof(short) * length);
	smry.adjusted_times = malloc(sizeof(double) * length);
	smry.new_node = malloc(sizeof(short) * length);
	smry.conditioning_index[0] = -1;
	smry.freedom_indicator[0] = -1;
	smry.adjusted_times[0] = -1;
	smry.new_node[0] = -1;
	smry.length = length;

	/* the limiting data itself */
	smry.limits = malloc(sizeof(struct LimitData) * length);
	smry.limits[0].inds = NULL;
	smry.limits[0].is_free = NULL;
	smry.limits[0].times = NULL;
	smry.limits[0].glb_inds = NULL;
	smry.limits[0].step = -1;

	/* iterate over the path and extract limiter data for each new node in the path */
	for (int i = 0; i < length - 1; i++) {

		o[0] = path.opers[i][0];
		o[1] = path.opers[i][1];
		if (o[0] == -1 && o[1] == -1) {
			// identity operation does not introduce a 'new node'
			smry.is_new[i + 1] = 0;
			smry.conditioning_index[i + 1] = -1;
			smry.freedom_indicator[i + 1] = -1;
			smry.limits[i + 1].inds = NULL;
			smry.limits[i + 1].is_free = NULL;
			smry.limits[i + 1].times = NULL;
			smry.limits[i + 1].glb_inds = NULL;
			smry.limits[i + 1].step = -1;
			smry.adjusted_times[i + 1] = -1;
			smry.new_node[i + 1] = -1;
			continue;
		}
		smry.is_new[i + 1] = 1;

		t_curr = path.tree_path[i];
		t_next = path.tree_path[i + 1];

		smry.new_node[i + 1] = determineNewNode(t_curr, t_next) + nl;
		new_node_gi = path.global_index[i + 1][smry.new_node[i + 1] - nl]; // global index
		new_node_time = path.tree_path[i + 1].times[smry.new_node[i + 1] - nl];

		smry.glb_times.t[smry.glb_times.n_times] = new_node_time;
		smry.glb_times.step[smry.glb_times.n_times++] = i + 1;

		smry.conditioning_index[i + 1] = findConditioningIndex(path.global_index[length - 1],
				new_node_gi, nb);
		smry.freedom_indicator[i + 1] = smry.conditioning_index[i + 1] < 0 ? 1 : 0;

		smry.limits[i + 1] = getLimiters(t_curr, new_node_time, i, path);

		smry.adjusted_times[i + 1] =
				smry.limits[i + 1].times[1] == DBL_MAX ?
						DBL_MAX : (smry.limits[i + 1].times[0] + smry.limits[i + 1].times[1]) / 2;
	}
	return smry;
}

short secondPass(struct LimitSummary smry, struct Smc path, struct Tree tr) {

	struct Smc path_in;
	struct Tree t_prev, t_curr;
	short nl, len = smry.length, new_node_idx, new_node_gi, i_curr, not_fixed, push_dir,
			need_to_push, ok = 1, it, consistency_check_done = 0;
	short *is_limiter, *is_limiter_tmp, *limit_inds;
	struct DoubleVector u_bounders, l_bounders;
	double trgt_time, incr, new_time, lim_time, upp_bound, low_bound;
	double *new_times;
	double limit_times[2];
	char infile[] = "/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/pathin.tex";
//	char intrmfile[] = "/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/interm.tex";

	nl = (tr.n_nodes + 1) / 2;

	assert(checkTreePathStructure(path) == 1);

	/* copy of the inut is required to make sure that the path structure has
	 * not changed while adjusting the times */
	path_in = createPathCopy(path);

	is_limiter = malloc(sizeof(short) * len);
	fillShortArray(is_limiter, 0, 1, len, 1);

	/* iterate over the path from RIGHT to LEFT */
	for (int i = len - 1; i >= 0; i--) {

		// to proceed, there must be a 'new node' and it must also be in the terminal tree
		if (!(smry.is_new[i] == 1 && smry.conditioning_index[i] >= 0))
			continue;

		assert(smry.limits[i].glb_inds!=NULL);

		trgt_time = tr.times[smry.conditioning_index[i] - nl]; // target time

		if (trgt_time < smry.limits[i].times[1] && trgt_time > smry.limits[i].times[0]) {
			// target time in range --> no need to push limits

			new_node_idx = smry.new_node[i];
			new_node_gi = path.global_index[i][new_node_idx - nl];

			forwardTimeAdjustment(path, i, trgt_time, new_node_idx, new_node_gi);

			is_limiter[i] = 1;
			smry.adjusted_times[i] = trgt_time;

			updateLimiters(smry, i, trgt_time, nl, new_node_gi);

		} else {

			// a limit has to be pushed, determine the direction
			if (trgt_time <= smry.limits[i].times[0] + smry.limits[i].times[0] * TOL) {
				incr = -EPS;
				push_dir = 0; // lower limit
//				printf("%8.2f not in [%8.2f,%8.2f]--> pushing lower limit\n", trgt_time,
//						smry.limits[i].times[0] + TOL, smry.limits[i].times[1]);
			} else {
				incr = EPS;
				push_dir = 1; // upper limit
//				printf("%8.2f not in [%8.2f,%8.2f]--> pushing upper limit\n", trgt_time,
//						smry.limits[i].times[0], smry.limits[i].times[1]);
			}

			/* Initiate a while loop, which moves all the chosen boundaries
			 * throughout the path */
			i_curr = i;
			new_time = trgt_time;

			// is the limit, that needs to be pushed, free?
			not_fixed = smry.limits[i_curr].is_free[push_dir];
			lim_time = smry.limits[i_curr].times[push_dir];

			is_limiter_tmp = malloc(sizeof(short) * smry.length);
			new_times = malloc(sizeof(double) * smry.length);
			fillShortArray(is_limiter_tmp, 0, 1, smry.length, 1);
			fillDoubleArray(new_times, 0, 1, smry.length, 1);
			new_times[i_curr] = trgt_time;
			is_limiter_tmp[i_curr] = 1;

			it = 0; // debug

			while (not_fixed == 1) {

				/* jump to the limiting step along the path */
				i_curr = findLimitingStep(path, smry.limits[i_curr].glb_inds[push_dir]);

				assert(smry.limits[i_curr].is_free!=NULL);

				not_fixed = smry.limits[i_curr].is_free[push_dir];
				lim_time = smry.limits[i_curr].times[push_dir];

				/* push the 'parent time' as little as possible */
				new_time = new_time + new_time * incr;

				/* see if the limit has to be pushed (might have been pushed further earlier) */
				need_to_push =
						smry.adjusted_times[i_curr] == DBL_MAX
								|| (push_dir == 0 && new_time < smry.adjusted_times[i_curr])
								|| (push_dir == 1 && new_time > smry.adjusted_times[i_curr]) ? 1 : 0;
				if (need_to_push == 0)
					break;

				is_limiter_tmp[i_curr] = 1;
				new_times[i_curr] = new_time;
				it++; // debug
				assert(it < 1000); // debug
			}

			ok =
					push_dir == 0 ? (lim_time < trgt_time ? 1 : 0) : (lim_time > trgt_time ? 1 : 0);

			if (ok == 0) {
				free(is_limiter_tmp);
				free(new_times);
				free(is_limiter);
				deallocatePath(path_in);
				return ok;
			}

			/* the actual adjustment of the times */
			for (int j = 0; j < smry.length; j++) {

				if (is_limiter_tmp[j] == 0)
					continue;

				new_node_idx = smry.new_node[j];
				new_node_gi = path.global_index[j][new_node_idx - nl];
				forwardTimeAdjustment(path, j, new_times[j], new_node_idx, new_node_gi);

				smry.adjusted_times[j] = new_times[j];

				updateLimiters(smry, j, new_times[j], nl, new_node_gi);

				is_limiter[j] =
						is_limiter[j] > is_limiter_tmp[j] ? is_limiter[j] : is_limiter_tmp[j];
			}
			free(is_limiter_tmp);
			free(new_times);
		}
	}

//	writeTikzTexFileForTreePath(infile, path_in.tree_path, path_in.opers, NULL, NULL,
//			path_in.global_index, path_in.path_len);

//	writeTikzTexFileForTreePath(intrmfile, path.tree_path, path.opers, NULL, NULL,
//			path.global_index, path.sites, path.path_len);

	/* Finally we make sure that the trees are consistent */
	for (int i = 1; i < path.path_len; i++) {
		if (checkNodeOrdering(path.tree_path[i]) == 0
				|| checkNodeTimeOrder(path.tree_path[i]) == 0) {
			consistency_check_done = 1;

			t_prev = path_in.tree_path[i - 1];
			t_curr = path_in.tree_path[i];

			new_node_idx = determineNewNode(t_prev, t_curr);
			// only a tree with a 'new node' can have been messed up by the operations above
			if (new_node_idx < 0) {
				printTree(t_prev);
				printTree(t_curr);
				assert(0 > 1);
			}
			new_node_gi = path_in.global_index[i][new_node_idx];

			limit_inds = findLimitingIndices(t_prev.times, t_curr.times[new_node_idx], nl - 1);
			limit_times[0] =
					limit_inds[0] < 0 ? 0 : path.tree_path[i - 1].times[limit_inds[0] - nl];
			limit_times[1] =
					limit_inds[1] < 0 ?
							limit_times[0] + 1 : path.tree_path[i - 1].times[limit_inds[1] - nl];
			free(limit_inds);

			u_bounders = postBounders(path, smry, i, 0, is_limiter);
			l_bounders = postBounders(path, smry, i, 1, is_limiter);

			upp_bound = fmin(limit_times[1], findMin(u_bounders));
			low_bound = fmax(limit_times[0], findMax(l_bounders));

			new_time = (upp_bound + low_bound) / (double) 2;
			forwardTimeAdjustment(path, i, new_time, new_node_idx + nl, new_node_gi);

			if (checkTree(path.tree_path[i]) != 1) {
				printLimitSummary(smry);
				printf("Initial limits %1.8f %1.8f\n", limit_times[0], limit_times[1]);
				printf("Upper:\n");
				printDoubleArray(u_bounders.v, 1, (int) u_bounders.length, 1);
				printf("Lower:\n");
				printDoubleArray(l_bounders.v, 1, (int) l_bounders.length, 1);
				printf("Limiter indicator:\n");
				printShortArray(is_limiter, 1, len, 1);
				writeTikzTexFileForTreePath(NULL, path.tree_path, path.opers, NULL, NULL,
						path.global_index, path.sites, path.rec_times, path.path_len);
				writeTikzTexFileForTreePath(infile, path_in.tree_path, path_in.opers, NULL, NULL,
						path_in.global_index, path_in.sites, path_in.rec_times, path_in.path_len);
				assert(false);
			}
			free(u_bounders.v);
			free(l_bounders.v);
		}
	}

//	writeTikzTexFileForTreePath(NULL, path.tree_path, path.opers, NULL, NULL,
//			path.global_index, path.path_len);

// final check
	for (int i = 0; i < path.path_len; i++) {
		if (checkTree(path.tree_path[i]) != 1) {
			writeTikzTexFileForTreePath(NULL, path.tree_path, path.opers, NULL, NULL,
					path.global_index, path.sites, path.rec_times, path.path_len);
			writeTikzTexFileForTreePath(infile, path_in.tree_path, path_in.opers, NULL, NULL,
					path_in.global_index, path_in.sites, path_in.rec_times, path_in.path_len);
			assert(false);
		}
	}

	deallocatePath(path_in);
	free(is_limiter);

	return ok;
}

struct DoubleVector postBounders(struct Smc path, struct LimitSummary smry, short step,
		short push_dir, short *is_limiter) {

	struct DoubleVector out;
	short *bounders, *processed, *newly_procd;
	short j, parent, count = 0, has_parent;
	bounders = malloc(sizeof(short) * path.path_len);
	processed = malloc(sizeof(short) * path.path_len);
	newly_procd = malloc(sizeof(short) * path.path_len);
	fillShortArray(processed, 0, 1, path.path_len, 1);
	fillShortArray(bounders, 0, 1, path.path_len, 1);
	fillShortArray(newly_procd, 0, 1, path.path_len, 1);

	/* iterate path from RIGHT to LEFT */
	for (int i = path.path_len - 1; i > step; i--) {

		if (processed[i] == 0) {

			fillShortArray(newly_procd, 0, 1, path.path_len, 1);
			j = i;
			newly_procd[j] = 1;

			has_parent = smry.limits[j].glb_inds == NULL ? 0 : 1;
			if (has_parent == 1)
				parent = findLimitingStep(path, smry.limits[j].glb_inds[push_dir]);
			else
				parent = -1;

			while (parent >= 0 && parent != step) {

				j = parent;
				newly_procd[j] = 1;

				has_parent = smry.limits[j].glb_inds == NULL ? 0 : 1;
				if (has_parent == 1)
					parent = findLimitingStep(path, smry.limits[j].glb_inds[push_dir]);
				else
					parent = -1;
			}

			if (has_parent == 1 && parent == step)
				for (int j = 0; j < path.path_len; j++)
					bounders[j] = newly_procd[j] == 1 ? 1 : bounders[j];

		}
		for (int j = 0; j < path.path_len; j++)
			processed[j] = processed[j] > newly_procd[j] ? processed[j] : newly_procd[j];
	}

	/* convert boolean indicator into index array */
	out.v = malloc(sizeof(double) * path.path_len);
	count = 0;
	for (int i = 0; i < path.path_len; i++)
		if (bounders[i] == 1 && is_limiter[i] == 1)
			out.v[count++] = smry.adjusted_times[i];
	out.length = count;
	out.v = realloc(out.v, sizeof(double) * count);

	free(bounders);
	free(processed);
	free(newly_procd);
	return out;
}

/* Find the step where the node with global index 'gi' appears for the first time.
 * Return -1, if not found at all. */
short findLimitingStep(struct Smc path, short gi) {

	short nb = (path.tree_path[0].n_nodes + 1) / 2 - 1;

	for (int i = 0; i < path.path_len; i++)
		for (int j = 0; j < nb; j++)
			if (path.global_index[i][j] == gi)
				return i;

	return -1;
}

void updateLimiters(struct LimitSummary smry, short step, double new_time, short nl,
		short new_node_gi) {

	for (int i = 0; i < smry.length; i++)
		if (smry.limits[i].glb_inds != NULL) {
			if (smry.limits[i].glb_inds[0] == new_node_gi) {
				smry.limits[i].times[0] = new_time;
			} else if (smry.limits[i].glb_inds[1] == new_node_gi) {
				smry.limits[i].times[1] = new_time;
			}
		}
}

void forwardTimeAdjustment(struct Smc path, short start_step, double new_time, short idx,
		short new_node_gi) {

	short nl = (path.tree_path[start_step].n_nodes + 1) / 2;
	int j;

	path.tree_path[start_step].times[idx - nl] = new_time;

	for (int i = start_step + 1; i < path.path_len; i++) {
		for (j = 0; j < nl - 1; j++)
			if (path.global_index[i][j] == new_node_gi) {
				path.tree_path[i].times[j] = new_time;
				break;
			}
		if (j == nl)
			// if adjustable node was not found, it won't be found in the
			// remaining trees either
			break;
	}
}

void printLimitSummary(struct LimitSummary summary) {

	short len = summary.length;

	printf("Limit summary\n-------------\n");
	printf("%15s", "Step:");
	for (int i = 0; i < len; i++)
		printf("%8hd ", summary.limits[i].step);
	printf("\n");
	printf("%15s", "Low lim:");
	for (int i = 0; i < len; i++)
		printf("%8hd ",
				(short) (summary.limits[i].inds == NULL ? -1 : summary.limits[i].inds[0]));
	printf("\n");
	printf("%15s", "Upp lim:");
	for (int i = 0; i < len; i++)
		printf("%8hd ",
				(short) (summary.limits[i].inds == NULL ? -1 : summary.limits[i].inds[1]));
	printf("\n");
	printf("%15s", "Free low:");
	for (int i = 0; i < len; i++)
		printf("%8hd ",
				(short) (summary.limits[i].is_free == NULL ? -1 : summary.limits[i].is_free[0]));
	printf("\n");
	printf("%15s", "Free upp:");
	for (int i = 0; i < len; i++)
		printf("%8hd ",
				(short) (summary.limits[i].is_free == NULL ? -1 : summary.limits[i].is_free[1]));
	printf("\n");
	printf("%15s", "Time low:");
	for (int i = 0; i < len; i++)
		printf("%8.2f ", summary.limits[i].times == NULL ? -1 : summary.limits[i].times[0]);
	printf("\n");
	printf("%15s", "Time upp:");
	for (int i = 0; i < len; i++)
		printf("%8.2f ", summary.limits[i].times == NULL ? -1 : summary.limits[i].times[1]);
	printf("\n");

	printf("%15s", "Is new:");
	for (int i = 0; i < len; i++)
		printf("%8hd ", summary.is_new[i]);
	printf("\n");
	printf("%15s", "Free:");
	for (int i = 0; i < len; i++)
		printf("%8hd ", summary.freedom_indicator[i]);
	printf("\n");
	printf("%15s", "Cond idx:");
	for (int i = 0; i < len; i++)
		printf("%8hd ", summary.conditioning_index[i]);
	printf("\n");
	printf("%15s", "Adj'd times:");
	for (int i = 0; i < len; i++)
		printf("%8.2f ", summary.adjusted_times[i]);
	printf("\n");
	printf("%15s", "New node idx:");
	for (int i = 0; i < len; i++)
		printf("%8hd ", summary.new_node[i]);
	printf("\n");

	printf("%15s", "GI low:");
	for (int i = 0; i < len; i++)
		printf("%8hd ",
				(short) (summary.limits[i].glb_inds == NULL ? -1 : summary.limits[i].glb_inds[0]));
	printf("\n");
	printf("%15s", "GI upp:");
	for (int i = 0; i < len; i++)
		printf("%8hd ",
				(short) (summary.limits[i].glb_inds == NULL ? -1 : summary.limits[i].glb_inds[1]));
	printf("\n");

	printf("\n");
	return;
}

short * findLimitingIndices(double *times, double new_time, short nb) {

	short *out;
	short nl = nb + 1;
	out = malloc(sizeof(short) * 2);

// find the largest smaller
	out[0] = -1;
	for (int i = 0; i < nb; i++) {
		if (times[i] >= new_time + new_time * TOL)
			break;
		out[0] = i + nl;
	}

	/* two special cases:
	 * 1) out[0] = -1, in which case out[1] = nl
	 * 2) out[0] = nb-1, in which case upper bound node does not exist */
	out[1] = out[0] == -1 ? nl : out[0] - nl == nb - 1 ? -1 : out[0] + 1;

	return out;
}

struct LimitData getLimiters(struct Tree t, double new_time, short step, struct Smc path) {

	struct LimitData lims;
	short nl = (t.n_nodes + 1) / 2;
	short nb = nl - 1;
	short *glb_idx;

	glb_idx = path.global_index[step];

	lims.is_free = malloc(sizeof(short) * 2);
	lims.times = malloc(sizeof(double) * 2);
	lims.glb_inds = malloc(sizeof(short) * 2);

	lims.inds = findLimitingIndices(t.times, new_time, nb);

//	printf("%hd %hd\n", lims.inds[0], lims.inds[1]);
//	printTree(t);
	lims.times[0] = lims.inds[0] < 0 ? 0 : t.times[lims.inds[0] - nl];
	lims.times[1] = lims.inds[1] < 0 ? DBL_MAX : t.times[lims.inds[1] - nl];

	lims.is_free[0] = lims.inds[0] < 0 ? 0 : is_unfixed(lims.inds[0] - nl, glb_idx, nb);
	lims.is_free[1] = lims.inds[1] < 0 ? 0 : is_unfixed(lims.inds[1] - nl, glb_idx, nb);

	lims.glb_inds[0] = lims.inds[0] < 0 ? -1 : path.global_index[step][lims.inds[0] - nl];
	lims.glb_inds[1] = lims.inds[1] < 0 ? -1 : path.global_index[step][lims.inds[1] - nl];

	lims.step = step;

	return lims;

}

/* If the global index is greater than the number of branch nodes -1 then node must be a
 * 'new node' indtroduces at some point and therefore not fixed
 */
short is_unfixed(short i, short *glb_idx, short nb) {
	return glb_idx[i] > nb - 1 ? 1 : 0;
}

short findConditioningIndex(short *gis, short gi, short nb) {

	for (int i = 0; i < nb; i++)
		if (gis[i] == gi)
			return i + nb + 1;
	return -1;

}

/* This function checks if a path contains extra recombinations */
short hasExtraRecombinations(struct Smc path, short *extra, struct Conditioning cond) {

	short k;

	for (int j = 0; j < cond.length - 1; j++) {
// iterate over extra recombination flags

		if (extra[j] == 0)
			continue;

// find the index of the tree at the destination site
// NB: destination site means the site of the tree that results from
// the recombination
		for (k = 0; k < path.path_len; k++)
			if (path.sites[k] == cond.sites[j + 1])
				break;

		assert(k > 0 && k < path.path_len);

// If extra combination is allowed and the operation is not
// no-op then the path has an extra recombination
		if (!(path.opers[k - 1][0] == -1 && path.opers[k - 1][1] == -1))
			return 1;
	}
	return 0;
}

