/*
 * pathutils.h
 *
 *  Created on: 17.11.2016
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

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "constants.h"
#include "data.h"
#include "datastructures.h"
#include "debugging.h"
#include "graph2tikz.h"
#include "MCMCutils.h"
#include "pathutils.h"
#include "randomness.h"
#include "sorting.h"
#include "timeadjustment.h"
#include "treeutils.h"
#include "utils.h"

short checkGlobalOrderConsistency(struct Smc path1, struct Smc path2) {

	struct Smc pths[2];
	short gi_max[] = { -1, -1 };
	double *times[2];
	int *i_sorts[2];
	short nb;

	nb = (path1.tree_path[0].n_nodes + 1) / 2 - 1;

	pths[0] = path1;
	pths[1] = path2;

	times[0] = malloc(sizeof(double) * (pths[0].path_len + nb));
	times[1] = malloc(sizeof(double) * (pths[0].path_len + nb));
	i_sorts[0] = malloc(sizeof(int) * (pths[0].path_len + nb));
	i_sorts[1] = malloc(sizeof(int) * (pths[0].path_len + nb));

	assert(path1.global_index!=NULL);
	assert(path2.global_index!=NULL);

	for (int k = 0; k < 2; k++)
		for (int i = 0; i < path1.path_len; i++)
			for (int j = 0; j < nb; j++)
				if (pths[k].global_index[i][j] > gi_max[k]) {
					gi_max[k] = pths[k].global_index[i][j];
					i_sorts[k][gi_max[k]] = gi_max[k];
					times[k][gi_max[k]] = pths[k].tree_path[i].times[j];
				}

	for (int i = 0; i <= gi_max[0]; i++)
		printf("%i %8.10f %i %8.10f\n", i_sorts[0][i], times[0][i], i_sorts[1][i],
				times[1][i]);

	sortDoubleArray(i_sorts[0], times[0], gi_max[0], 0);
	sortDoubleArray(i_sorts[1], times[1], gi_max[1], 0);

	for (int i = 0; i <= gi_max[0]; i++)
		printf("%3i", i_sorts[0][i]);
	printf("\n");
	for (int i = 0; i <= gi_max[0]; i++)
		printf("%3i", i_sorts[1][i]);
	printf("\n");

	free(times[0]);
	free(times[1]);

	for (int i = 0; i < gi_max[0]; i++)
		if (i_sorts[0][i] != i_sorts[1][i]) {
			free(i_sorts[0]);
			free(i_sorts[1]);
			return 0;
		}

	free(i_sorts[0]);
	free(i_sorts[1]);
	return 1;
}

struct Smc removeNoOps(struct Smc path) {

	struct Smc out;
	short n_rec = 0, j, new_length, n_global = 0;

	// count recombinations
	for (int i = 0; i < path.path_len - 1; i++)
		n_rec += path.opers[i][0] == -1 ? 0 : 1;
	new_length = n_rec + 1;

	out.tree_path = malloc(sizeof(struct Tree) * new_length);
	out.opers = malloc(sizeof(short*) * n_rec);
	out.sites = malloc(sizeof(int) * new_length);
	out.rec_times = malloc(sizeof(double) * n_rec);
	out.global_index = malloc(sizeof(short*) * new_length);
	out.is_free = malloc(sizeof(short) * new_length);
	out.path_len = new_length;
	out.selector_length = 0;
	out.tree_selector = NULL;

	j = 0;
	for (int i = 0; i < new_length - 1; i++) {

		out.sites[i] = path.sites[j];

		// skip over no-ops
		while (path.opers[j][0] == -1) {
//			deleteTree(path.tree_path[j]);
//			free(path.opers[j]);
			j++;
		}

//		out.tree_path[i] = path.tree_path[j];
		out.tree_path[i] = createCopy(path.tree_path[j]);
//		out.opers[i] = path.opers[j];
		out.opers[i] = createOperationCopy(path.opers[j]);
		out.rec_times[i] = path.rec_times[j];
		j++;
	}
//	out.tree_path[new_length - 1] = path.tree_path[j];
	out.tree_path[new_length - 1] = createCopy(path.tree_path[j]);
	out.sites[new_length - 1] = path.sites[j];

	for (int i = 0; i < out.path_len; i++)
		assignGlobalIndices(out.global_index, &n_global, i, out.tree_path);

	deallocatePath(path);
//	for (int i = j + 1; i < path.path_len; i++)
//		deleteTree(path.tree_path[i]);
//	for (int i = j + 1; i < path.path_len - 1; i++)
//		free(path.opers[i]);
//	free(path.tree_path);
//	free(path.opers);
//	free(path.sites);
//	free(path.is_free);
//	free(path.tree_selector);
//	free(path.rec_times);
//	for (int i = 0; i < path.path_len; i++)
//		free(path.global_index[i]);
//	free(path.global_index);

	assert(checkTreePathCompletely(out)==1);
	return out;
}

struct Smc extractToSites(struct Smc path, struct LongVector sites) {

	struct Smc out;
	struct LongVector site_union;
	short n_global = 0, i_tree, i_tree_next;

	site_union = siteUnion(sites, path);

	out.tree_path = malloc(sizeof(struct Tree) * site_union.length);
	out.opers = malloc(sizeof(short *) * (site_union.length - 1));
	out.sites = malloc(sizeof(int) * site_union.length);
	out.rec_times = malloc(sizeof(double) * (site_union.length - 1));
	out.is_free = malloc(sizeof(short) * site_union.length);
	out.global_index = malloc(sizeof(short *) * site_union.length);
	out.path_len = site_union.length;

	fillShortArray(out.is_free, -1, 1, (int) site_union.length, 1);

	/* trees and sites */
	for (int i = 0; i < site_union.length; i++) {
		i_tree = path.tree_selector[site_union.values[i] - 1];
		out.tree_path[i] = createCopy(path.tree_path[i_tree]);
		out.sites[i] = (int) site_union.values[i];
	}

	/* operations and recombination times */
	for (int i = 0; i < site_union.length - 1; i++) {

		i_tree = path.tree_selector[site_union.values[i] - 1];
		i_tree_next = path.tree_selector[site_union.values[i + 1] - 1];

		if (i_tree == i_tree_next) {
			out.opers[i] = createNoop();
			out.rec_times[i] = -1;
		} else {
			out.opers[i] = createOperationCopy(path.opers[i_tree]);
			out.rec_times[i] = path.rec_times[i_tree];
		}
	}
	free(site_union.values);

	for (int i = 0; i < out.path_len; i++)
		assignGlobalIndices(out.global_index, &n_global, i, out.tree_path);

	out.selector_length = 0;
	out.tree_selector = NULL;

	return out;
}

struct LongVector siteUnion(struct LongVector sites, struct Smc path) {

	short j = 0, cnt = 0;
	struct LongVector out;
	out.values = malloc(sizeof(long) * (sites.length + path.path_len));

	for (int i = 0; i < sites.length - 1; i++) {
		out.values[cnt++] = sites.values[i];
		while (j < path.path_len && path.sites[j] < sites.values[i + 1]) {
			if (path.sites[j] > sites.values[i] && path.sites[j] < sites.values[i + 1])
				out.values[cnt++] = path.sites[j];
			j++;
		}

	}
	out.values[cnt++] = sites.values[sites.length - 1];

	out.values = realloc(out.values, sizeof(long) * cnt);
	out.length = cnt;

	return out;
}

short *createNoop() {
	short *out;
	out = malloc(sizeof(short) * 2);
	out[0] = -1;
	out[1] = -1;
	return out;
}

struct Smc createPath(long length) {

	struct Smc out;

	out.tree_path = malloc(sizeof(struct Tree) * length);
	out.opers = malloc(sizeof(short *) * (length - 1));
	out.global_index = malloc(sizeof(short *) * length);
	out.is_free = malloc(sizeof(short) * length);
	out.rec_times = malloc(sizeof(double) * (length - 1));
	out.sites = malloc(sizeof(int) * length);
	out.tree_selector = NULL;

	return out;
}

struct Smc createPathCopy(const struct Smc path) {

	short nb = (path.tree_path[0].n_nodes + 1) / 2 - 1;
	struct Smc copy;

	copy.path_len = path.path_len;
	copy.selector_length = path.selector_length;

	/* trees */
	copy.tree_path = malloc(sizeof(struct Tree) * path.path_len);
	for (int i = 0; i < path.path_len; i++)
		copy.tree_path[i] = createCopy(path.tree_path[i]);

	/* global indices */
	if (path.global_index != NULL) {
		copy.global_index = malloc(sizeof(short *) * path.path_len);
		for (int i = 0; i < path.path_len; i++) {
			copy.global_index[i] = malloc(sizeof(short) * nb);
			for (int j = 0; j < nb; j++)
				copy.global_index[i][j] = path.global_index[i][j];
		}
	} else {
		copy.global_index = NULL;
	}

	/* opers */
	if (path.opers != NULL) {
		copy.opers = malloc(sizeof(short *) * (path.path_len - 1));
		for (int i = 0; i < path.path_len - 1; i++) {
			copy.opers[i] = malloc(sizeof(short) * 2);
			copy.opers[i][0] = path.opers[i][0];
			copy.opers[i][1] = path.opers[i][1];
		}
	} else {
		copy.opers = NULL;
	}

	/* recombination times */
	if (path.rec_times != NULL) {
		copy.rec_times = malloc(sizeof(double) * (path.path_len - 1));
		for (int i = 0; i < path.path_len - 1; i++)
			copy.rec_times[i] = path.rec_times[i];
	} else {
		copy.rec_times = NULL;
	}
	/* recombination sites */
	if (path.sites != NULL) {
		copy.sites = malloc(sizeof(int) * path.path_len);
		for (int i = 0; i < path.path_len; i++)
			copy.sites[i] = path.sites[i];
	} else {
		copy.sites = NULL;
	}

	/* tree selector */
	if (path.tree_selector != NULL) {
		copy.tree_selector = malloc(sizeof(short) * path.selector_length);
		for (int i = 0; i < copy.selector_length; i++)
			copy.tree_selector[i] = path.tree_selector[i];
	} else {
		copy.tree_selector = NULL;
	}

	if (path.is_free != NULL) {
		copy.is_free = malloc(sizeof(short) * path.path_len);
		for (int i = 0; i < path.path_len; i++)
			copy.is_free[i] = path.is_free[i];
	} else {
		copy.is_free = NULL;
	}
	return copy;
}
struct Smc createPathCopy2(const struct Smc path) {

	short nb = (path.tree_path[0].n_nodes + 1) / 2 - 1;
	struct Smc copy;

	copy.path_len = path.path_len;
	copy.selector_length = path.selector_length;

	/* trees */
	copy.tree_path = malloc(sizeof(struct Tree) * path.path_len);
	for (int i = 0; i < path.path_len; i++)
		copy.tree_path[i] = createCopy(path.tree_path[i]);

	/* global indices */
	if (path.global_index != NULL) {
		copy.global_index = malloc(sizeof(short *) * path.path_len);
		for (int i = 0; i < path.path_len; i++) {
			copy.global_index[i] = malloc(sizeof(short) * nb);
			for (int j = 0; j < nb; j++)
				copy.global_index[i][j] = path.global_index[i][j];
		}
	} else {
		copy.global_index = NULL;
	}

	/* opers */
	if (path.opers != NULL) {
		copy.opers = malloc(sizeof(short *) * (path.path_len - 1));
		for (int i = 0; i < path.path_len - 1; i++) {
			copy.opers[i] = malloc(sizeof(short) * 2);
			copy.opers[i][0] = path.opers[i][0];
			copy.opers[i][1] = path.opers[i][1];
		}
	} else {
		copy.opers = NULL;
	}

	/* recombination times */
	if (path.rec_times != NULL) {
		copy.rec_times = malloc(sizeof(double) * (path.path_len - 1));
		for (int i = 0; i < path.path_len - 1; i++)
			copy.rec_times[i] = path.rec_times[i];
	} else {
		copy.rec_times = NULL;
	}
	/* recombination sites */
	if (path.sites != NULL) {
		copy.sites = malloc(sizeof(int) * path.path_len);
		for (int i = 0; i < path.path_len; i++)
			copy.sites[i] = path.sites[i];
	} else {
		copy.sites = NULL;
	}

	/* tree selector */
	if (path.tree_selector != NULL) {
		copy.tree_selector = malloc(sizeof(short) * path.selector_length);
		for (int i = 0; i < copy.selector_length; i++)
			copy.tree_selector[i] = path.tree_selector[i];
	} else {
		copy.tree_selector = NULL;
	}

	if (path.is_free != NULL) {
		copy.is_free = malloc(sizeof(short) * path.path_len);
		for (int i = 0; i < path.path_len; i++)
			copy.is_free[i] = path.is_free[i];
	} else {
		copy.is_free = NULL;
	}
	return copy;
}

long *segmentEndIndices(struct Smc path, struct Conditioning cond) {

	long i_start, i_end;
	long *out;
	int i_left, i_right;
	out = malloc(sizeof(long) * 2);

	if (cond.backward == 1) {
		i_left = cond.sites[cond.length - 1];
		i_right = cond.sites[0];
	} else {
		i_left = cond.sites[0];
		i_right = cond.sites[cond.length - 1];
	}

	/* find the greatest site index smaller than or equal to first conditioning site */
	i_start = 0;
	for (long i = 0; i < path.path_len; i++)
		if (path.sites[i] > i_left)
			break;
		else
			i_start = i;

	/* find the greatest site index smaller than or equal to the last conditioning site */
	i_end = 0;
	for (long i = i_start; i < path.path_len; i++)
		if (path.sites[i] > i_right)
			break;
		else
			i_end = i;

	out[0] = i_start;
	out[1] = i_end;

	return out;
}

void printRecombinationTimes(struct Smc path) {
	printf("Recombination times: ");
	for (int i = 0; i < path.path_len - 1; i++)
		printf("%8.4f ", path.rec_times[i]);
	printf("\n");
}

struct SimplePath createNewPath(struct Smc path) {

	struct SimplePath new_path;
	short len = path.path_len;

	new_path.trees = malloc(sizeof(struct Tree) * len);
	for (int i = 0; i < len; i++)
		new_path.trees[i] = createCopy(path.tree_path[i]);
	new_path.opers = malloc(sizeof(short *) * (len - 1));
	for (int i = 0; i < len - 1; i++) {
		new_path.opers[i] = malloc(sizeof(short) * 2);
		new_path.opers[i][0] = path.opers[i][0];
		new_path.opers[i][1] = path.opers[i][1];
	}
	new_path.length = 1;

	assert(path.sites!=NULL);

	new_path.sites = malloc(sizeof(long) * len);
	for (int i = 0; i < len; i++)
		new_path.sites[i] = path.sites[i];

	return new_path;
}

void printGlobalIndices(struct Smc path) {

	short nb = (path.tree_path[0].n_nodes + 1) / 2 - 1;
	printf("Global indices\n");
	for (int j = 0; j < nb; j++) {
		for (int i = 0; i < path.path_len; i++)
			printf("%4hd ", path.global_index[i][j]);
		printf("\n");
	}
}

void printTreeSelector(struct Smc path) {
	for (int i = 0; i < path.selector_length; i++)
		printf("%hd ", path.tree_selector[i]);
	printf("\n");
}

void printOperations(struct Smc path) {
	printf("\n");
	for (int i = 0; i < path.path_len - 1; i++)
		printf("(%3hd,%3hd)\n", path.opers[i][0], path.opers[i][1]);
}

short countRecombinations(struct Smc path) {
	short n_rec = 0;
	for (int i = 0; i < path.path_len - 1; i++)
		if (path.opers[i][0] != -1)
			n_rec++;
	return n_rec;
}

enum Node_Colors **colorPath(struct Smc path, struct Data data) {

	enum Node_Colors **colors;
	short col;

	colors = malloc(sizeof(enum Node_Colors*) * path.path_len);
	for (int i = 0; i < path.path_len; i++) {
		col = getDataColumn(data, path.sites[i]);
		colors[i] = colorTree(path.tree_path[i], data.M, col);
	}
	return colors;
}

struct MRCA timesToMRCA(struct Smc path) {

	struct MRCA out;
	short nl = (path.tree_path[0].n_nodes + 1) / 2, j2;
	struct ShortVector ancs1, ancs2;
	double time_to_MRCA;

	/* create the time array */
	out.times = malloc(sizeof(double*) * (nl - 1));
	for (int i = 0; i < nl - 1; i++){
		out.times[i] = malloc(sizeof(double) * (nl - 1 - i));
		fillDoubleArray(out.times[i],DBL_MAX,1,(nl - 1 - i),1);
	}
	out.nl = nl;

	for (int k = 0; k < path.path_len; k++) {
		// iterate over leaf pairs
		for (int i = 0; i < nl - 1; i++)
			for (int j = i + 1; j < nl; j++) {
				ancs1 = findAncestors(path.tree_path[k], i);
				ancs2 = findAncestors(path.tree_path[k], j);
				time_to_MRCA = path.tree_path[k].times[firstCommon(ancs1, ancs2)];
				j2 = j - i - 1;
				out.times[i][j2] =
						out.times[i][j2] < time_to_MRCA ? out.times[i][j2] : time_to_MRCA;
				free(ancs1.v);
				free(ancs2.v);
			}
	}

	return out;
}

void deallocateMRCA(struct MRCA mrca) {
	for(int i = 0; i < (mrca.nl-1);i++)
		free(mrca.times[i]);
	free(mrca.times);
}

void deallocatePath(struct Smc path) {

// free trees and operations
	for (int i = 0; i < path.path_len - 1; i++) {
		deleteTree(path.tree_path[i]);
		free(path.opers[i]);
	}
	free(path.tree_path[path.path_len - 1].C);
	free(path.tree_path[path.path_len - 1].times);
	free(path.tree_path);
	free(path.opers);

// global indices
	if (path.global_index != NULL) {
		for (int i = 0; i < path.path_len; i++)
			free(path.global_index[i]);
		free(path.global_index);
	}

	if (path.tree_selector != NULL)
		free(path.tree_selector);

	if (path.rec_times != NULL)
		free(path.rec_times);

	if (path.sites != NULL)
		free(path.sites);

	if (path.is_free != NULL)
		free(path.is_free);

}

struct Smc generateTimes(struct Smc path, struct Parameters parm) {

	struct Smc new_path;
	struct SamplingSet recreated;
	short nb = (path.tree_path[0].n_nodes + 1) / 2 - 1;

	/* replace the times in the first tree by the randomly sampled ones */
	free(path.tree_path[0].times);
	path.tree_path[0].times = createCoalescenceTimes(nb, parm);

	/* recreate the path using the operators to update all times */
	recreated = recreatePath(path);

	new_path = createPathCopy(recreated.paths[0]);

	deallocatePath(path);
	deallocateSamplingSet(recreated);

	return new_path;
}

double *createCoalescenceTimes(short nb, struct Parameters parm) {

	double *times, *cumtimes;

	times = malloc(sizeof(double) * nb);
	cumtimes = malloc(sizeof(double) * nb);
	for (int i = 0; i < nb; i++)
		times[nb - 1 - i] = exprnd((double) (2 * parm.n_eff) / (double) nchoosek(i + 2, 2));
	cumtimes[0] = times[0];
	for (int i = 1; i < nb; i++)
		cumtimes[i] = cumtimes[i - 1] + times[i];

	free(times);
	return cumtimes;
}

double *calculateRegraftTimes(struct Tree *t, short path_len) {

	double *regraft_times;
	struct Tree t_curr, t_next;
	short nni;

	regraft_times = malloc(sizeof(double) * path_len);
	fillDoubleArray(regraft_times, -1.0, 1, path_len, 1);

	for (int i = 0; i < path_len - 1; i++) {

		t_curr = t[i];
		t_next = t[i + 1];

		if ((nni = determineNewNode(t_curr, t_next)) >= 0)
			regraft_times[i] = t_next.times[nni];
	}

	return regraft_times;
}
