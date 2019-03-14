/*
 * backtracking.c
 *
 *  Created on: 16.11.2016
 *      Author: heine
 *
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
#include "constants.h"
#include "backtracking.h"
#include "graph2tikz.h"
#include "utils.h"
#include "timeadjustment.h"

void printLinkedSetArray(struct LinkedSetArray ls) {

	printf("Tree counts: ");
	for (int i = 0; i < ls.length; i++)
		printf("%ld ", ls.sets[i].n_trees);

}

struct SamplingSet simpleBacktrack(long idx, struct Tree init_tree, long init_tree_site,
		struct LinkedSetArray array, struct Data data) {

	struct SamplingSet sampling_set;
	struct LinkedSetArray extended_array;
	struct SimplePath path;
	short step;
	long path_count = 0, n_trees;

	sampling_set.paths = malloc(sizeof(struct Smc) * MAX_SAMPLING_SET_SIZE);
	sampling_set.set_selector = malloc(sizeof(short) * MAX_SAMPLING_SET_SIZE);
	sampling_set.n_paths = malloc(sizeof(long));
	sampling_set.set_cardinalities = malloc(sizeof(long) * 2);
	*sampling_set.n_paths = 0;

	// path is one entry longer than the linked set array as we will include the
	// initial tree in the path too
	path.length = array.length + 1;
	path.trees = malloc(sizeof(struct Tree) * path.length);
	path.sites = malloc(sizeof(long) * path.length);
	path.opers = malloc(sizeof(short *) * (path.length - 1));

	// create a search path starting from the initial tree
	for (int i = 0; i < path.length; i++)
		path.trees[i] = createTree(data.n_seq);
	copyTree(&path.trees[0], init_tree);

	for (int i = 0; i < path.length - 1; i++)
		path.opers[i] = malloc(sizeof(short) * 2);

	extended_array = includeInitialTree(array, init_tree, init_tree_site);
	step = extended_array.length - 1;
	n_trees = extended_array.sets[step].n_trees;

//	printf("Sites\n");
//	for (int i = 0; i < extended_array.length; i++)
//		printf("%ld %ld\n", extended_array.sets[i].sites[0], extended_array.sets[i].sites[1]);

	if (idx < 0) {
		// idx < 0 means that terminal tree is not specified and we do back tracking for
		// all terminal trees
		for (int i = 0; i < n_trees; i++) {
			copyTree(&path.trees[step], extended_array.sets[step].trees[i + n_trees]);
			path.sites[step] = extended_array.sets[step].sites[1];
			backwardStep(extended_array, step, i, path, &path_count, data, sampling_set);
		}
	} else {
		// do backtracking only for the specified index
		copyTree(&path.trees[step], extended_array.sets[step].trees[idx + n_trees]);
		path.sites[step] = extended_array.sets[step].sites[1];
		backwardStep(extended_array, step, idx, path, &path_count, data, sampling_set);
	}

	free(extended_array.sets[0].sites);
	deleteTree(extended_array.sets[0].trees[0]);
	deleteTree(extended_array.sets[0].trees[1]);
	free(extended_array.sets[0].trees);
	free(extended_array.sets);

	for (int i = 0; i < path.length; i++)
		deleteTree(path.trees[i]);
	free(path.trees);
	free(path.sites);
	for (int i = 0; i < path.length - 1; i++)
		free(path.opers[i]);
	free(path.opers);

	// trim the sampling set to actual size
	sampling_set.paths = realloc(sampling_set.paths, sizeof(struct Smc) * path_count);

	return sampling_set;

}

struct LinkedSetArray includeInitialTree(struct LinkedSetArray array,
		struct Tree init_tree, long init_tree_site) {

	struct LinkedSetArray out;
	struct LinkedSet init;

	out.sets = malloc(sizeof(struct LinkedSet) * (array.length + 1));
	out.length = array.length + 1;

	init.length = 2;
	init.links = NULL;
	init.n_trees = 1;
	init.oper = NULL;
	init.set_size = NULL;
	init.sites = malloc(sizeof(long) * 2);
	init.sites[0] = -1;
	init.sites[1] = init_tree_site;
	init.trees = malloc(sizeof(struct Tree) * 2);
	init.trees[0] = createCopy(init_tree);
	init.trees[1] = createCopy(init_tree);

	out.sets[0] = init;
	for (int i = 0; i < array.length; i++)
		out.sets[i + 1] = array.sets[i];

	return out;
}

void backwardStep(struct LinkedSetArray array, short step, long tree_idx,
		struct SimplePath path, long *path_count, struct Data data, struct SamplingSet ss) {

	struct Smc new_sampling_set_entry;
	long n_parents, i_prnt, next_step, n_trees;
	short n_global = 0;

	// termination
	if (step == 0) {

//		checkBacktrackRecursion(path, data, array);
		new_sampling_set_entry.path_len = path.length;

		new_sampling_set_entry.tree_path = malloc(sizeof(struct Tree) * path.length);
		new_sampling_set_entry.sites = malloc(sizeof(int) * path.length);
		for (int i = 0; i < path.length; i++) {
			new_sampling_set_entry.tree_path[i] = createCopy(path.trees[i]);
			new_sampling_set_entry.sites[i] = (int) path.sites[i];
		}

		new_sampling_set_entry.opers = malloc(sizeof(short *) * (path.length - 1));
		for (int i = 0; i < path.length - 1; i++) {
			new_sampling_set_entry.opers[i] = malloc(sizeof(short) * 2);
			new_sampling_set_entry.opers[i][0] = path.opers[i][0];
			new_sampling_set_entry.opers[i][1] = path.opers[i][1];
		}
		new_sampling_set_entry.rec_times = malloc(sizeof(double) * (path.length - 1));
		fillDoubleArray(new_sampling_set_entry.rec_times, -1, 1, path.length - 1, 1);
		new_sampling_set_entry.tree_selector = NULL;
		new_sampling_set_entry.global_index = malloc(sizeof(short*) * path.length);
		for (int i = 0; i < path.length; i++)
			assignGlobalIndices(new_sampling_set_entry.global_index, &n_global, i,
					new_sampling_set_entry.tree_path);
		new_sampling_set_entry.is_free = malloc(sizeof(short) * path.length);

		ss.paths[*path_count] = new_sampling_set_entry;
		*path_count = *path_count + 1;
		*ss.n_paths = *path_count;

	} else { // deepen recursion

		n_parents = array.sets[step].set_size[tree_idx];
		next_step = step - 1;
		for (int i = 0; i < n_parents; i++) {
			i_prnt = array.sets[step].links[tree_idx][i];
			n_trees = array.sets[next_step].n_trees;
			assert(i_prnt < n_trees);

			// insert tree new tree in the path
			copyTree(&path.trees[next_step], array.sets[next_step].trees[i_prnt + n_trees]);

			// insert new operation to the path
			path.opers[next_step][0] = array.sets[step].oper[tree_idx][i];
			path.opers[next_step][1] = array.sets[step].oper[tree_idx][i + n_parents];

			// deepen the recursion
			path.sites[next_step] = array.sets[next_step].sites[1];
			backwardStep(array, step - 1, i_prnt, path, path_count, data, ss);
		}

	}

}

void checkBacktrackRecursion(struct SimplePath path, struct Data data,
		struct LinkedSetArray array) {

	enum Node_Colors **colors;
	short n_nodes = path.trees[0].n_nodes;
	short n_leaves = (n_nodes + 1) / 2;

	if (checkPath(path, data) != 1) {

		colors = malloc(sizeof(enum Node_Colors *) * array.length);
		for (int i = 0; i < array.length; i++) {
			colors[i] = malloc(sizeof(enum Node_Colors) * n_nodes);
			for (int k = 0; k < n_leaves; k++)
				colors[i][k] = data.M[k + i * n_leaves] == 1 ? black : white;
			for (int k = 0; k < n_leaves - 1; k++)
				colors[i][k + n_leaves] = white;
		}

		writeTikzTexFileForTreePath(NULL, path.trees, NULL, colors, NULL, NULL, NULL, NULL,
				array.length);
		for (int i = 0; i < array.length; i++)
			free(colors[i]);
		free(colors);

		assert(1 < 0);
	}
}

