/*
 * initialisation.c
 *
 *  Created on: 5.12.2016
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
#include <stdlib.h>
#include <stdio.h>
#include "ARGdecomposer.h"
#include "data.h"
#include "datastructures.h"
#include "debugging.h"
#include "fileaccess.h"
#include "graph2tikz.h"
#include "initialisation.h"
#include "pathutils.h"
#include "recombinationTimes.h"
#include "shrub.h"
#include "timeadjustment.h"
#include "treeutils.h"

void writeInitialisationTikzTexFile(struct Tree *trees, short **opers, struct Data data,
		int *op_sites, short length) {

	enum Node_Colors **colors;
	char filepath[] = "/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/init.tex"; // Fix the path

	colors = malloc(sizeof(enum Node_Colors *) * length);

	for (int i = 0; i < length; i++)
		if (intIncluded(op_sites[i], data.segregating_sites, data.n_sites) == 1)
			colors[i] = colorTree(trees[i], data.M, getDataColumn(data, op_sites[i]));
		else {
			colors[i] = malloc(sizeof(enum Node_Colors) * trees[i].n_nodes);
			for (int j = 0; j < trees[i].n_nodes; j++)
				colors[i][j] = white;
		}

	writeTikzTexFileForTreePath(filepath, trees, opers, colors, NULL, NULL, NULL, NULL,
			length);

	for (int i = 0; i < length; i++)
		free(colors[i]);
	free(colors);
}

struct Smc initialisation(struct Data data, struct Parameters parm) {

	struct ShortVector selector;
	struct Smc path;
	struct ShrubARG arg;
	struct Tree *trees;
	short **opers;
	int *op_sites;
	short n_branch, length, n_global = 0;

	// Initialise the chain with the SHRUB type algorithm
	arg = shrub(data);
//	printARG(arg);
	n_branch = arg.n_nodes - arg.n_seq;
	length = arg.n_rec + 1;

	trees = malloc(sizeof(struct Tree) * length);
	opers = malloc(sizeof(short *) * (length - 1));
	for (int i = 0; i < (length - 1); i++)
		opers[i] = malloc(sizeof(short) * 3);
	op_sites = malloc(sizeof(int) * length);

	/* Decompose the initial ARG into a coalescent tree sequence */
	decomposeARG(arg, trees, opers, data, n_branch, op_sites);

	path.opers = opers;
	path.path_len = length;
	path.rec_times = malloc(sizeof(double) * (length - 1));
	path.sites = op_sites;
	path.tree_path = trees;
	path.tree_selector = NULL;
	path.selector_length = 0;
	path.is_free = NULL;

	sortChildren(path.tree_path, path.path_len);

	/* ensure that the path has global indexing */
	path.global_index = malloc(sizeof(short *) * path.path_len);
	for (int i = 0; i < path.path_len; i++)
		assignGlobalIndices(path.global_index, &n_global, i, path.tree_path);

	path = generateTimes(path, parm);
	recombinationTimes(path, 1);
	assert(checkRecombinationTimes(path)==1);

	selector = createTreeSelector(path, data);
	path.tree_selector = selector.v;
	path.selector_length = (int) selector.length;

	/* delete diagnostics file */
	removeDiagnosticsFile();

	free(arg.C);
	return path;
}

// This function is simply for debugging purposes
void manualInitialPathTimes(struct Smc path) {

	// Coalescence times
	path.tree_path[0].times[0] = 0.362767047375791e04;
	path.tree_path[0].times[1] = 0.448543112159897e04;
	path.tree_path[0].times[2] = 0.449506478908682e04;
	path.tree_path[0].times[3] = 0.455735711198971e04;
	path.tree_path[0].times[4] = 0.667206250588489e04;
	path.tree_path[0].times[5] = 1.361301330113255e04;
	path.tree_path[0].times[6] = 1.653037808501839e04;
	path.tree_path[0].times[7] = 3.983450341571479e04;

	path.tree_path[1].times[0] = 0.362767047375791e04;
	path.tree_path[1].times[1] = 0.394841166658902e04;
	path.tree_path[1].times[2] = 0.449506478908682e04;
	path.tree_path[1].times[3] = 0.455735711198971e04;
	path.tree_path[1].times[4] = 0.667206250588489e04;
	path.tree_path[1].times[5] = 1.361301330113255e04;
	path.tree_path[1].times[6] = 1.653037808501839e04;
	path.tree_path[1].times[7] = 3.983450341571479e04;

	path.tree_path[2].times[0] = 0.394841166658902e04;
	path.tree_path[2].times[1] = 0.434429093235772e04;
	path.tree_path[2].times[2] = 0.449506478908682e04;
	path.tree_path[2].times[3] = 0.455735711198971e04;
	path.tree_path[2].times[4] = 0.667206250588489e04;
	path.tree_path[2].times[5] = 1.361301330113255e04;
	path.tree_path[2].times[6] = 1.653037808501839e04;
	path.tree_path[2].times[7] = 3.983450341571479e04;

	path.tree_path[3].times[0] = 0.332316855196124e04;
	path.tree_path[3].times[1] = 0.394841166658902e04;
	path.tree_path[3].times[2] = 0.449506478908682e04;
	path.tree_path[3].times[3] = 0.455735711198971e04;
	path.tree_path[3].times[4] = 0.667206250588489e04;
	path.tree_path[3].times[5] = 1.361301330113255e04;
	path.tree_path[3].times[6] = 1.653037808501839e04;
	path.tree_path[3].times[7] = 3.983450341571479e04;

	path.tree_path[4].times[0] = 0.072260648346975e04;
	path.tree_path[4].times[1] = 0.332316855196124e04;
	path.tree_path[4].times[2] = 0.449506478908682e04;
	path.tree_path[4].times[3] = 0.455735711198971e04;
	path.tree_path[4].times[4] = 0.667206250588489e04;
	path.tree_path[4].times[5] = 1.361301330113255e04;
	path.tree_path[4].times[6] = 1.653037808501839e04;
	path.tree_path[4].times[7] = 3.983450341571479e04;

	path.tree_path[5].times[0] = 0.332316855196124e04;
	path.tree_path[5].times[1] = 0.377241023823807e04;
	path.tree_path[5].times[2] = 0.449506478908682e04;
	path.tree_path[5].times[3] = 0.455735711198971e04;
	path.tree_path[5].times[4] = 0.667206250588489e04;
	path.tree_path[5].times[5] = 1.361301330113255e04;
	path.tree_path[5].times[6] = 1.653037808501839e04;
	path.tree_path[5].times[7] = 3.983450341571479e04;

	path.tree_path[6].times[0] = 0.136100333544548e04;
	path.tree_path[6].times[1] = 0.332316855196124e04;
	path.tree_path[6].times[2] = 0.377241023823807e04;
	path.tree_path[6].times[3] = 0.449506478908682e04;
	path.tree_path[6].times[4] = 0.455735711198971e04;
	path.tree_path[6].times[5] = 1.361301330113255e04;
	path.tree_path[6].times[6] = 1.653037808501839e04;
	path.tree_path[6].times[7] = 3.983450341571479e04;

	path.tree_path[7].times[0] = 0.332316855196124e04;
	path.tree_path[7].times[1] = 0.377241023823807e04;
	path.tree_path[7].times[2] = 0.449506478908682e04;
	path.tree_path[7].times[3] = 0.455735711198971e04;
	path.tree_path[7].times[4] = 0.605969575685963e04;
	path.tree_path[7].times[5] = 1.361301330113255e04;
	path.tree_path[7].times[6] = 1.653037808501839e04;
	path.tree_path[7].times[7] = 3.983450341571479e04;

	// Recombination times
	path.rec_times[0] = 1.976843871131627e03;
	path.rec_times[1] = 1.697994380839845e03;
	path.rec_times[2] = 1.664184740585989e03;
	path.rec_times[3] = 0.107797109547817e03;
	path.rec_times[4] = 0.592738154670294e03;
	path.rec_times[5] = 0.437377518912925e03;
	path.rec_times[6] = 1.060458408955579e03;

}

