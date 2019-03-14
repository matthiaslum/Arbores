/*
 * debugging.c
 *
 *  Created on: 14.11.2016
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
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "datastructures.h"
#include "debugging.h"
#include "fileaccess.h"
#include "graph2tikz.h"
#include "pathutils.h"
#include "randomness.h"
#include "sorting.h"
#include "treeutils.h"

short checkNodeOrdering(struct Tree t) {

	short nl = (t.n_nodes + 1) / 2;
	short nb = nl - 1;
	short out = 1;

	for (int i = 0; i < nb; i++) {
		out = out * (t.C[i] < i + nl ? 1 : 0);
		out = out * (t.C[i + nb] < i + nl ? 1 : 0);
	}
	return out;
}

short checkNodeTimeOrder(struct Tree t) {

	short nl = (t.n_nodes + 1) / 2;
	short nb = nl - 1;
	short out = 1;

	for (int i = 0; i < nb - 1; i++)
		out *= t.times[i + 1] > t.times[i] ? 1 : 0;

	return out;
}

short checkTopology(struct Tree t) {

	short nb = (t.n_nodes + 1) / 2 - 1;
	short out = 1;
	int *indices, *array;

	array = malloc(sizeof(int) * 2 * nb);
	indices = malloc(sizeof(int) * 2 * nb);
	for (int i = 0; i < 2 * nb; i++) {
		array[i] = (int) t.C[i];
		indices[i] = i;
	}
	sortIntArray(indices, array, 2 * nb);

	for (int i = 0; i < 2 * nb - 1; i++)
		out *= array[indices[i + 1]] > array[indices[i]] ? 1 : 0;

	for (int i = 0; i < nb; i++)
		out *= t.C[i] < t.C[i + nb] ? 1 : 0;

	free(array);
	free(indices);
	return out;
}

short checkTreePathStructure(const struct Smc path) {
	short ok = 1;
	for (int i = 0; i < path.path_len; i++)
		ok *= checkNodeOrdering(path.tree_path[i]) * checkTopology(path.tree_path[i]);
	return ok;
}

short checkTreePathCompletely(const struct Smc path) {
	short ok = 1;
	for (int i = 0; i < path.path_len; i++)
		ok *= checkNodeOrdering(path.tree_path[i]) * checkTopology(path.tree_path[i])
				* checkNodeTimeOrder(path.tree_path[i]);
	return ok;
}

short checkTree(struct Tree t) {
	return checkNodeOrdering(t) * checkTopology(t) * checkNodeTimeOrder(t);
}

short checkOperations(struct Smc path) {

//	char dbfile[] = "/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/db.tex";

//	writeTikzTexFileForTreePath(dbfile, path.tree_path, path.opers, NULL, NULL,
//			path.global_index, path.sites, path.path_len);

	for (int i = 0; i < path.path_len - 1; i++)
		if (path.opers[i][0] != -1) {
			if (countNewNodes(path.tree_path[i], path.tree_path[i + 1]) != 1)
				return 0;
		} else if (path.opers[i][0] == -1) {
			if (countNewNodes(path.tree_path[i], path.tree_path[i + 1]) != 0)
				return 0;
		}

	return 1;
}

short checkRecombinationTimes(struct Smc path) {

	struct Tree t_curr, t_next;
	short nl = (path.tree_path[0].n_nodes + 1) / 2, nni, prune_parent;
	double prune_node_time, prune_parent_time, regraft_time;

	for (int i = 0; i < path.path_len - 1; i++) {

		if (path.rec_times[i] <= 0) {
			if (path.opers[i][0] >= 0)
				return 0;
			else
				continue;
		}

		t_curr = path.tree_path[i];
		t_next = path.tree_path[i + 1];

		prune_node_time = path.opers[i][0] < nl ? 0 : t_curr.times[path.opers[i][0] - nl];
		prune_parent = findParent(t_curr.C, nl - 1, nl, path.opers[i][0]);
		prune_parent_time = t_curr.times[prune_parent - nl];
		nni = determineNewNode(t_curr, t_next);
		regraft_time = t_next.times[nni];

		if (path.rec_times[i] < prune_node_time) {
			printf("Step %i, operation %hd %hd, prune node time conflict %.10e %.10e\n", i,
					path.opers[i][0], path.opers[i][1], path.rec_times[i], prune_node_time);
			return 0;
		}
		if (path.rec_times[i] > prune_parent_time) {
			printf("Prune parent time conflict\n");
			return 0;
		}
		if (path.rec_times[i] > regraft_time) {
			printf("Regraft time conflict \n");
			return 0;
		}
	}

	return 1;
}

short countNewNodes(struct Tree t_curr, struct Tree t_next) {

	short nb = (t_curr.n_nodes + 1) / 2 - 1;
	short new_node_count = 0;
	short j;

	for (int i = 0; i < nb; i++) {
		for (j = 0; j < nb; j++)
			if (fabs(t_curr.times[j] - t_next.times[i]) < fabs(t_next.times[i]) * RELTOL)
				break;
		if (j == nb) // not match
			new_node_count++;
	}

	return new_node_count;
}

void testExponentialDisotribution() {
	FILE *file;
	file = fopen("exptest.txt", "w");
	for (long i = 0; i < 100000; i++)
		fprintf(file, "%8.2f\n", exprnd(20000));
	fclose(file);
}

struct Smc synchronisePath(struct SegmentSamplerSynchro sss) {

	struct Smc path;
	short nb = sss.nl - 1, i_last;

	path = createPath(sss.length);
	path.path_len = sss.length;

	for (int i = 0; i < path.path_len - 1; i++) {

		path.tree_path[i] = createTree(sss.nl);
		path.sites[i] = sss.sites[i];

		for (int j = 0; j < nb; j++) {
			path.tree_path[i].C[j] = sss.C[i][j];
			path.tree_path[i].C[j + nb] = sss.C[i][j + nb];
			path.tree_path[i].times[j] = sss.times[i][j];
		}
		path.opers[i] = malloc(sizeof(short) * 2);
		path.opers[i][0] = sss.opers[i][0];
		path.opers[i][1] = sss.opers[i][1];
		path.global_index[i] = malloc(sizeof(short) * nb);
	}
	i_last = path.path_len - 1;

	path.tree_path[i_last] = createTree(sss.nl);
	for (int j = 0; j < nb; j++) {
		path.tree_path[i_last].C[j] = sss.C[i_last][j];
		path.tree_path[i_last].C[j + nb] = sss.C[i_last][j + nb];
		path.tree_path[i_last].times[j] = sss.times[i_last][j];
	}
	path.global_index[i_last] = malloc(sizeof(short) * nb);
	path.sites[i_last] = sss.sites[i_last];

	return path;
}

long findSynchronisedSegment(struct SamplingSet set, struct Smc path) {

	FILE *file;
	char line[1500];
	const char sep[1] = " ";
	struct Smc tmp;
	short nb = (path.tree_path[0].n_nodes + 1) / 2 - 1, set_selector;
	int j, k;

	/* See if the sampled path is from the forced set or not */
	file = fopen("/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/setselectdebug.txt",
			"r");
	fgets(line, sizeof line, file);
	set_selector = (short) atoi(strtok(line, sep));

	/* Change canonical sorting to the time sorting */
	for (int i = 0; i < path.path_len - 1; i++)
		timeSorting(path.tree_path[i], path.opers[i]);
	timeSorting(path.tree_path[path.path_len - 1], NULL);

	/* find match in the generated set */

	for (long i = 0; i < *set.n_paths; i++) {
		if (set.set_selector[i] != set_selector)
			continue;
		tmp = set.paths[i];
		for (k = 0; k < tmp.path_len; k++) {
			for (j = 0; j < nb; j++) {
				if (tmp.tree_path[k].C[j] != path.tree_path[k].C[j]
						|| tmp.tree_path[k].C[j + nb] != path.tree_path[k].C[j + nb])
					break;
			}
			if (j < nb)
				break;
		}
		if (k < tmp.path_len)
			continue;
		return i;
	}
	assert(false);
	return -1;
}

void deallocateSegmentSamplingSynchro(struct SegmentSamplerSynchro sss) {

	for (int i = 0; i < sss.length; i++) {
		free(sss.C[i]);
		free(sss.times[i]);
	}
	for (int i = 0; i < sss.length - 1; i++)
		free(sss.opers[i]);
	free(sss.C);
	free(sss.times);
	free(sss.opers);
	free(sss.sites);

}

struct FreeTimeSynchro synchroniseFreeTimes() {

	FILE *file;
	char line[1500];
	const char sep[1] = " ";
	char *token;
	struct FreeTimeSynchro out;
	short cnt = 0;

	file = fopen("/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/freetimesdebug.txt",
			"r");

	out.times = malloc(sizeof(double) * 50);

	fgets(line, sizeof line, file);
	token = strtok(line, sep);
	while (token != NULL) {
		out.times[cnt++] = (double) atof(token);
		token = strtok(NULL, sep);
	}
	out.length = cnt;

	return out;
}

struct InitTimeSynchro synchroniseInitTimeProposal() {

	FILE *file;
	char line[1500];
	const char sep[1] = " ";
	char *token;
	struct InitTimeSynchro out;
	short cnt = 0;

	file = fopen("/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/inittimesdebug.txt",
			"r");

	out.indices = malloc(sizeof(int) * 50);
	out.times = malloc(sizeof(double) * 50);

	while (fgets(line, sizeof line, file) != NULL) {
		token = strtok(line, sep);
		out.times[cnt] = (double) atof(token);
		if (out.times[cnt] == 0.0)
			break;
		out.indices[cnt++] = (int) atoi(token);
		token = strtok(NULL, sep);
	}

	out.length = cnt;

	return out;
}

struct RecTimeSynchro synchroniseRecombinationTimes() {

	FILE *file;
	char line[1500];
	const char sep[1] = " ";
	char *token;
	struct RecTimeSynchro out;
	short cnt = 0;

	file = fopen("/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/rectimesdebug.txt",
			"r");

	out.times = malloc(sizeof(double) * 50);

	fgets(line, sizeof line, file);
	token = strtok(line, sep);
	while (token != NULL) {
		out.times[cnt] = (double) atof(token);
		if (out.times[cnt] == 0.0)
			break;
		cnt++;
		token = strtok(NULL, sep);
	}
	out.length = cnt;

	return out;
}

struct AcceptanceSynchro synchroniseAcceptance() {
	FILE *file;
	char line[1500];
	const char sep[1] = " ";
	char *token;
	struct AcceptanceSynchro out;
	file = fopen("/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/acceptdebug.txt",
			"r");
	fgets(line, sizeof line, file);
	token = strtok(line, sep);
	out.accept = (short) atoi(token);
	return out;
}
