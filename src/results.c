/*
 * results.c
 *
 *  Created on: 29.11.2016
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
#include "datastructures.h"
#include "debugging.h"
#include "fileaccess.h"
#include "graph2tikz.h"
#include "MCMCutils.h"
#include "pathutils.h"
#include "timeadjustment.h"
#include "utils.h"

void map(struct MCMCSummary *chain, long chain_length, long full_scan_count,
		struct Data data) {

	struct Smc path, map_path;
	struct CanonicalMatrix cmtx;
	struct UniqueRowsWithCounts unique_rows;
	struct ShortRowSortedMtx sorted_cmtx;
	struct Smc *paths;
	short nb = (chain[0].path.tree_path[0].n_nodes + 1) / 2 - 1, max_len = 0, len,
			count = 0, n_global = 0;
	short *cmtx_tmp;
	long max_count = 0, max_path = -1, avg_N = 0;
	double *avg_times, *avg_rec_times;
	enum Node_Colors **colors;

	// find maximum path length
	for (int i = 0; i < chain_length; i++) {
		len = chain[i].path.path_len;
		max_len = (chain[i].full_scan == 1 && len > max_len) ? len : max_len;
	}

	/* Construct a matrix whose each row is a canonical representation of a path */
	paths = malloc(sizeof(struct Smc) * full_scan_count);

	cmtx.n = max_len * (2 * nb + 2);
	cmtx.m = full_scan_count;
	cmtx.data = malloc(sizeof(short) * cmtx.m * cmtx.n);

	for (long i = 0; i < chain_length; i++) {

		if (chain[i].full_scan != 1)
			continue;

		path = chain[i].path;
		paths[count] = path;

		cmtx_tmp = canonisePath(path);

		for (int j = 0; j < path.path_len * (2 * nb + 2); j++)
			cmtx.data[count + j * cmtx.m] = cmtx_tmp[j];
		// pad the end with zeros
		for (int j = path.path_len * (2 * nb + 2); j < max_len * (2 * nb + 2); j++)
			cmtx.data[count + j * cmtx.m] = 0;

		count++;
		free(cmtx_tmp);
	}

	/* Sort and unique the paths */
	sorted_cmtx = sortCanonicalMtx(cmtx);
	unique_rows = uniqueRowsWithCounts(sorted_cmtx);

	// find the path with maximum number of occurrences
	for (int i = 0; i < unique_rows.m; i++)
		if (unique_rows.counts[i] > max_count) {
			max_count = unique_rows.counts[i];
			max_path = i;
		}

	// average coalescence times for the MAP ARG
	avg_times = malloc(sizeof(double) * nb * max_len);
	fillDoubleArray(avg_times, 0.0, 1, nb * max_len, 1);
	avg_rec_times = malloc(sizeof(double) * (max_len - 1));
	fillDoubleArray(avg_rec_times, 0.0, 1, max_len - 1, 1);

	map_path.tree_path = NULL;

	for (int i = 0; i < sorted_cmtx.m; i++)
		// iterate over the whole path sequence

		if (unique_rows.inds[i] == max_path) {
			// if current paths coincides with the path with maximum count

			// create the map path if not created already
			if (map_path.tree_path == NULL)
				map_path = createPathCopy(paths[sorted_cmtx.i_sort[i]]);

			// average the coalescence times
			for (int j = 0; j < paths[sorted_cmtx.i_sort[i]].path_len; j++)
				for (int k = 0; k < nb; k++) {
					avg_times[k + j * nb] += paths[sorted_cmtx.i_sort[i]].tree_path[j].times[k];
				}

			for (int j = 0; j < paths[sorted_cmtx.i_sort[i]].path_len - 1; j++) {
				avg_rec_times[j] += paths[sorted_cmtx.i_sort[i]].rec_times[j];
			}
			avg_N++;
		}

	// TODO Clean outputs
	// coalescence times

	// set the times to the averages
	for (int i = 0; i < map_path.path_len; i++)
		for (int j = 0; j < nb; j++)
			map_path.tree_path[i].times[j] = avg_times[j + i * nb] / (double) avg_N;

	for (int i = 0; i < map_path.path_len - 1; i++)
		map_path.rec_times[i] = avg_rec_times[i] / (double) avg_N;

	colors = colorPath(map_path, data);

	for (int i = 0; i < map_path.path_len; i++)
		free(map_path.global_index[i]);
	free(map_path.global_index);

	map_path.global_index = malloc(sizeof(short*) * map_path.path_len);
	for (int i = 0; i < map_path.path_len; i++)
		assignGlobalIndices(map_path.global_index, &n_global, i, map_path.tree_path);

	writeTikzFileForMapPath(map_path);

	assert(checkTreePathCompletely(map_path) == 1);
	assert(checkCompatibility(map_path, data) == 1);
	assert(checkRecombinationTimes(map_path) == 1);

	for (int i = 0; i < map_path.path_len; i++)
		free(colors[i]);
	free(colors);

	deallocatePath(map_path);

	free(avg_times);
	free(avg_rec_times);
	free(sorted_cmtx.data);
	free(sorted_cmtx.i_sort);
	free(unique_rows.counts);
	free(unique_rows.data);
	free(unique_rows.inds);
	free(unique_rows.selector);
	free(paths);
}
