/*
 * freeTimes.c
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

#include <float.h>
#include <math.h>
#include <assert.h>
#include "treeutils.h"
#include "datastructures.h"
#include "debugging.h"
#include "freeTimes.h"
#include "utils.h"
#include "constants.h"
#include "randomness.h"
#include "pathutils.h"
#include "graph2tikz.h"

double freetimes(struct Smc path, struct Parameters parm, struct Conditioning cond,
		short draw) {

	struct Smc path_in;
	struct Tree t_prev, t_curr;
	struct DoubleVector upper, lower;
	short nn_idx, nl, nb, nn_gi;
	double p = 1, free_time, t_min, t_max, t_sampled = -1, increment;
	int k;
	short *pinned;

	nl = (path.tree_path[0].n_nodes + 1) / 2;
	nb = nl - 1;
	assert(path.is_free!=NULL);
	path_in = createPathCopy(path);

	assert(path.is_free[0] == 0);
	// because the first tree never contains a free time we can start the iteration for
	// the second tree
	for (int i = 1; i < path.path_len; i++) {
		/* jump to next iteration if not a free time */
		if (path.is_free[i] == 0)
			continue;

		t_prev = path_in.tree_path[i - 1];
		t_curr = path_in.tree_path[i];
		nn_idx = determineNewNode(t_prev, t_curr);
		if (!(nn_idx >= 0)) {
			for (int j = 0; j < path_in.path_len; j++)
				printTree(path_in.tree_path[j]);
			assert(false);
		}
		free_time = path_in.tree_path[i].times[nn_idx];
		nn_gi = path_in.global_index[i][nn_idx];

		upper.v = malloc(sizeof(double) * (path.path_len - i));
		lower.v = malloc(sizeof(double) * (path.path_len - i));
		upper.length = path.path_len - i;
		lower.length = path.path_len - i;
		fillDoubleArray(upper.v, DBL_MAX, 1, path.path_len - i, 1);
		fillDoubleArray(lower.v, -DBL_MAX, 1, path.path_len - i, 1);

		/* the bounds in the tree where the free node is introduced are always taken from the
		 * modified tree */
		lower.v[0] = nn_idx == 0 ? 0 : path.tree_path[i].times[nn_idx - 1];
		upper.v[0] = nn_idx == nb - 1 ? DBL_MAX : path.tree_path[i].times[nn_idx + 1];

		if (cond.twosided == 1) {
			// iterate over the later trees
			for (int j = i + 1; j < path.path_len; j++) {
				pinned = pinnedTimes(path.global_index[j], path.global_index[path.path_len - 1],
						nb);
				// upper bound
				for (int kk = 0; kk < nb; kk++)
					/* if node is pinned and in the unmodified tree has a time greater
					 * than the original free time, then it is an upper bounder */
					if (pinned[kk] == 1 && path_in.tree_path[j].times[kk] > free_time) {
						/* bounding time itself is taken from the modified tree */
						upper.v[j - i] = path.tree_path[j].times[kk];
						break;
					}
				// lower bound
				for (int kk = nb - 1; kk >= 0; kk--)
					if (pinned[kk] == 1 && path_in.tree_path[j].times[kk] < free_time) {
						lower.v[j - i] = path.tree_path[j].times[kk];
						break;
					}

				free(pinned);
			}
		}

		t_min = findMax(lower);
		t_max = findMin(upper);

		if (t_max < DBL_MAX) {
			if (draw == 1) {
				t_sampled = t_min + 2 * t_min * RELTOL
						+ (t_max - 2 * t_max * RELTOL - t_min - 2 * t_min * RELTOL) * genrand_real3();
			}
			p *= (double) 1 / (t_max - t_min);
		} else {
			if (draw == 1) {
				increment = exprnd((double) (2 * parm.n_eff));
				t_sampled = t_min + increment + 2 * t_min * RELTOL;
			} else {
				increment = free_time - t_min;
			}
			p *= expPdf(increment, (double) (2 * parm.n_eff));
		}

		if (draw == 1) {
			assert(t_sampled >= 0);
			for (int j = 0; j < path.path_len; j++)
				for (k = 0; k < nb; k++)
					if (path.global_index[j][k] == nn_gi) {
						path.tree_path[j].times[k] = t_sampled;
						break;
					}
		}

		free(upper.v);
		free(lower.v);

	}

	assert(checkTreePathCompletely(path)==1);
	deallocatePath(path_in);

	return p;
}

short *pinnedTimes(short *gi, short *gir, short nb) {

	short *out;
	out = malloc(sizeof(short) * nb);
	fillShortArray(out, 0, 1, nb, 1);

	for (int i = 0; i < nb; i++)
		for (int j = 0; j < nb; j++)
			if (gi[i] == gir[j]) {
				out[i] = 1;
				break;
			}

	return out;
}

