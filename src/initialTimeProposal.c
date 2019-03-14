/*
 * initialTimeProposal.c
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

#include <float.h>
#include <assert.h>
#include "constants.h"
#include "datastructures.h"
#include "debugging.h"
#include "graph2tikz.h"
#include "treeutils.h"
#include "utils.h"
#include "randomness.h"
#include "timeadjustment.h"

double initialTimeProposal(struct Smc path, struct Parameters parm, short draw) {

	struct Tree t_curr, t_last;
	short nl = (path.tree_path[0].n_nodes + 1) / 2;
	short nb = nl - 1, n_free_times = 0;
	short **affixed;
	double old_time, upper, lower, t_cmp, mu, c, p = 1, new_time;
	int jj, kk;

	t_last = path.tree_path[path.path_len - 1];

	/* flag affixed nodes */
	affixed = malloc(sizeof(short*) * path.path_len);
	for (int i = 0; i < path.path_len; i++) {

		t_curr = path.tree_path[i];
		affixed[i] = malloc(sizeof(short) * nb);
		for (int j = 0; j < nb; j++)
			affixed[i][j] = findMatchingNode(t_last, t_curr.times[j]) >= 0 ? 1 : 0;

	}

	//	printf("Affix indicators\n");
	//	for (int j = 0; j < nb; j++) {
	//		for (int i = 0; i < path.path_len; i++)
	//			printf("%2hd", affixed[i][j]);
	//		printf("\n");
	//	}

	for (int i = 0; i < path.path_len; i++) {
		for (int j = 0; j < nb; j++) {

			if (affixed[i][j] == 1)
				continue;

			old_time = path.tree_path[i].times[j];
//			printf("Old time %1.5f\n", old_time);

			/*
			 * Upper limits
			 */
			upper = DBL_MAX;
			jj = j;
			while (jj < nb && affixed[i][jj] == 0) {

				// iterate over all following trees
				for (int ii = i; ii < path.path_len; ii++) {

					// see if the free upper bounding time is found the the following tree
					for (kk = 0; kk < nb; kk++)
						if (path.global_index[ii][kk] == path.global_index[i][jj])
							break;
					if (kk < nb) { // found
						// find the first greater affixed time
						for (int k = kk; k < nb; k++) {
							t_cmp = path.tree_path[ii].times[k];
							if (affixed[ii][k] == 1) {
								upper = t_cmp < upper ? t_cmp : upper;
								break;
							}
						}
					}
				}
				jj++;
			}
//			printf("upper=%1.5f\n", upper);

			/*
			 * Lower limits
			 */
			lower = 0;
			jj = j;
			while (jj >= 0 && affixed[i][jj] == 0) {

				// iterate over all following trees
				for (int ii = i; ii < path.path_len; ii++) {

					// see if the free lower bounding time is found the the following tree
					for (kk = 0; kk < nb; kk++)
						if (path.global_index[ii][kk] == path.global_index[i][jj])
							break;
					if (kk < nb) { // found
						// find the first smaller affixed time
						for (int k = kk; k >= 0; k--) {
							t_cmp = path.tree_path[ii].times[k];
							if (affixed[ii][k] == 1) {
								lower = t_cmp > lower ? t_cmp : lower;
								break;
							}
						}
					}
				}
				jj--;
			}
//			printf("lower=%1.5f\n", lower);

			assert(upper > lower);

			mu = (double) (2 * parm.n_eff) / (double) nchoosek(nl - j, 2);
			if (draw == 1) {
				c = exprnd(mu);
				while (c >= upper - lower)
					c = exprnd(mu);
			} else {
				c = path.tree_path[i].times[j] - lower;
			}
			p *= expPdf(c, mu) / expCdf(upper - lower, mu);
			n_free_times++;
			new_time = lower + c + 2 * lower * RELTOL;

			/* set new time throughout the path */
			path.tree_path[i].times[j] = new_time;
			affixed[i][j] = 1;
			for (int ii = i + 1; ii < path.path_len; ii++) {
				for (jj = 0; jj < nb; jj++)
					if (path.global_index[ii][jj] == path.global_index[i][j]) {
						path.tree_path[ii].times[jj] = new_time;
						affixed[ii][jj] = 1;
						break;
					}
				if (jj == nb)
					// if adjustable node was not found, it won't be found in the
					// remaining trees either
					break;
			}
		}
	}

	for (int i = 0; i < path.path_len; i++)
		free(affixed[i]);
	free(affixed);

	for (int i = 0; i < path.path_len; i++)
		assert(checkTree(path.tree_path[i]) == 1);

	return p;
}

