/*
 * recombinationTimes.c
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

#include <assert.h>
#include <math.h>
#include "constants.h"
#include "datastructures.h"
#include "treeutils.h"
#include "randomness.h"

struct Recombinations recombinationTimes(struct Smc path, short draw) {

	struct Recombinations rec;
	short pr_node, rg_node, nl, new_node_idx, pr_prnt, nb;
	double rg_time, start_time, end_time, pr_prnt_time;

	rec.log_p = 0;
	rec.n_recombinations = 0;

	nl = (path.tree_path[0].n_nodes + 1) / 2;
	nb = nl - 1;

	for (int i = 0; i < path.path_len - 1; i++) {

		pr_node = path.opers[i][0];
		rg_node = path.opers[i][1];

		if (pr_node == -1 && rg_node == -1)
			continue;

		start_time = pr_node < nl ? 0 : path.tree_path[i].times[pr_node - nl];
		pr_prnt = findParent(path.tree_path[i].C, nb, nl, pr_node);
		assert(pr_prnt >= 0);
		pr_prnt_time = path.tree_path[i].times[pr_prnt - nl];
		new_node_idx = determineNewNode(path.tree_path[i], path.tree_path[i + 1]);
		rg_time = path.tree_path[i + 1].times[new_node_idx];

		end_time = rg_time < pr_prnt_time ? rg_time : pr_prnt_time;

		if (draw == 1)
			path.rec_times[i] = start_time + 2 * start_time * RELTOL
					+ (end_time - 2 * end_time * RELTOL - start_time - 2 * start_time * RELTOL)
							* genrand_real3();

		rec.log_p += log((double) 1 / (end_time - start_time));
		rec.n_recombinations++;
	}

	return rec;
}
