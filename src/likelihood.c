/*
 * likelihood.c
 *
 *  Created on: 13.11.2016
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
#include <stdio.h>
#include "datastructures.h"
#include "likelihood.h"
#include "data.h"
#include "treeutils.h"

struct LikelihoodData likelihood(struct Smc path, struct Data data,
		struct Parameters parm) {

	struct LikelihoodData out;
	struct Tree t_curr;
	short is_mutation, data_col;
	struct BranchLength branch;
	double incr;

	out.increments = malloc(sizeof(double) * path.selector_length);
	out.total_branch_length = malloc(sizeof(double) * path.selector_length);
	out.mutated_branch_lengths = malloc(sizeof(double) * path.selector_length);
	out.length = path.selector_length;
	out.log_likelihood = 0;

	for (long i = 0; i < path.selector_length; i++) {
		t_curr = path.tree_path[path.tree_selector[i]];

		data_col = getDataColumn(data, i + 1);
		is_mutation = isMutationSite(data_col, data);

		branch = branchLengths(t_curr, data, data_col);

		if (is_mutation == 1)
			incr = log((double) 1 - exp(-parm.mu * branch.total_length))
					+ log(branch.mutated_length / branch.total_length);
		else
			incr = -parm.mu * branch.total_length;

		out.increments[i] = incr;
		out.total_branch_length[i] = branch.total_length;
		out.mutated_branch_lengths[i] = branch.mutated_length;
		out.log_likelihood += incr;
		free(branch.lenghts);
	}

	return out;
}

/* Returns 1 if column 'col' in 'data' has any non zero entries. 0 is returned if
 * col < -1 or the column contains only zeros. */
short isMutationSite(long col, struct Data data) {

	if (col == -1)
		return 0;

	for (int i = 0; i < data.n_seq; i++)
		if (data.M[i + col * data.n_seq] == 1)
			return 1;

	return 0;
}

void deallocateLikelihood(struct LikelihoodData like) {
	free(like.increments);
	free(like.mutated_branch_lengths);
	free(like.total_branch_length);
}
