/*
 * smcPrior.c
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

#include <math.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include "constants.h"
#include "datastructures.h"
#include "graph2tikz.h"
#include "pathutils.h"
#include "smcPrior.h"
#include "treeutils.h"
#include "utils.h"

struct SmcPriorData smcprior(struct Smc path, struct Parameters parm, struct Data data) {

	struct SmcPriorData out;
	struct BranchLength br;
	struct Tree t_curr, t_next;
	struct RegraftTimeData rg_data;
	double pr_branch_len, pr_node_time, pr_prnt_time, rg_time, p_rec, p_rec_time,
			p_rg_branch, p = 0, div, lam;
	struct NonIdentityOps operations;
	short *rs_indicator;
	short nl, nb, was_recombination = 1, rec_count = 0, pr_node, new_node_idx;
	double *components;
	double *incrs;

	incrs = malloc(sizeof(double) * path.selector_length);
	components = malloc(sizeof(double) * path.selector_length);
	nl = (path.tree_path[0].n_nodes + 1) / 2;
	nb = nl - 1;
	br.lenghts = NULL;
	div = (double) (2 * parm.n_eff);
	/* Initial density:
	 * Because in the coalescent model the coalescing pairs are chosen in a
	 * uniformly random manner, every tree topology is equally probable.
	 * Therefore there is in fact no need to evaluate the probability of the tree
	 * topology as it will cancel out in the acceptance probability
	 * calculation.*/
	t_curr = path.tree_path[0];
	for (int i = 1; i < nb; i++) {
		lam = (double) nchoosek(i + 1, 2) / div;
		p += log(lam) - lam * (t_curr.times[nb - i] - t_curr.times[nb - i - 1]);
	}
	lam = (double) nchoosek(nb + 1, 2) / div;
	p += log(lam) - lam * t_curr.times[0];

	components[0] = p;
	rs_indicator = recombinationSiteIndicator(path);
	operations = getNonIdentityOperations(path);

	/* the remaining tree sequence (iterate over the full length sequence)*/
	for (int i = 0; i < path.selector_length - 1; i++) {
		// tree has changed, so we need to recalculate the branch lengths
		if (was_recombination == 1) {
			if (br.lenghts != NULL)
				free(br.lenghts);
			br = branchLengths(path.tree_path[path.tree_selector[i]], data, -1);
		}
		was_recombination = rs_indicator[i];

		if (rs_indicator[i] == 1) {

			t_curr = path.tree_path[path.tree_selector[i]];
			t_next = path.tree_path[path.tree_selector[i + 1]];

			pr_node = operations.ops[rec_count][0];
			pr_branch_len = br.lenghts[pr_node];
			pr_node_time = pr_node < nl ? 0 : t_curr.times[pr_node - nl];
			pr_prnt_time = pr_node_time + pr_branch_len;

			new_node_idx = determineNewNode(t_curr, t_next);
			rg_time = t_next.times[new_node_idx];

			/* joint probability of "recombination occurs" AND "it happens on 'prune_node'
			 * branch " */
			p_rec = pr_branch_len / br.total_length
					* ((double) 1 - exp(-parm.rho * br.total_length));

			assert(
					fabs(pr_prnt_time - pr_node_time - pr_branch_len)
							< (pr_prnt_time - pr_node_time + pr_branch_len) * TOL);
			assert(fmin(rg_time, pr_prnt_time) - pr_node_time <= pr_branch_len);

			/* uniform density over the prune branch */
			p_rec_time = (double) 1 / (fmin(rg_time, pr_prnt_time) - pr_node_time);

			/* regraft time density */
			rg_data = rgTimePdf(rg_time, t_curr, operations.rec_times[rec_count], pr_prnt_time,
					br.lenghts, parm);

			p_rg_branch = (double) 1 / (double) rg_data.n_active;

			incrs[i + 1] = log(p_rec) + log(p_rec_time) + log(p_rg_branch)
					+ log(rg_data.p_rg_time);
			p += incrs[i + 1];

			assert(isnan(incrs[i+1]) == 0);
			assert(isinf(incrs[i+1]) == 0);
			rec_count++;

		} else {
			incrs[i + 1] = -parm.rho * br.total_length;
			p += incrs[i + 1];
		}
		components[i + 1] = p;
	}

	if (br.lenghts != NULL)
		free(br.lenghts);
	free(rs_indicator);

	deallocateOperations(operations);

	out.density = p;
	out.components = components;
	out.increments = incrs;
	out.length = path.selector_length;
	out.number_of_recombinations = rec_count;
	assert(isnan(p)==0);

	return out;
}

void deallocatePriorData(struct SmcPriorData priordata) {
	free(priordata.components);
	free(priordata.increments);
}

void testrgpdf(struct Tree t, double rec_time, double pr_prnt_time, double *lengths,
		struct Parameters parm) {

	FILE *file;
	int n = 10000;
	short nb = (t.n_nodes + 1) / 2 - 1;
	double x;
	struct RegraftTimeData rg_data;

	file = fopen("rgtest.txt", "w");

	for (int i = 0; i <= n; i++) {
		x = i * (t.times[nb - 1] / (double) n);
		rg_data = rgTimePdf(x, t, rec_time, pr_prnt_time, lengths, parm);
		fprintf(file, "%15.5e %15.5e %hd\n", x, rg_data.p_rg_time, rg_data.n_active);
	}
	fclose(file);
}

struct RegraftTimeData rgTimePdf(double x, struct Tree t, double rec_time,
		double pr_prnt_time, double *br_lengths, struct Parameters parm) {

	struct RegraftTimeData out;
	double *dt, *t_low, *t_upp, *t_low_trunc;
	double t_mid, t_tmp, y, div, tmp, normaliser;
	short *active_counts;
	short nl = (t.n_nodes + 1) / 2, hit;

	t_low = malloc(sizeof(double) * nl);
	t_upp = malloc(sizeof(double) * nl);
	t_low[0] = (double) 0;
	t_upp[nl - 1] = DBL_MAX;

	for (int i = 0; i < nl - 1; i++) {
		t_upp[i] = t.times[i];
		t_low[i + 1] = t.times[i];
	}

	out.n_active = activeBranchCount(x, nl, br_lengths, t.times, pr_prnt_time);
	active_counts = malloc(sizeof(short) * nl);
	for (int i = 0; i < nl; i++) {
		t_mid = (t_low[i] + t_upp[i]) / (double) 2;
		active_counts[i] = activeBranchCount(t_mid, nl, br_lengths, t.times, pr_prnt_time);
	}

	dt = malloc(sizeof(double) * nl);
	t_low_trunc = malloc(sizeof(double) * nl);
	for (int i = 0; i < nl; i++) {
		t_low_trunc[i] = t_low[i] > rec_time ? t_low[i] : rec_time;
		t_tmp = t_upp[i] - t_low_trunc[i];
		dt[i] = t_tmp < 0 ? 0 : t_tmp;
	}

	/* find the time interval x hits */
	for (hit = 0; hit < nl; hit++)
		if (x >= t_low_trunc[hit] - t_low_trunc[hit] * TOL && x < t_upp[hit])
			break;

// evaluate the density
	if (hit < nl) {
		y = evalDensity(x, (double) out.n_active, t_low_trunc, hit, (double) parm.n_eff,
				active_counts, dt);
	} else {
		y = 0;
	}

// normalise the density
	div = (double) (2 * parm.n_eff);
	normaliser = 0;
	for (int i = 0; i < nl; i++) {
		tmp = 0;
		for (int j = 0; j < i; j++)
			tmp += (double) active_counts[j] * dt[j];
		normaliser += ((double) 1 - exp(-(double) active_counts[i] * dt[i] / div))
				* exp(-tmp / div);
//		printf("Normaliser cum %.10f %.10f  %.10f\n",normaliser,tmp,active_counts[i] * dt[i]);
	}
	out.p_rg_time = y / normaliser;
//	printf("%.10f\n",normaliser);
//	printShortArray(active_counts,1,nl,1);
//	printDoubleArray(dt,nl,1,nl);

	free(t_low);
	free(t_upp);
	free(active_counts);
	free(dt);
	free(t_low_trunc);
	return out;
}

struct NonIdentityOps getNonIdentityOperations(struct Smc path) {

	struct NonIdentityOps out;
	short** ops;
	double *rec_times;
	short nops = 0;

	ops = malloc(sizeof(short*) * (path.path_len - 1));
	rec_times = malloc(sizeof(double) * (path.path_len - 1));

	for (int i = 0; i < path.path_len - 1; i++)
		if (path.opers[i][0] != -1) {
			ops[nops] = createOperationCopy(path.opers[i]);
			rec_times[nops++] = path.rec_times[i];
		}

	out.ops = realloc(ops, sizeof(short*) * nops);
	out.rec_times = realloc(rec_times, sizeof(double) * nops);
	out.nops = nops;
	return out;
}

void deallocateOperations(struct NonIdentityOps ops) {
	for (int i = 0; i < ops.nops; i++)
		free(ops.ops[i]);
	free(ops.ops);
	free(ops.rec_times);
}

double evalDensity(double x, double n_act, double *t_lo_trnc, short hit, double n_eff,
		short *active_counts, double *dt) {

	double tmp_sum = 0, out, div = (double) (2 * n_eff);

	for (int i = 0; i < hit; i++)
		tmp_sum += dt[i] * (double) active_counts[i];

	out = n_act * exp(-n_act * (x - t_lo_trnc[hit]) / div - tmp_sum / div) / div;
	return out;
}

short activeBranchCount(double x, short nl, double *branch_lengths, double *times,
		double pr_prnt_time) {

	short out = 0;

	/* x is beyond the root */
	if (x > times[nl - 2])
		return 1;

	for (int i = 0; i < 2 * nl - 2; i++)
		if (i < nl)
			out += x < branch_lengths[i] ? 1 : 0;
		else
			out +=
					x >= times[i - nl] - times[i - nl] * TOL
							&& x < times[i - nl] + branch_lengths[i] ? 1 : 0;

	if (x <= pr_prnt_time + pr_prnt_time * TOL)
		out--;

	return out;
}
