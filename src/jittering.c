/*
 * jittering.c
 *
 *  Created on: 15.11.2016
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

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "constants.h"
#include "datastructures.h"
#include "debugging.h"
#include "freeTimes.h"
#include "fileaccess.h"
#include "graph2tikz.h"
#include "jittering.h"
#include "likelihood.h"
#include "MCMCutils.h"
#include "pathutils.h"
#include "randomness.h"
#include "recombinationTimes.h"
#include "smcPrior.h"
#include "sorting.h"
#include "timeadjustment.h"
#include "treeutils.h"
#include "utils.h"

long jittering(struct MCMCSummary *chain, long it, long max_it, struct Parameters parm,
		struct Data data) {

	struct TimeParms gtimes;
	struct Smc path;
	struct MCMCDiagnostics dgn;
	struct LikelihoodData like;
	struct SmcPriorData prior;
	short nl, n_times = 0, nb;
	double t_min, t_max, t_curr, incr, old_time, lam, q_prop, q_curr;

	nl = (chain[it - 1].path.tree_path[0].n_nodes + 1) / 2;
	nb = nl - 1;

	gtimes = collectTimes(chain[it - 1].path, &n_times, parm);

	if (parm.verb == 1)
		printf("*****************************************\n");

	for (int i = 0; i < n_times; i++) {

		if (parm.verb == 1)
			printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bJitter: %5i/%-5i", i + 1, n_times);

		path = createPathCopy(chain[it - 1].path);

		t_min = -DBL_MAX;
		t_max = DBL_MAX;
		for (int j = 0; j < path.path_len; j++)
			for (int k = 0; k < nb; k++) {
				t_curr = path.tree_path[j].times[k];
				t_min = t_curr < gtimes.times[i] && t_curr > t_min ? t_curr : t_min;
				t_max = t_curr > gtimes.times[i] && t_curr < t_max ? t_curr : t_max;
			}
		t_min = t_min < 0 ? 0 : t_min;

		incr = -1;
		q_prop = -1;

		/* rejection sampling */
		if (t_max == DBL_MAX)
			lam = gtimes.parms[i] / (double) 4;
		else
			lam = gtimes.parms[i];

		while (!(incr > RELTOL * t_min && incr + t_min < t_max))
			incr = exprnd(lam);
		incr += t_min;
		old_time = gtimes.times[i];

		adjustTimeGlobally(path, incr, old_time);
		q_prop = expPdf(incr - t_min, lam) / expCdf(t_max - t_min, lam);
		assert(q_prop >= 0);

		/* proposal density for the current time */
		q_curr = expPdf(old_time - t_min, lam) / expCdf(t_max - t_min, lam);

		dgn = calculateJitterAlpha(path, chain[it - 1].path, data, parm, q_curr / q_prop);
		dgn.jitter_step = 1;
		dgn.proposed_number_of_recombinations = path.path_len - 1;
		dgn.current_number_of_recombinations = path.path_len - 1;
		dgn.u = genrand_real3();
		dgn.accept_indicator = dgn.u < dgn.alpha ? 1 : 0;

		if (dgn.accept_indicator == 1) {
			// ACCPT: Insert the jittered copy of the previous state
			chain[it].full_scan = 0;
			chain[it++].path = path;
		} else {
			// REJECT: Create copy of the current path
			chain[it].full_scan = 0;
			chain[it].path = createPathCopy(chain[it - 1].path);
			it++;
			deallocatePath(path);
		}

		like = likelihood(chain[it - 1].path, data, parm);
		dgn.log_likelihood = like.log_likelihood;
		prior = smcprior(chain[it - 1].path, parm, data);
		dgn.log_prior = prior.density;
		deallocatePriorData(prior);
		dgn.log_posterior = dgn.log_likelihood + dgn.log_prior;
		deallocateLikelihood(like);

		if (dgn.accept_indicator == 0)
			if (fabs(dgn.log_likelihood - chain[it - 2].data.log_likelihood)
					> (fabs(chain[it - 2].data.log_likelihood) + fabs(dgn.log_likelihood)) * TOL)
				// reject sample, but likelihood changes
				assert(false);

		writeDiagnosticsFile(dgn);
		chain[it - 1].data = dgn;

//		if (parm.verb == 0) {
//			printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%10ld/%-10ld", it, max_it);
//			fflush(stdout);
//		}

//		appendChainTikzFile(chain[it - 1].path, data, 0);

		if (it >= max_it) {
			free(gtimes.parms);
			free(gtimes.times);
			return it;
		}
	}

	free(gtimes.times);
	free(gtimes.parms);

	if (parm.verb == 1) {
		printf("Done\n");
		printf("*****************************************\n");
	}
	return it;
}

void adjustTimeGlobally(struct Smc path, double new_time, double old_time) {

	short nl = (path.tree_path[0].n_nodes + 1) / 2;
	short nb = nl - 1;
	int j = -1, i;

	for (i = 0; i < path.path_len; i++) {
		for (j = 0; j < nb; j++)
			if (fabs(path.tree_path[i].times[j] - old_time)
					< (old_time + path.tree_path[i].times[j]) * TOL)
				break;
		if (j < nb)
			break;
	}
	assert(j >= 0);
	forwardTimeAdjustment(path, i, new_time, j + nl, path.global_index[i][j]);
}

struct TimeParms collectTimes(struct Smc path, short *n_times, struct Parameters parm) {

	struct TimeParms out;
	double *all_times;
	double *exp_parms;
	short nl = (path.tree_path[0].n_nodes + 1) / 2;
	short nb = nl - 1;
	struct Tree tmp = createTree(nl);
	int k, n;

	/* collect all node times in the path */
	all_times = malloc(sizeof(double) * path.path_len * nb);
	exp_parms = malloc(sizeof(double) * path.path_len * nb);
	*n_times = 0;
	for (int i = 0; i < path.path_len; i++) {
		copyTree(&tmp, path.tree_path[i]);
		canonise(tmp);
		for (int j = 0; j < nb; j++) {
			for (k = 0; k < n_times[0]; k++)
				if (fabs(all_times[k] - tmp.times[j]) < (all_times[k] + tmp.times[j]) * RELTOL)
					break;
			if (k == n_times[0]) { // add
				all_times[n_times[0]] = tmp.times[j];
				for (n = 0; n < nb; n++)
					if (fabs(path.tree_path[i].times[n] - tmp.times[j])
							< (path.tree_path[i].times[n] + tmp.times[j]) * RELTOL)
						break;
				exp_parms[n_times[0]++] = (double) (2 * parm.n_eff)
						/ (double) nchoosek(nl - n, 2);
			}
		}
	}

	out.times = all_times;
	out.parms = exp_parms;
	out.length = *n_times;

	deleteTree(tmp);

	return out;
}

struct MCMCDiagnostics calculateJitterAlpha(struct Smc path_p, struct Smc path_c,
		struct Data data, struct Parameters parm, double extra) {

	struct MCMCDiagnostics dgn;
	struct LikelihoodData like_p, like_c;
	struct SmcPriorData prior_p, prior_c;
	struct Recombinations rec_p, rec_c;
	double like, prior, rec_time;

	/* recombination times need to be resampled because they may have become
	 * incompatible */
	rec_p = recombinationTimes(path_p, 1);
	assert(checkRecombinationTimes(path_p)==1);
//	js = synchroniseJitter();
//	for (int i = 0; i < js.length; i++)
//		path_p.rec_times[i] = js.rec_times[i];
//	rec_p = recombinationTimes(path_p, 0);
//	free(js.rec_times);
	rec_c = recombinationTimes(path_c, 0);
	dgn.proposed_recombination_density = rec_p.log_p;
	dgn.current_recombination_density = rec_c.log_p;
	rec_time = exp(rec_c.log_p - rec_p.log_p);

	/* free times do not apply */
	dgn.proposed_free_time_density = -1;
	dgn.proposed_free_time_density = -1;

	/* likelihood */
	like_p = likelihood(path_p, data, parm);
	like_c = likelihood(path_c, data, parm);
	like = exp(like_p.log_likelihood - like_c.log_likelihood);
	dgn.proposed_log_likelihood = like_p.log_likelihood;
	dgn.current_log_likelihood = like_c.log_likelihood;
	deallocateLikelihood(like_p);
	deallocateLikelihood(like_c);

	/* prior */
	prior_p = smcprior(path_p, parm, data);
	prior_c = smcprior(path_c, parm, data);
	prior = exp(prior_p.density - prior_c.density);
	dgn.proposed_log_prior = prior_p.density;
	dgn.current_log_prior = prior_c.density;
	deallocatePriorData(prior_p);
	deallocatePriorData(prior_c);

	dgn.alpha = like * prior * rec_time * extra;
	dgn.alpha = dgn.alpha > 1 ? 1 : dgn.alpha;
	dgn.indicators[0] = 'N';
	dgn.indicators[1] = 'A';
	dgn.indicators[2] = '\0';

	return dgn;
}

struct JitterSychro synchroniseJitter() {

	FILE *file;
	struct JitterSychro out;
	char line[500];
	const char sep[1] = " ";
	char *token;
	double *rec_times;
	short rec_count = 0;

	rec_times = malloc(sizeof(double) * 50);

	file = fopen("/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/jitterdebug.txt",
			"r");

	fgets(line, sizeof(line), file);

	token = strtok(line, sep);
	out.new_time = atof(token);
	token = strtok(NULL, sep);
	out.accept = (short) atoi(token);
	token = strtok(NULL, sep);
	while (token != NULL) {
		rec_times[rec_count++] = atof(token);
		token = strtok(NULL, sep);
	}

	out.rec_times = realloc(rec_times, sizeof(double) * rec_count);
	out.length = rec_count;
	fclose(file);

	return out;
}
