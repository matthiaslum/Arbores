/*
 * smcPrior.h
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

#ifndef SMCPRIOR_H_
#define SMCPRIOR_H_

#include "datastructures.h"

struct NonIdentityOps getNonIdentityOperations(struct Smc path);
struct SmcPriorData smcprior(struct Smc path, struct Parameters parm, struct Data data);
struct RegraftTimeData rgTimePdf(double x, struct Tree t, double rec_time,
		double pr_prnt_time, double *br_lengths, struct Parameters parm);
double evalDensity(double x, double n_act, double *t_lo_trnc, short hit, double n_eff,
		short *active_counts, double *dt);
short activeBranchCount(double x, short nl, double *branch_lengths, double *times,
		double pr_prnt_time);
void deallocateOperations(struct NonIdentityOps ops);
void deallocatePriorData(struct SmcPriorData priordata);
void testrgpdf(struct Tree t, double rec_time, double pr_prnt_time, double *lengths,
		struct Parameters parm);
#endif /* SMCPRIOR_H_ */
