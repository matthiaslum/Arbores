/*
 * timeadjustment.h
 *
 *  Created on: 8.11.2016
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

#ifndef TIMEADJUSTMENT_H_
#define TIMEADJUSTMENT_H_

#include "datastructures.h"

short adjustTimes(struct Smc path, struct Tree tr);
short checkFixedTimeCondition(struct Smc path, struct Tree tr);
short findConditioningIndex(short *gis, short gi, short nb);
short hasExtraRecombinations(struct Smc path, short *extra, struct Conditioning cond);
short is_unfixed(short i, short *glb_idx, short nb);
short secondPass(struct LimitSummary smry, struct Smc path, struct Tree tr);
short findLimitingStep(struct Smc path, short gi);
short *findLimitingIndices(double *times, double new_time, short nb);
struct LimitData getLimiters(struct Tree t, double new_time, short step, struct Smc path);
struct LimitSummary firstPass(struct Smc path, struct Tree tr);
struct DoubleVector postBounders(struct Smc path, struct LimitSummary smry, short step,
		short push_dir, short *is_limiter);
struct SamplingSet recreatePath(struct Smc path);
struct ShortVector topologicalMatching(struct Tree *trees, long n_trees, struct Tree t);
void assignGlobalIndices(short **global_idx, short *n_global, short step,
		struct Tree *trees);
void assignFreeTimeIndicators(struct Smc path, struct Conditioning cond);
void deallocateSummary(struct LimitSummary summary);
void forwardTimeAdjustment(struct Smc path, short step, double trgt, short idx, short gi);
void printLimitSummary(struct LimitSummary summary);
void recreationRecursion(struct SimplePath path, struct Smc orig_path,
		struct SamplingSet set);
long timeAdjustment(struct LinkedSetArray array, struct Conditioning cond,
		struct Data data, struct SamplingSet sampling_set, short *extra, short verb);
void updateLimiters(struct LimitSummary smry, short step, double trgt, short nl, short gi);

#endif /* TIMEADJUSTMENT_H_ */
