/*
 * pathutils.h
 *
 *  Created on: 17.11.2016
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

#ifndef PATHUTILS_H_
#define PATHUTILS_H_

#include "datastructures.h"

double *calculateRegraftTimes(struct Tree *t, short path_len);
double *createCoalescenceTimes(short nb, struct Parameters parm);
enum Node_Colors **colorPath(struct Smc path, struct Data data);
long *segmentEndIndices(struct Smc path, struct Conditioning cond);
short checkGlobalOrderConsistency(struct Smc path1, struct Smc path2);
short countRecombinations(struct Smc path);
short *createNoop();
struct Smc_array_version convertPathToArray(struct Smc path, struct Smc_array_version * array_version);
struct Smc convertPathArrayToSmc (struct Smc_array_version array_version, struct Data data);
struct MRCA timesToMRCA(struct Smc path);
struct SimplePath createNewPath(struct Smc path);
struct Smc createPathCopy(const struct Smc path);
struct Smc createPathCopy2(const struct Smc path);
struct Smc createPath(long length);
struct Smc generateTimes(struct Smc path, struct Parameters parm);
struct Smc extractToSites(struct Smc path, struct LongVector sites);
struct Smc removeNoOps(struct Smc path);
struct LongVector siteUnion(struct LongVector sites, struct Smc path);
void deallocateMRCA(struct MRCA mrca);
void deallocatePath(struct Smc path);
void printGlobalIndices(struct Smc path);
void printOperations(struct Smc path);
void printRecombinationTimes(struct Smc path);
void printTreeSelector(struct Smc path);

#endif /* PATHUTILS_H_ */
