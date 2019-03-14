/*
 * debugging.h
 *
 *  Created on: 27.10.2016
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

#ifndef DEBUGGING_H_
#define DEBUGGING_H_

long findSynchronisedSegment(struct SamplingSet set, struct Smc path);
short checkNodeOrdering(struct Tree t);
short checkNodeTimeOrder(struct Tree t);
short checkOperations(struct Smc path);
short checkRecombinationTimes(struct Smc path);
short checkTopology(struct Tree t);
short checkTree(struct Tree t) ;
short checkTreePathCompletely(const struct Smc path);
short checkTreePathStructure(const struct Smc path);
short countNewNodes(struct Tree t_curr, struct Tree t_next);
struct AcceptanceSynchro synchroniseAcceptance() ;
struct FreeTimeSynchro synchroniseFreeTimes();
struct InitTimeSynchro synchroniseInitTimeProposal();
struct RecTimeSynchro synchroniseRecombinationTimes();
struct Smc synchronisePath(struct SegmentSamplerSynchro sss);
void testExponentialDisotribution();
void deallocateSegmentSamplingSynchro(struct SegmentSamplerSynchro sss);
#endif /* DEBUGGING_H_ */
