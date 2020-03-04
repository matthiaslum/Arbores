/*
 * MCMCutils.h
 *
 *  Created on: 21.10.2016
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

#ifndef MCMCUTILS_H_
#define MCMCUTILS_H_

double calculateCardinalityRatio(struct SamplingSetData ss_data, long pick);
long drawPathIndex(struct SamplingSet set);
long findCurrentSegment(struct Smc path, struct SamplingSet set);
short checkAdjacencySets(struct AdjacencySets a, struct Conditioning cond, int j);
short checkCompatibility(struct Smc path, struct Data data);
short checkLastSet(struct AdjacencySets a, struct LinkedSet *l, struct Data d, int j);
short *canonisePath(struct Smc path);
void shareWorkloadAdjacencySets(long domain_size);
short *extract2DOperations(short *op3D, long m, short length);
short *noopOperExtension();
struct AdjacencySet uniqueOperationSequences(struct AdjacencySet aSet);
struct AdjacencySets createAdjacencySetsRecursively(struct Tree *domain, long domain_size,
		struct Conditioning cond, int j, short *forced);
struct AdjacencySets createFinalAdjSets(struct Tree *domain, long domain_size,
		struct Conditioning cond, int j);
struct AdjacencySets createNonTrivialAdjSets(struct Tree *domain, long domain_size,
		struct Conditioning cond, int j);
struct AdjacencySets createNoopAdjacencySets(struct Conditioning cond,
		struct Tree *domain, long domain_size, int j, short n_l);
struct AdjacencySets createTrivialAdjSets(struct Tree *domain, long domain_size,
		struct Conditioning cond, int j);
struct CanonicalMatrix canoniseTreeSet(struct Tree *trees, long n_trees);
struct CanonicalMatrix uniqueRows(struct ShortRowSortedMtx canon);
struct Conditioning prepareConditioning(int sgmnt, struct BridgePoints bps, struct Data d,
		struct Smc path);
struct Data augmentWithNonSegregatingSites(struct Data d, struct Smc path);
struct FindResult findSegmentSites(struct Data d, struct BridgePoint bp);
struct ForcingConfig createEmergencyForcingConfiguration(int length);
struct ForcingConfig createForcingConfigurations(int length);
struct LinkedSetArray exhaustivePathFinder(struct Conditioning cond, short *force_conf,
		struct Data d);
struct LinkedSet *linkSets(struct AdjacencySets A);
struct LongVector getSites(struct Conditioning cond);
double calculateAlphaModified(struct Smc segment_p, struct Smc segment_c, struct Smc segment_c_with_times, struct Data data, struct Parameters parm,
                              struct Conditioning condition, double extra);
struct MCMCDiagnostics calculateAlpha(struct Smc segment_p, struct Smc path_p, struct Smc segment_c,
		struct Smc path_c, struct Data data, struct Parameters parm,
		struct Conditioning condition, double extra);
struct MCMCSummary segmentSampler(struct Smc path_c, int seg, struct BridgePoints bps,
		struct Data data, struct Parameters parm);
struct segment_output truncatedSegmentSampler(struct Smc path_c, int seg, struct BridgePoints bps,
                                              struct Data data, struct Parameters parm);
struct SamplingSetData generateSamplingSet(struct Smc segment_c, struct Conditioning cond,
		struct Data d);
struct ShortVector createTreeSelector(struct Smc path, struct Data d);
struct Smc createCompleteProposal(struct Smc segment, struct Smc path,
		struct Parameters parm, struct Conditioning condition, struct Data data);
struct Conditioning_array_version copyConditioningToArray(struct Conditioning condition);
struct Conditioning copyConditioning(struct Conditioning_array_version condition);
struct arraySegmentOutput convertSegmentOutputToArray (struct segment_output output,  struct arraySegmentOutput * array_version);
struct Smc combineSegments( struct arraySegmentOutput *temp_summary, int len, struct Data data);
struct Smc extendToFullSequence(struct Smc segment, struct Smc path,
		struct Conditioning cond, struct Data data);
struct Smc extractSegment(struct Smc path, struct Conditioning cond);
struct Smc extractAndReverseSegment(struct Smc path, struct Conditioning cond);
struct Tree *noopTreeExtension(struct Tree t);
struct LongVector equivalence(struct ShortRowSortedMtx cmtx, short *entry, long start_row);
struct LinkedSet plainLinkedSet(struct AdjSetUnion u, long *sites);
struct LinkedSet linkedSetWithEquivalenceClasses(struct AdjSetUnion u, long *sites,
		short offset);
struct Tree* getTerminalTrees(struct AdjSetUnion u);
struct SegmentSamplerSynchro synchroniseSegmentSampler();
struct ShortRowSortedMtx sortedCanonMtxOfTerminalTrees(struct AdjSetUnion u);
struct ShortRowSortedMtx sortCanonicalMtx(struct CanonicalMatrix canon_mtx);
struct AdjacencySets createEmptyAdjacencySets(long domain_size);
struct BridgePoints createPhaseOneBridgePoints(struct Data d);
struct BridgePoints createPhaseTwoBridgePoints(struct Data d);
struct BridgePoints createPhaseThreeBridgePoints(struct Data d);
struct BridgePoints createBridgePoints(struct Data d, const int len);
struct Smc reversePathSegment(struct Smc path);
struct UniqueRowsWithCounts uniqueRowsWithCounts(struct ShortRowSortedMtx canon);
void calculateSetSizes(struct SamplingSet set);
void deallocateDomain(struct Tree *domain, long domain_size);
void deallocateForcingConfiguration(struct ForcingConfig config);
void printForcingConfigurations(struct ForcingConfig fc);
void printAdjSets(struct AdjacencySets a);
void jitterStep();
void deleteAdjacencySets(struct AdjacencySets a);
void deleteLinkedSetArray(struct LinkedSetArray ls);
void deallocateSamplingSet(struct SamplingSet s);
void printSegmentPreamble(struct Conditioning condition, struct Data data);

#endif /* MCMCUTILS_H_ */
