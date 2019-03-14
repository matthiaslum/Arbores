/*
 * exhaustiveSearch.h
 *
 *  Created on: 24.10.2016
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

#ifndef EXHAUSTIVESEARCH_H_
#define EXHAUSTIVESEARCH_H_

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "datastructures.h"
#include "treeutils.h"
#include "graph2tikz.h"
#include "debugging.h"

struct AdjacencySet exhaustiveSearch(struct Tree t, struct Conditioning cond,
		short site_offset, short forced, short length);
struct SearchPath createSearchPath(struct Tree t, short length);
void deleteSearchPath(struct SearchPath path);
void step(enum Step_Type st, struct SearchPath path, struct AdjacencySet a,
		struct Conditioning cond, short col, short max_len, short extra);
struct OperationSet mergeSearchSpace(struct SearchPath path, struct Conditioning cond,
		short col, short extra_rec);
struct RegraftTimes getRegraftTimeRange(struct Tree t, short *op);
struct OperationSet pruneSearchSpace(struct SearchPath path, struct Conditioning cond,
		short col);
struct Roots findSingleColorSubtreeRoots(struct Tree t, enum Node_Colors *colors,
		enum Node_Colors clr);
//enum Node_Colors *colorTree(struct Tree t, short *M, short col);
struct Descendants descendants(struct Tree t, short node, struct Descendants descs);
short *subtreeComplement(struct Descendants desc, short n_n);
struct Removed_Nodes removeRecursion(struct Tree t, struct Removed_Nodes remove);
short isroot(struct Tree t, short node, struct Roots r);
short *pickOp(short *ops, long row, long n_ops);
short checkOperationValidity(struct Tree t, short *o);
short determineNewNode(struct Tree t_curr, struct Tree t_next);
short checkNodeOrdering(struct Tree t);
short checkNodeTimeOrder(struct Tree t);

#endif /* EXHAUSTIVESEARCH_H_ */
