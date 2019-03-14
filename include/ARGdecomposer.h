/*
 * ARGdecomposer.h
 *
 *  Created on: 12.10.2016
 *      Author: heine
 *
 * Written by Kari Heine.
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

#ifndef ARGDECOMPOSER_H_
#define ARGDECOMPOSER_H_

#include <math.h>
#include "datastructures.h"
#include "treeutils.h"

int decomposeARG(struct ShrubARG arg, struct Tree *tree_path, short **opers,
		struct Data data, short n_nodes, int *op_sites);
void binarize(short *tree, short *opers, short n_nodes, short n_leaves, short *times);
void shiftIndices(short *tree, short *oper, short n_nodes, short i, short n_alloc);
int jumpGeneration(short *tree, short *parent_rows, short i, short *ch, short n_nodes);
void extractRecombinationSites(struct ShrubARG arg, struct Data data, int *op_sites);
void extractTreeSequence(struct ShrubARG ARG, struct Data data, short **trees,
		int *op_sites);
void extractTree(short *C, int site, short n_nodes, short n_leaves);
void extractOperations(struct ShrubARG arg, short **opers, int *op_sites);
void createARGtimes(struct ShrubARG arg, short *times);
void setTime(short *C, short* times, short n, short n_nodes, short n_leaves);
void resolveParallelRecombinations(int *op_sites, int n_op, struct Data data);


#endif /* ARGDECOMPOSER_H_ */
