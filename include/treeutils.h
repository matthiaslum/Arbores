/*
 * treeutils.h
 *
 *  Created on: 5.10.2016
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

#ifndef TREEUTILS_H_
#define TREEUTILS_H_

#include "datastructures.h"

int compareTreeMtxs(short *T1, short *T2, short n1, short n2);
enum Node_Colors *colorTree(struct Tree t, short *M, short col);
short checkFeasibility(struct Tree t, short *M, int col);
short checkPath(struct SimplePath path, struct Data data);
short countChildren(short *tree, short i, short n_rows, short *ch);
short determineNewNode(struct Tree t_curr, struct Tree t_next);
short findMatchingNode(struct Tree tree, double t);
short findParentRows(short *C, short n_rows, short node, short *parent_rows);
short findParent(short *C, short n_rows, short n_leaves, short node);
short findRoot(short *C, short n_nodes, short n_leaves);
short SPR(struct Tree t, const short *o);
short *calculateOperation(struct Tree t, struct Tree t_next);
short *createOperationCopy(short *op);
short *recombinationSiteIndicator(struct Smc path);
struct BranchLength branchLengths(struct Tree t, struct Data data, long col);
struct ShortVector findAncestors(struct Tree t, short node);
//void copyTreeToTreeArray (struct Tree t, struct Tree_array_version * tree_array);
//struct Tree createCopyFromTreeArray(struct Tree_array_version t);
struct Tree createCopy(struct Tree t);
struct Tree createTree(short n_leaves);
struct Tree duplicateTree(struct Tree t);
struct Tree treeFromData(short *C, double *times, short n_nodes, short m);
void canonise(struct Tree t);
void canonicalOrdering(const short *C, int *order, short n_nodes);
void computeBranchwidth(struct Tree t, short node, short *bw);
void copyTree(struct Tree *t_dst, struct Tree t_src);
void deletePath(struct Smc smc);
void deleteTree(struct Tree t);
void getParentCounts(short *C, short n_nodes, short *parent_cnts);
void leafXCoordinates(struct Tree t, double *x);
void printARG(struct ShrubARG ARG);
void printTreePath(struct Tree *t, short **opers, short path_len, int *op_sites);
void printTree(struct Tree t);
void reorderSparseMtx(short *C, const int *order, short n_nodes);
void sortChildren(struct Tree *trees, long length);
void timeSorting(struct Tree t, short *o);

#endif /* TREEUTILS_H_ */
