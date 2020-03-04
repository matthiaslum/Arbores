/*
 * treeutils.c
 *
 *  Created on: 13.11.2016
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
#include <assert.h>
#include <limits.h>
#include <math.h>
#include "datastructures.h"
#include "sorting.h"
#include "debugging.h"
#include "data.h"
#include "constants.h"
#include "graph2tikz.h"
#include "treeutils.h"
#include "utils.h"

int time_i_sort[15];
int time_i_sort_inv[15];
short C_tmp[2*15];
double times_tmp[15];

short branch_width[2*15-1];
short selected_branch_width[2*15-1];

short *createOperationCopy(short *op) {
	short *out;
	out = malloc(sizeof(short) * 2);
	out[0] = op[0];
	out[1] = op[1];
	return out;
}

/*
 * NOTE: This function returns the indices of the parent rows. If searching parents
 * in a tree matrix not including rows for the leaves, one must add the number of
 * leaves to the output to get the right parent labels.
 *
 * C			- The graph matrix
 * n_rows	-	number of rown in C
 * node		- the label of the node whos parents are looked for
 * parent_rows - the row indices where parents are found
 */
short findParentRows(short *C, short n_rows, short node, short *parent_rows) {

	short j = 0;
	for (int i = 0; i < n_rows; i++)
    if (C[i] == node || C[i + (int) n_rows] == node) {
			parent_rows[j++] = (short) i;
    }
	return j;
}

short findParent(short *C, short n_rows, short n_leaves, short node) {
	for (int i = 0; i < n_rows; i++)
		if (C[i] == node || C[i + (int) n_rows] == node)
			return (short) (i + n_leaves);
	return -1;
}

/* C 			- The array defining the graph
 * n_rows - The number of rows in C
 * node 	- Label of the node whose parents are looked for */
void printARG(struct ShrubARG ARG) {

	int collab[] = { 0, 1, 2, 4 };
	int *rowlab, *Cprint;
	short *parent_cnts, *C;
	short n_branch = ARG.n_nodes - ARG.n_seq;

	rowlab = malloc(sizeof(int) * ARG.n_nodes);
	for (int i = 0; i < ARG.n_nodes; i++)
		rowlab[i] = i;

	C = malloc(sizeof(short) * 3 * n_branch);
	for (int i = 0; i < 3 * n_branch; i++)
		C[i] =
				ARG.C[i + ARG.n_seq + ((int) floor((double) i / (double) n_branch)) * ARG.n_seq];

	/* Get the number of parents for each node */
	parent_cnts = malloc(sizeof(short) * ARG.n_nodes);
	getParentCounts(C, ARG.n_nodes, parent_cnts);
	free(C);

	Cprint = malloc(sizeof(int) * 4 * ARG.n_nodes);
	for (int i = 0; i < ARG.n_nodes; i++)
		for (int j = 0; j < 3; j++)
			Cprint[i + j * ARG.n_nodes] = (int) ARG.C[i + j * ARG.n_nodes];

	for (int i = 0; i < ARG.n_nodes; i++)
		Cprint[i + 3 * ARG.n_nodes] = (int) parent_cnts[i];

	printf("ARG with:\n\t%i recombinations\n\t%i leaves\n\t%i nodes\n", ARG.n_rec,
			ARG.n_seq, ARG.n_nodes);
	printIntWithLabelsW(Cprint, collab, rowlab, 6, 4, ARG.n_nodes, ARG.n_nodes);

	free(rowlab);
	free(parent_cnts);
	free(Cprint);

}

short countParents(short *C, short n_rows, short node) {

	short j = 0;
	for (int i = 0; i < n_rows; i++)
		if (C[i] == node || C[i + (int) n_rows] == node)
			j++;

	return j;
}

/* Returns the number of children and the child indices in the argument
 * array ch.
 * */
short countChildren(short *tree, short i, short n_rows, short *ch) {

	short n_childs = 0;

	if (tree[i] >= 0)
		ch[n_childs++] = tree[i];
	if (tree[i + n_rows] >= 0)
		ch[n_childs++] = tree[i + n_rows];

	return n_childs;
}

void getParentCounts(short *C, short n_nodes, short *parent_cnts) {
	for (int i = 0; i < n_nodes; i++)
		parent_cnts[i] = countParents(C, n_nodes, i);
}

short findRoot(short *C, short n_nodes, short n_leaves) {

	int i;
	for (i = n_leaves; i < n_nodes; i++)
		if (countParents(C, n_nodes, i) == 0)
			break;
	return i;
}

void printTree(struct Tree t) {

	short n_branch, n_leaves;
	n_leaves = (t.n_nodes + 1) / 2;
	n_branch = n_leaves - 1;
	for (int i = 0; i < n_branch; i++)
		printf("%4i: %4i %4i %2.12f %p\n", i + n_leaves, t.C[i], t.C[i + n_branch],
				t.times[i], &t.times[i]);
	printf("\n");
}

void printAdjacencySets(struct AdjacencySets a) {

	for (int i = 0; i < a.n_sets; i++)
		if (a.set_size[i] > 0) {
			for (int k = 0; k < a.set_size[i]; k++)
				for (int j = 0; j < a.length; j++) {
					printf("Set %i, Entry %i, Step %i:\n", i, k, j);
					printTree(a.tree_adj_set[i][k + j * a.set_size[i]]);
				}
		}
}

void printLinkedSets(struct LinkedSet *ls, short n_sets) {
	for (int i = 0; i < n_sets; i++) {
		assert(ls[i].length == 2);
		printf("Trees after step %i\n", i);
		for (int j = 0; j < ls[i].n_trees; j++) {
			printTree(ls[i].trees[j + ls[i].n_trees]);
		}
	}
}

void printTreePath(struct Tree *t, short **opers, short path_len, int *op_sites) {

	short n_branch, n_leaves;
	n_leaves = (t[1].n_nodes + 1) / 2;
	n_branch = n_leaves - 1;

	if (op_sites != NULL) {
		for (int i = 0; i < path_len; i++)
			printf("Site: %4i               |", op_sites[i]);
		printf("\n");
	}
	if (opers != NULL) {
		for (int i = 0; i < path_len - 1; i++)
			printf("SPR: [%4i, %4i]        |", opers[i][0], opers[i][1]);
		printf("\n");
	}
	for (int j = 0; j < n_branch; j++) {
		for (int i = 0; i < path_len; i++)
			printf("%4i: %4i %4i %8.2f |", j + n_leaves, t[i].C[j], t[i].C[j + n_branch],
					t[i].times[j]);
		printf("\n");
	}
}

enum Node_Colors *colorTree(struct Tree t, short *M, short col) {

	enum Node_Colors *colors;
	enum Node_Colors L, R;
	short n_l, n_b;

	n_l = (t.n_nodes + 1) / 2;
	n_b = n_l - 1;

	colors = malloc(sizeof(enum Node_Colors) * t.n_nodes);

	for (int i = 0; i < n_l; i++) // leaves
		colors[i] = M[i + col * n_l] == 1 ? black : white;

	for (int i = 0; i < n_b; i++) {

		L = t.C[i] >= 0 ? colors[t.C[i]] : undef;
		R = t.C[i + n_b] >= 0 ? colors[t.C[i + n_b]] : undef;

		if ((L == black && R == black) || (L == black && R == undef)
				|| (L == undef && R == black))
			colors[i + n_l] = black; // black
		else if ((L == white && R == white) || (L == white && R == undef)
				|| (L == undef && R == white))
			colors[i + n_l] = white; // white
		else
			colors[i + n_l] = gray; // grey
	}
	return colors;
}

/* NB: Returns the index of the new node in the range [0,...,nb].
 * If not found returns -1;
 * */
short determineNewNode(struct Tree t_curr, struct Tree t_next) {

	short nb = (t_curr.n_nodes + 1) / 2 - 1;
	double time_next;
	int j;

	for (int i = 0; i < nb; i++) {
		time_next = t_next.times[i];
		for (j = 0; j < nb; j++)
			if (fabs(time_next - t_curr.times[j]) < (time_next + t_curr.times[j]) * TOL)
				break;
		if (j == nb)
			return i;
	}
	return -1;
}

int compareTreeMtxs(short *T1, short *T2, short n1, short n2) {

	if (n1 != n2) // Numbers of nodes don't not match
		return 0;

	for (int i = 0; i < n1; i++)
		if (T1[i] != T2[i] || T1[i + n1] != T2[i + n1]) // graph matrices dno't match
			return 0;
	return 1;

}

void canonise(struct Tree t) {

	double *new_times;
	int *order;
	short n_leaves, n_branch;

	n_leaves = (t.n_nodes + 1) / 2;
	n_branch = n_leaves - 1;
	new_times = malloc(sizeof(double) * n_branch);
	order = malloc(sizeof(int) * n_branch);

	for (int i = 0; i < n_branch; i++)
		order[i] = i;

	canonicalOrdering(t.C, order, t.n_nodes);

// Reorder times according to the ordering
	for (int i = 0; i < n_branch; i++)
		new_times[i] = t.times[order[i]];

// Reorder the adjacency matrix
	reorderSparseMtx(t.C, order, t.n_nodes);

// Replace the reordered contents in T
	for (int i = 0; i < n_branch; i++)
		t.times[i] = new_times[i];

	free(order);
	free(new_times);

}

void canonicalOrdering(const short *C, int *order, short n_nodes) {

	short *array;
	int n_cols, n_rows;
	short nl, nb, bw_left, bw_right, se_left, se_right;
	nl = (n_nodes + 1) / 2;
	nb = nl - 1;

	// initialise the array for sorting
	array = malloc(sizeof(short) * 2 * nb);
	n_cols = 2;
	n_rows = nb;

	for (int i = 0; i < nb; i++) {

		if (C[i] < nl) {
			bw_left = 1;
			se_left = C[i];
		} else {
			bw_left = array[C[i] - nl];
			se_left = array[C[i] - nl + nb];
		}

		if (C[i + nb] < nl) {
			bw_right = 1;
			se_right = C[i + nb];
		} else {
			bw_right = array[C[i + nb] - nl];
			se_right = array[C[i + nb] - nl + nb];
		}
		// branch widths
		array[i] = bw_left + bw_right;

		// smallest elements
		array[i + nb] = se_left < se_right ? se_left : se_right;
	}

	shortSortRows(array, order, 0, n_rows, n_cols, n_rows);

	free(array);
}

void reorderSparseMtx(short *C, const int *order, short n_nodes) {
	short *h, *tmp;
	short n_l, n_b, v1, v2;
	n_l = (n_nodes + 1) / 2;
	n_b = n_l - 1;

	tmp = malloc(sizeof(short) * 2 * n_b);
	h = malloc(sizeof(short) * n_b);
	for (int i = 0; i < n_b; i++)
		h[order[i]] = (short) i;

	for (int i = 0; i < n_b; i++) {
		if (C[i] >= n_l)
			tmp[i] = h[C[i] - n_l] + n_l;
		else
			tmp[i] = C[i];
		if (C[i + n_b] >= n_l)
			tmp[i + n_b] = h[C[i + n_b] - n_l] + n_l;
		else
			tmp[i + n_b] = C[i + n_b];
	}
	free(h);
	for (int i = 0; i < n_b; i++) {
		v1 = tmp[order[i]];
		v2 = tmp[order[i] + n_b];
		if (v1 < v2) {
			C[i] = v1;
			C[i + n_b] = v2;
		} else {
			C[i] = v2;
			C[i + n_b] = v1;
		}
	}
	free(tmp);
}

struct Tree createTree(short n_leaves) {

	short n_b = n_leaves - 1;
	struct Tree t;
	t.C = malloc(sizeof(short) * 2 * n_b); //All the left leaves first, followed by the right leaves.
	t.times = malloc(sizeof(double) * n_b);
	t.n_nodes = 2 * n_leaves - 1;
	return t;
}

void deleteTree(struct Tree t) {
	free(t.C);
	free(t.times);
}

struct Tree duplicateTree(struct Tree t) {
	struct Tree out;
	short n_l = (t.n_nodes + 1) / 2;
	short n_b = n_l - 1;
	out = createTree(n_l);
	shortArrayCopy(out.C, t.C, n_b, 2, n_b);
	doubleArrayCopy(out.times, t.times, n_b, 1, n_b);
	return out;
}

short checkFeasibility(struct Tree t, short *M, int col) {

	short *bw, *sbw;
	short n_n, n_l, n_b, n_mutated = 0, out;
	int j;
	n_n = t.n_nodes;
	n_l = (n_n + 1) / 2;
	n_b = n_l - 1;

//	bw = malloc(sizeof(short) * n_n); // branch widths
//	sbw = malloc(sizeof(short) * n_n); // branch widths (mutated leaves only)
	bw = branch_width;
	sbw = selected_branch_width;

	for (int i = 0; i < n_l; i++) {
		sbw[i] = M[i + col * n_l];
		bw[i] = 1;
		n_mutated += sbw[i];
	}
	if (n_mutated <= 1) {
//		free(bw);
//		free(sbw);
		return 1;
	}

	for (j = 0; j < n_b; j++) {
		sbw[j + n_l] = sbw[t.C[j]] + sbw[t.C[j + n_b]];
		bw[j + n_l] = bw[t.C[j]] + bw[t.C[j + n_b]];
		if (sbw[j + n_l] >= n_mutated)
			break;
	}

	out = sbw[j + n_l] == bw[j + n_l] ? 1 : 0;

//	free(bw);
//	free(sbw);
	return out;
}

void timeSorting(struct Tree t, short *o) {

	short nl, nb;
	nl = (t.n_nodes + 1) / 2;
	nb = nl - 1;

	for (int i = 0; i < nb; i++)
		time_i_sort[i] = i;
	sortDoubleArray(time_i_sort, t.times, nb, 0);

	for (int i = 0; i < nb; i++) {
		C_tmp[i] = t.C[i];
		C_tmp[i + nb] = t.C[i + nb];
		times_tmp[i] = t.times[i];
		time_i_sort_inv[time_i_sort[i]] = i;
	}

// rearrange the rows
	for (int i = 0; i < nb; i++) {
		t.C[i] = C_tmp[time_i_sort[i]];
		t.C[i + nb] = C_tmp[time_i_sort[i] + nb];
		t.times[i] = times_tmp[time_i_sort[i]];
	}

// update the matrix entries
	for (int i = 0; i < nb; i++) {
		t.C[i] = t.C[i] >= nl ? time_i_sort_inv[t.C[i] - nl] + nl : t.C[i];
		t.C[i + nb] = t.C[i + nb] >= nl ? time_i_sort_inv[t.C[i + nb] - nl] + nl : t.C[i + nb];
	}

// update the operation, if provided
	if (o != NULL) {
		o[0] = o[0] >= nl ? time_i_sort_inv[o[0] - nl] + nl : o[0];
		o[1] = o[1] >= nl ? time_i_sort_inv[o[1] - nl] + nl : o[1];
	}

	sortChildren(&t, 1);
}

struct Tree treeFromData(short *C, double *times, short n_nodes, short m) {
	struct Tree out;
	short n_l = (n_nodes + 1) / 2;
	short n_b = n_l - 1;
	out = createTree(n_l);
	shortArrayCopy(out.C, C, n_b, 2, m);
	doubleArrayCopy(out.times, times, n_b, 1, m);
	return out;
}

void copyTree(struct Tree *t_dst, struct Tree t_src) {
	short n_b = (t_src.n_nodes + 1) / 2 - 1;
	for (int i = 0; i < n_b; i++) {
		t_dst->C[i] = t_src.C[i];
		t_dst->C[i + n_b] = t_src.C[i + n_b];
		t_dst->times[i] = t_src.times[i];
	}
	t_dst->n_nodes = t_src.n_nodes;
}


void copyTreeToTreeArray (struct Tree t, struct Tree_array_version * tree_array){
    short n_b = (t.n_nodes + 1) / 2 - 1;
    tree_array->n_nodes = t.n_nodes;
    for (int j = 0; j < n_b; j++) {
        tree_array->C[j] = t.C[j];
        tree_array->C[j + n_b] = t.C[j + n_b];
        tree_array->times[j] = t.times[j];
    }
}

//version 2 assigns a new tree array in memory whereas version 1 uses a reference.
struct Tree_array_version copyTreeToTreeArray_v2 (struct Tree t){
    struct Tree_array_version out;
    short n_b = (t.n_nodes + 1) / 2 - 1;
    out.n_nodes = t.n_nodes;
    for (int j = 0; j < n_b; j++) {
        out.C[j] = t.C[j];
        out.C[j + n_b] = t.C[j + n_b];
        out.times[j] = t.times[j];
    }
    return out;
}


struct Tree createCopyFromTreeArray(struct Tree_array_version t) {

    struct Tree out;
    short n_b = (t.n_nodes + 1) / 2 - 1;

    out.C = malloc(sizeof(short) * n_b * 2);
    out.times = malloc(sizeof(double) * n_b);
    out.n_nodes = t.n_nodes;

    for (int i = 0; i < n_b; i++) {
        out.C[i] = t.C[i];
        out.C[i + n_b] = t.C[i + n_b];
        out.times[i] = t.times[i];
    }
    return out;
}

struct Tree createCopy(struct Tree t) {

	struct Tree out;
	short n_b = (t.n_nodes + 1) / 2 - 1;

	out.C = malloc(sizeof(short) * n_b * 2);
	out.times = malloc(sizeof(double) * n_b);
	out.n_nodes = t.n_nodes;

	for (int i = 0; i < n_b; i++) {
		out.C[i] = t.C[i];
		out.C[i + n_b] = t.C[i + n_b];
		out.times[i] = t.times[i];
	}
	return out;
}


void computeBranchwidth(struct Tree t, short node, short *bw) {

	short n_b, n_l;
	n_l = (t.n_nodes + 1) / 2;
	n_b = n_l - 1;

	if (node < 0)
		bw[node] = 0;
	else if (node < n_l)
		bw[node] = 1;
	else {
		computeBranchwidth(t, t.C[node - n_l], bw);
		computeBranchwidth(t, t.C[node - n_l + n_b], bw);
		bw[node] = bw[t.C[node - n_l]] + bw[t.C[node - n_l + n_b]];
	}
}

/* This is a low level auxiliary function for calculating x-coordinates for
 * nodes in a tree for nice graphical output.
 */
void childXCoordinates(struct Tree t, short node, double *x, short *bw, double dist) {

	short bw_l, bw_r, nl, nb;
	nl = (t.n_nodes + 1) / 2;
	nb = nl - 1;

	bw_l = bw[t.C[node - nl]];
	bw_r = bw[t.C[node - nl + nb]];

	if (t.C[node - nl] >= 0)
		x[t.C[node - nl]] = x[node] - dist;
	else
		x[t.C[node - nl]] = x[node];
	if (t.C[node - nl + nb] >= 0)
		x[t.C[node - nl + nb]] = x[node] + dist;
	else
		x[t.C[node - nl + nb]] = x[node];

	if (t.C[node - nl] >= nl)
		childXCoordinates(t, t.C[node - nl], x, bw, dist / (double) 2);
	if (t.C[node - nl + nb] >= nl)
		childXCoordinates(t, t.C[node - nl + nb], x, bw, dist / (double) 2);

}

/* This function computes coordinates for the leaf nodes in a manner that the
 * resulting tree does not have crossing branches.
 * */
void leafXCoordinates(struct Tree t, double *x) {

	short nl;
	short *bw;
	int *order;
	double dist = 1;

	bw = malloc(sizeof(short) * t.n_nodes);
	computeBranchwidth(t, t.n_nodes - 1, bw);
	x[t.n_nodes - 1] = 0;
	childXCoordinates(t, t.n_nodes - 1, x, bw, dist);
	free(bw);

	nl = (t.n_nodes + 1) / 2;
	order = malloc(sizeof(int) * nl);
	for (int i = 0; i < nl; i++)
		order[i] = i;
	sortDoubleArray(order, x, nl, 0);

	for (int i = 0; i < nl; i++)
		x[order[i]] = (double) i;

	free(order);

}

short SPR(struct Tree t, const short *o) {

	short pr, rg, pr_prnt, nl, nb, sibling, pr_gparent, rg_parent, i_ins, tmp;

	nl = (t.n_nodes + 1) / 2;
	nb = nl - 1;

	pr = o[0];
	rg = o[1];

// no-op
	if (pr == -1 && rg == -1)
		return -1;

	pr_prnt = findParent(t.C, nb, nl, pr);

	if ((pr_gparent = findParent(t.C, nb, nl, pr_prnt)) != -1) {
		// if prune node has a grand parent, attach it to the prune node's sibling

		sibling = t.C[pr_prnt - nl] == pr ? t.C[pr_prnt - nl + nb] : t.C[pr_prnt - nl];
		if (t.C[pr_gparent - nl] == pr_prnt)
			t.C[pr_gparent - nl] = sibling;
		else
			t.C[pr_gparent - nl + nb] = sibling;
	}

// remove prune parent
	for (int i = pr_prnt - nl; i < nb - 1; i++) {
		t.C[i] = t.C[i + 1];
		t.C[i + nb] = t.C[i + 1 + nb];
		t.times[i] = t.times[i + 1];
	}

// downshift indices
	for (int i = 0; i < 2 * nb; i++)
		if (t.C[i] >= pr_prnt)
			t.C[i]--;
	rg = rg > pr_prnt ? rg - 1 : rg;
	i_ins = rg < pr ? pr - nl + 1 : rg - nl + 1; // new node is inserted before this row
	i_ins = i_ins < 0 ? 0 : i_ins;

// regraft (insert)
	rg_parent = findParent(t.C, nb, nl, rg);
	for (int i = nb - 1; i > i_ins; i--) {
		t.C[i] = t.C[i - 1];
		t.C[i + nb] = t.C[i - 1 + nb];
		t.times[i] = t.times[i - 1];
	}
	t.C[i_ins] = pr;
	t.C[i_ins + nb] = rg;

	if (i_ins > 0) {
		if (i_ins == nb - 1) {
			t.times[i_ins] = t.times[i_ins - 1] + 1;
		} else {
			t.times[i_ins] = (t.times[i_ins] + t.times[i_ins - 1]) / (double) 2;
		}
	} else {
		t.times[i_ins] = t.times[i_ins] / (double) 2;
	}
// upshift indices
	for (int i = 0; i < 2 * nb; i++)
		if (t.C[i] >= i_ins + nl)
			t.C[i]++;
// update the rg parents entry for the inserted node
	if (rg_parent != -1) {
		if (t.C[rg_parent + 1 - nl] == rg) {
			t.C[rg_parent + 1 - nl] = i_ins + nl;
		} else {
			t.C[rg_parent + 1 - nl + nb] = i_ins + nl;
		}
	}

// make sure that first column is smaller than the second
	for (int i = 0; i < nb; i++) {
		tmp = t.C[i];
		if (tmp > t.C[i + nb]) {
			t.C[i] = t.C[i + nb];
			t.C[i + nb] = tmp;
		}
	}

//	if (checkTopology(t) != 1) {
//		printf("Failing tree\n");
//		printTree(t);
//		assert(false);
//	}

// this is the index of the new node
	return i_ins + nl;
}

short findMatchingNode(struct Tree tree, double t) {

	short nb = (tree.n_nodes + 1) / 2 - 1;
	for (int i = 0; i < nb; i++)
		if (fabs(tree.times[i] - t) < (tree.times[i] + t) * TOL)
			return i;

	return -1;
}

void deletePath(struct Smc smc) {
	for (int i = 0; i < smc.path_len - 1; i++)
		free(smc.opers[i]);
	free(smc.opers);
	free(smc.rec_times);
	free(smc.tree_path);
	free(smc.sites);
}

short checkPath(struct SimplePath path, struct Data data) {
	short out = 1;
	short dataCol;
	for (int i = 0; i < path.length; i++) {
		dataCol = getDataColumn(data, path.sites[i]);
		out *= dataCol == -1 ? 1 : checkFeasibility(path.trees[i], data.M, dataCol);
	}
	return out;
}

short structureMatch(struct Tree t1, struct Tree t2) {
	printf("NOT IMPLEMENTED\n");
	assert(false);
}

void sortChildren(struct Tree *trees, long length) {

	short nb, tmp;

	if (length == 0)
		return;

	nb = (trees[0].n_nodes + 1) / 2 - 1;

	for (int i = 0; i < length; i++) {
		for (int j = 0; j < nb; j++) {
			tmp = trees[i].C[j];
			if (tmp > trees[i].C[j + nb]) {
				trees[i].C[j] = trees[i].C[j + nb];
				trees[i].C[j + nb] = tmp;
			}
		}
	}
}

enum Node_Colors **createColorsForPath(struct Tree *trees, struct Data data, int *sites,
		short length) {

	enum Node_Colors **colors;

	colors = malloc(sizeof(enum Node_Colors *) * length);
	for (int i = 0; i < length; i++)
		colors[i] = colorTree(trees[i], data.M, getDataColumn(data, sites[i]));

	return colors;
}

/*
 * Calculates the total branch length and the total branch length of the mutated subtree,
 * if argument 'col' >= 0. Also returns the array of individual branch lengths for each
 * branch node.
 */
struct BranchLength branchLengths(struct Tree t, struct Data data, long col) {

	struct BranchLength branches;
	enum Node_Colors *colors;
	short parent;
	short nl = (t.n_nodes + 1) / 2;
	short nb = nl - 1;
	double increment, t_lower;

	branches.mutated_length = 0;
	branches.total_length = 0;
	branches.lenghts = malloc(sizeof(double) * (t.n_nodes - 1));
	branches.n = t.n_nodes - 1;

	if (col >= 0)
		colors = colorTree(t, data.M, col);
	else {
		branches.mutated_length = -1;
		colors = malloc(sizeof(enum Node_Colors) * t.n_nodes);
	}

	for (int i = 0; i < t.n_nodes - 1; i++) {

		parent = findParent(t.C, nb, nl, i);
		if (parent == -1) {
			writeTikzTexFileForTreePath(NULL, &t, NULL, NULL, NULL, NULL, NULL, NULL, 1);
			printTree(t);
			assert(false);
		}

		t_lower = i < nl ? 0.0 : t.times[i - nl];
		increment = t.times[parent - nl] - t_lower;
		branches.total_length += increment;
		branches.lenghts[i] = increment;
		if (col >= 0 && colors[i] == black)
			branches.mutated_length = increment;
	}
	free(colors);
	return branches;
}

struct ShortVector findAncestors(struct Tree t, short node) {

	struct ShortVector out;
	short *ancestors = malloc(sizeof(short) * t.n_nodes);
	short count = 0, nl = (t.n_nodes + 1) / 2, nb = nl - 1, current_node = node;

	while ((current_node = findParent(t.C, nb, nl, current_node)) != -1)
		ancestors[count++] = current_node - nl;

	out.v = ancestors;
	out.length = count;
	return out;
}

short *recombinationSiteIndicator(struct Smc path) {

	short *rs_ind;

	rs_ind = malloc(sizeof(short) * (path.selector_length));
	fillShortArray(rs_ind, 0, 1, path.selector_length, 1);

	for (int i = 0; i < path.path_len - 1; i++)
		if (path.opers[i][0] != -1 && path.opers[i][1] != -1)
			rs_ind[path.sites[i + 1] - 2] = 1;

	return rs_ind;
}

short *calculateOperation(struct Tree t, struct Tree t_next) {

	short *op;
	short iL, iR, nn_idx, nl, nb, prntL, rg_prnt, ch;

	nl = (t.n_nodes + 1) / 2;
	nb = nl - 1;

	op = malloc(sizeof(short) * 2);

	nn_idx = determineNewNode(t, t_next);
// identity operation
	if (nn_idx == -1) {
		op[0] = -1;
		op[1] = -1;
		return op;
	}

	/* find nodes in t that match the children of the new node */
	ch = t_next.C[nn_idx];
	iL = ch < nl ? ch : findMatchingNode(t, t_next.times[ch - nl]) + nl;
	ch = t_next.C[nn_idx + nb];
	iR = ch < nl ? ch : findMatchingNode(t, t_next.times[ch - nl]) + nl;

	prntL = findParent(t.C, nb, nl, iL);
	rg_prnt = findParent(t_next.C, nb, nl, nn_idx + nl);

	if (prntL == -1) {
		op[0] = iR;
		op[1] = iL;
		return op;
	} else if (rg_prnt == -1) {
		op[0] = iL;
		op[1] = iR;
		return op;
	} else if (fabs(t.times[prntL - nl] - t_next.times[rg_prnt - nl])
			< (t.times[prntL - nl] + t_next.times[rg_prnt - nl]) * TOL) {
		op[0] = iR;
		op[1] = iL;
		return op;
	} else {
		op[0] = iL;
		op[1] = iR;
		return op;
	}
}
