/*
 * shrub.c
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

#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include "shrub.h"

struct ShrubARG shrub(struct Data data) {

	struct ShrubARG *ARG;
	struct ShrubARG outputARG;
	struct ShrubData sd;
	int n_ARG = 0, min_recs = INT_MAX, argmin = 0, n_leaves = data.n_seq;

	sd = createShrubData(data);
	ARG = malloc(sizeof(struct ShrubARG) * MAX_ARGS);

	printf("Constructing ARG's...           ");

	sd.nni = n_leaves;
	shrub_recursion(sd, ARG, &n_ARG);

	/* Find the minimum recombination ARG */
	printf("\n");
	minARG(ARG, n_ARG, &min_recs, &argmin);
//	printf("\nConstructed %i ARGs\n", n_ARG);

	outputARG.n_seq = n_leaves;
	outputARG.n_nodes = ARG[argmin].n_nodes;
	outputARG.n_rec = ARG[argmin].n_rec;
	outputARG.C = malloc(sizeof(short) * 3 * ARG[argmin].n_nodes);
	shortArrayCopy(outputARG.C, ARG[argmin].C, ARG[argmin].n_nodes, 3, ARG[argmin].n_nodes);

	for(int i = 0; i < n_ARG;i++)
		free(ARG[i].C);
	free(ARG);

	printf("\n");

	deleteShrubData(sd);
	return outputARG;
}

struct ShrubData createShrubData(struct Data data) {
	struct ShrubData sd;

	sd.C = malloc(sizeof(short) * 3 * MAX_NODES);
	sd.M = malloc(sizeof(short) * data.n_seq * data.n_sites);
	sd.rowlab = malloc(sizeof(int) * data.n_seq);
	sd.collab = malloc(sizeof(int) * data.n_sites);
	sd.m = data.n_seq;
	sd.n = data.n_sites;
	sd.nni = 0;
	sd.rec_counts = malloc(sizeof(int) * MAX_REC_ITERATIONS);
	sd.rec_index = 0;

	fillShortArray(sd.C, -1, MAX_NODES, 3, MAX_NODES);
	shortArrayCopy(sd.M, data.M, sd.m, sd.n, sd.m);
	intSimpleSeq(sd.rowlab, sd.m, 0);
	intSimpleSeq(sd.collab, sd.n, 1);

	return sd;
}

void deleteShrubData(struct ShrubData sd) {
	free(sd.C);
	free(sd.M);
	free(sd.collab);
	free(sd.rowlab);
	free(sd.rec_counts);
}

void shrub_recursion(struct ShrubData sd, struct ShrubARG *ARG, int *n_ARG) {

	short *C;
	int i_der;
	struct ShrubData sd1, sd2;

	if (sd.m <= 2) {
		// Terminate recursion
		printf("\b\b\b\b\b\b\b\b\b\b%10i", *n_ARG);

		/* Make a copy of the graph array to be stored */
		C = malloc(sizeof(short) * 3 * sd.nni);
		shortArrayCopy(C, sd.C, sd.nni, 3, MAX_NODES);
		assert(*n_ARG < MAX_ARGS);
		ARG[*n_ARG].C = C;
		ARG[*n_ARG].n_rec = intSum(sd.rec_counts, sd.rec_index);
		ARG[*n_ARG].n_nodes = (short) sd.nni;
		ARG[*n_ARG].n_seq = (short) -1;

		*n_ARG = *n_ARG + 1;
	} else {

		sd1 = mutateAndCoalesce(sd);

		if (sd1.m > 1) {
			// if M-matrix has rows, iterate over all possible derived rows
			for (i_der = 0; i_der < sd1.m; i_der++) {
				if (shortSingleRowSum(sd1.M, i_der, sd1.n, sd1.m) != 0) {
					// don't allow the zero row to be derived
					sd2 = recombine(sd1, i_der);
					shrub_recursion(sd2, ARG, n_ARG);
					deleteShrubData(sd2);
				}
			}
		} else {
			// if the M-matrix has no rows jump to termination
			shrub_recursion(sd1, ARG, n_ARG);
		}
		deleteShrubData(sd1);
	}
}

struct ShrubData recombine(struct ShrubData sd, int i_der) {

	struct ShrubData sd_out;
	short *M, *C;
	int *rl, *cl, *n_recs;
	int n_rows, n_cols, nni, ri;

	int *R;
	short *match_matrix, *alive, *lifespan, *derived;
	int n_R = 0, der_lab, argmax;

	der_lab = sd.rowlab[i_der];

	M = malloc(sizeof(short) * sd.m * sd.n);
	shortArrayCopy(M, sd.M, sd.m, sd.n, sd.m);
	C = malloc(sizeof(short) * 3 * MAX_NODES);
	shortArrayCopy(C, sd.C, MAX_NODES, 3, MAX_NODES);
	rl = malloc(sizeof(int) * sd.m);
	intArrayCopy(rl, sd.rowlab, sd.m, 1, sd.m);
	cl = malloc(sizeof(int) * sd.n);
	intArrayCopy(cl, sd.collab, 1, sd.n, 1);
	n_recs = malloc(sizeof(int) * MAX_REC_ITERATIONS);
	intArrayCopy(n_recs, sd.rec_counts, MAX_REC_ITERATIONS, 1, MAX_REC_ITERATIONS);

	n_cols = sd.n;
	n_rows = sd.m;

	/* Store the derived row */
	derived = malloc(sizeof(short) * sd.n);
	for (int i = 0; i < sd.n; i++)
		derived[i] = sd.M[i_der + i * sd.m];

	/* Remove the derived entry */
	intRemoveCol(rl, i_der, 1, sd.m);
	shortRemoveRow(M, i_der, sd.m, sd.n);
	n_rows--;

	/* Match matrix */
	match_matrix = malloc(sizeof(short) * n_rows * n_cols);
	createMatchMatrix(M, match_matrix, derived, n_rows, n_cols, i_der, sd.m);
	free(derived);

	lifespan = malloc(sizeof(short) * n_rows);
	alive = malloc(sizeof(short) * n_rows);
	for (int i = 0; i < n_rows; i++) {
		alive[i] = match_matrix[i];
		lifespan[i] = match_matrix[i];
	}

	R = malloc(sizeof(int) * 4 * MAX_RECOMBINATIONS);
	fillArray(R, -1, MAX_RECOMBINATIONS, 4, MAX_RECOMBINATIONS);
	for (int j = 1; j < n_cols; j++) {
		lifespanExtension(match_matrix, alive, lifespan, n_rows, j);
		if (shortSum(alive, n_rows) == 0) {
			shortMax(lifespan, n_rows, &argmax);
			updateRecombinations(R, n_R, rl[argmax], argmax, cl[j - 1], j - 1);
			n_R++;
			resetLifeSpan(match_matrix, alive, lifespan, n_rows, j);
		}
	}
	free(match_matrix);
	free(alive);

	/* Final recombination is not really a recombination */
	shortMax(lifespan, n_rows, &argmax);
	updateRecombinations(R, n_R, rl[argmax], argmax, SHRT_MAX, SHRT_MAX);
	n_R++;
	free(lifespan);

	/* Find the index of the derived label and remove from labels*/
	nni = sd.nni;
	ri = sd.rec_index;
	updateGraph(C, R, rl, n_rows, n_R, der_lab, &nni, n_recs, &ri);
	free(R);

	/* create the output*/
	sd_out = shrubDataFromData(C, M, rl, cl, n_rows, n_cols, nni, n_recs, ri, sd.m);
	free(C);
	free(M);
	free(rl);
	free(cl);
	free(n_recs);
	return sd_out;
}

struct ShrubData shrubDataFromData(short *C, short *M, int *rl, int *cl, int m, int n,
		int nni, int *rec_cnts, int rec_index, int m_orig) {

	struct ShrubData sd;

	sd.C = malloc(sizeof(short) * 3 * MAX_NODES);
	shortArrayCopy(sd.C, C, MAX_NODES, 3, MAX_NODES);
	sd.M = malloc(sizeof(short) * m * n);
	shortArrayCopy(sd.M, M, m, n, m_orig);
	sd.rowlab = malloc(sizeof(int) * m);
	intArrayCopy(sd.rowlab, rl, m, 1, m);
	sd.collab = malloc(sizeof(int) * n);
	intArrayCopy(sd.collab, cl, 1, n, 1);
	sd.m = m;
	sd.n = n;
	sd.nni = nni;
	sd.rec_counts = malloc(sizeof(int) * MAX_REC_ITERATIONS);
	intArrayCopy(sd.rec_counts, rec_cnts, MAX_REC_ITERATIONS, 1, MAX_REC_ITERATIONS);
	sd.rec_index = rec_index;

	return sd;
}

struct ShrubData mutateAndCoalesce(struct ShrubData sd) {

	struct ShrubData sd_out;
	short *M, *C;
	int *rl, *cl;
	int n_rows, n_cols, nni, n_rm_col, n_rm_row = 0;

	short *colsums;
	bool *keep;
	bool mutation_fail, coalescence_fail;
	int *i_sort;

	M = malloc(sizeof(short) * sd.m * sd.n);
	shortArrayCopy(M, sd.M, sd.m, sd.n, sd.m);
	C = malloc(sizeof(short) * 3 * MAX_NODES);
	shortArrayCopy(C, sd.C, MAX_NODES, 3, MAX_NODES);
	rl = malloc(sizeof(int) * sd.m);
	intArrayCopy(rl, sd.rowlab, sd.m, 1, sd.m);
	cl = malloc(sizeof(int) * sd.n);
	intArrayCopy(cl, sd.collab, 1, sd.n, 1);
	n_cols = sd.n;
	n_rows = sd.m;
	nni = sd.nni;

	colsums = malloc(sizeof(short) * n_cols);
	keep = malloc(sizeof(bool) * n_cols);
	shortZeros(colsums, 1, n_cols, 1);
	bool cont = true;

	while (cont) {
		n_rm_col = 0;

		/* Mutation i.e. removal of columns */
		shortColumnSum(M, colsums, n_rows, n_cols, sd.m);
		mutation_fail = true;
		for (int i = 0; i < n_cols; i++) {
			keep[i] = colsums[i] > 1;
			mutation_fail = mutation_fail && keep[i];
			n_rm_col += keep[i] ? 0 : 1;
		}
		shortRemoveColumns(M, keep, sd.m, n_cols);
		intRemoveColumns(cl, keep, 1, n_cols);
		n_cols -= n_rm_col; // book keeping for the effective number of columns

		/* Coalescence */
		coalescence_fail = true;
		if (n_rows > 1) {

			i_sort = malloc(sizeof(int) * n_rows);
			shortSortRows(M, i_sort, 0, n_rows, n_cols, sd.m); // just return the sorting index array
			n_rm_row = 0;
			for (int i = 0; i < n_rows - 1; i++) {
				if (equalRows(M, i_sort[i], i_sort[i + 1], n_cols, sd.m)) {

					/* Update tree array */
					C[nni] = rl[i_sort[i]];
					C[nni + MAX_NODES] = rl[i_sort[i + 1]];
					C[nni + 2 * MAX_NODES] = -1;

					rl[i_sort[i]] = nni;

					/* Remove row */
					shortRemoveRow(M, i_sort[i + 1], sd.m, n_cols);
					intRemoveCol(rl, i_sort[i + 1], 1, sd.m);
					nni++;
					n_rm_row++;

					coalescence_fail = false;
					break;
				}
			}
			free(i_sort);
		}
		n_rows -= n_rm_row;
		cont = !(mutation_fail && coalescence_fail);
	}
	free(keep);
	free(colsums);

	sd_out = shrubDataFromData(C, M, rl, cl, n_rows, n_cols, nni, sd.rec_counts,
			sd.rec_index, sd.m);
	free(C);
	free(M);
	free(rl);
	free(cl);

	return sd_out;
}

bool equalRows(short *M, int i1, int i2, int n, int me) {
	bool equal_rows = true;
	for (int j = 0; j < n; j++) {
		equal_rows = equal_rows && (M[i1 + j * me] == M[i2 + j * me]);
		if (!equal_rows)
			break;
	}
	return equal_rows;
}

void lifespanExtension(short *mm, short *alive, short *lifespan, int m, int col) {
	for (int i = 0; i < m; i++) // iterate over the rows
		if ((mm[i + col * m] == 1) && (alive[i] == 1)) {
			alive[i] = 1;
			lifespan[i]++;
		} else {
			alive[i] = 0;
		}
}

void resetLifeSpan(short* mm, short *alive, short *lifespan, int m, int col) {
	for (int i = 0; i < m; i++) {
		alive[i] = mm[i + col * m];
		lifespan[i] = alive[i];
	}
}

void checkCalculations(short *M, int m, int n, int me) {
	short *rowsums;
	/* Double check calculations */
	rowsums = malloc(sizeof(short) * m);
	shortRowSum(M, rowsums, m, n, me);
	for (int i = 0; i < m; i++)
		if (rowsums[i] != 0)
			printf("Row sums should be zero! Something has gone wrong.\n");
	free(rowsums);
}

void createMatchMatrix(short *M, short *mm, short *derived, int m, int n, int i_der,
		int me) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			mm[i + j * m] = derived[j] == M[i + j * me] ? 1 : 0;
}

void updateGraph(short *C, int *R, int *rowlab, int m, int n_R, int der_lab,
		int *next_node_index, int *n_recs, int *ri) {
	int nni = *next_node_index;
	int *tmp_lab;

	/* Copy the row labels */
	tmp_lab = malloc(sizeof(int) * m);
	intArrayCopy(tmp_lab, rowlab, 1, m, 1);

	for (int i = 0; i < n_R; i++) {
		C[nni + MAX_NODES * 0] = der_lab;
		C[nni + MAX_NODES * 1] = (short) tmp_lab[R[i + MAX_RECOMBINATIONS]];
		C[nni + MAX_NODES * 2] = R[i + 2 * MAX_RECOMBINATIONS];
		tmp_lab[R[i + MAX_RECOMBINATIONS]] = nni;
		nni++;
	}

	*next_node_index = nni;

	if (*ri > MAX_REC_ITERATIONS)
		printf("Index exceeding allocation\n");
	n_recs[*ri] = n_R - 1;
	*ri = *ri + 1;

	intArrayCopy(rowlab, tmp_lab, 1, m, 1);
	free(tmp_lab);
}

void updateRecombinations(int *R, int n_R, int v1, int v2, int v3, int v4) {
	R[n_R + MAX_RECOMBINATIONS * 0] = v1;
	R[n_R + MAX_RECOMBINATIONS * 1] = v2;
	R[n_R + MAX_RECOMBINATIONS * 2] = v3;
	R[n_R + MAX_RECOMBINATIONS * 3] = v4;
}

void minARG(struct ShrubARG *ARG, int n_ARGs, int *min_recs, int *argmin) {

	for (int i = 0; i < n_ARGs; i++)
		if (ARG[i].n_rec < *min_recs) {
			*min_recs = ARG[i].n_rec;
			*argmin = i;
		}
}

