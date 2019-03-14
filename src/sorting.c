/*
 * sorting.c
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

#include <stdlib.h>
#include "sorting.h"
#include "utils.h"

static int *intArray;
static short *shortArray;
static double *doubleArray;
static short *array;
static long sort_n_cols;
static long sort_n_rows;

int intComparator(const void *p1, const void *p2) {
	int a = *((int *) p1), b = *((int *) p2);
	return intArray[a] < intArray[b] ? -1 : intArray[a] > intArray[b] ? 1 : 0;
}

int shortComparator(const void *p1, const void *p2) {
	int a = *((int *) p1), b = *((int *) p2);
	return shortArray[a] < shortArray[b] ? -1 : shortArray[a] > shortArray[b] ? 1 : 0;
}

int doubleComparator(const void *p1, const void *p2) {
	int a = *((int *) p1), b = *((int *) p2);
	return doubleArray[a] < doubleArray[b] ? -1 : doubleArray[a] > doubleArray[b] ? 1 : 0;
}

int shortRowComparator(const void *p1, const void *p2) {
	int a = *((int *) p1);
	int b = *((int *) p2);
	int retvalue = 0;
	int cnt = 0;

	while (retvalue == 0 && cnt < sort_n_cols) {
		if (array[a + cnt * sort_n_rows] < array[b + cnt * sort_n_rows])
			retvalue = -1;
		else if (array[a + cnt * sort_n_rows] > array[b + cnt * sort_n_rows])
			retvalue = 1;
		cnt++;
	}
	return retvalue;
}

int longRowComparator(const void *p1, const void *p2) {
	long a = *((long *) p1);
	long b = *((long *) p2);
	int retvalue = 0;
	int cnt = 0;

	while (retvalue == 0 && cnt < sort_n_cols) {
		if (array[a + cnt * sort_n_rows] < array[b + cnt * sort_n_rows])
			retvalue = -1;
		else if (array[a + cnt * sort_n_rows] > array[b + cnt * sort_n_rows])
			retvalue = 1;
		cnt++;
	}
	return retvalue;
}

/**
 Sort the rows in a 2D array of short integers,

 @param M  Array to be sorted
 @param indices  array for returning the sorting indices
 @param apply indicate whether the array M is actually sorted or is only the sorting indices returned.
 @param m  number of rows
 @param n  number of columns
 @param me number of rows used for allocating the array
 */
void shortSortRows(short *M, int *indices, int apply, int m, int n, int me) {

	array = malloc(sizeof(short) * m * n);
	sort_n_cols = n;
	sort_n_rows = m;

	/* Create row index array */
	for (int i = 0; i < m; i++)
		indices[i] = i;

	shortArrayCopy(array, M, m, n, me);

	qsort(indices, m, sizeof(int), shortRowComparator);

	if (apply == 1)
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				M[i + j * me] = array[indices[i] + j * m];

	free(array);
}

void shortLongSortRows(short *M, long *indices, int apply, long m, long n, long me) {

	array = malloc(sizeof(short) * m * n);
	sort_n_cols = n;
	sort_n_rows = m;

	/* Create row index array */
	for (long i = 0; i < m; i++)
		indices[i] = i;

	// copy
	for (long i = 0; i < m; i++)
		for (long j = 0; j < n; j++)
			array[i + j * m] = M[i + j * me];

	qsort(indices, m, sizeof(long), longRowComparator);

	if (apply == 1)
		for (long i = 0; i < m; i++)
			for (long j = 0; j < n; j++)
				M[i + j * me] = array[indices[i] + j * m];

	free(array);
}

void sortIntArray(int *i_sort, int *array, int length) {

	intArray = malloc(sizeof(int) * length);
	for (int i = 0; i < length; i++)
		intArray[i] = array[i];
	qsort(i_sort, length, sizeof(int), intComparator);
	free(intArray);
}

void sortDoubleArray(int *i_sort, double *array, int length, short apply) {

	doubleArray = malloc(sizeof(double) * length);
	for (int i = 0; i < length; i++)
		doubleArray[i] = array[i];
	qsort(i_sort, length, sizeof(int), doubleComparator);

	if(apply==1)
		for(int i = 0; i <length;i++)
			array[i] = doubleArray[i_sort[i]];

	free(doubleArray);
}

