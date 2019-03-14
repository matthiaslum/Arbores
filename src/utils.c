/*
 * utils.c
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

#include "utils.h"

/**
 Copy 2D array of short integers into another.

 @param dst  Destination array
 @param src  Source array
 @param mdst Number of rows in the destination array
 @param n    Number of columns
 @param msrc Number of rows used for allocating source array
 */
void shortArrayCopy(short *dst, short *src, int mdst, int n, int msrc) {
	for (int i = 0; i < mdst; i++)
		for (int j = 0; j < n; j++)
			dst[i + j * mdst] = src[i + j * msrc];
}

void doubleArrayCopy(double *dst, double *src, int mdst, int n, int msrc) {
	for (int i = 0; i < mdst; i++)
		for (int j = 0; j < n; j++)
			dst[i + j * mdst] = src[i + j * msrc];
}

/**
 Copy 2D array of integers into another.

 @param dst  Destination array
 @param src  Source array
 @param mdst Number of rows
 @param n    Number of columns
 @param msrc Number of rows
 */
void intArrayCopy(int *dst, int *src, int mdst, int n, int msrc) {
	for (int i = 0; i < mdst; i++)
		for (int j = 0; j < n; j++)
			dst[i + j * mdst] = src[i + j * msrc];
}

/**
 Remove row by moving rows with higher index upwards by one element.

 @param M   Modified matrix
 @param row Row to be removed
 @param m   Number of rows
 @param n   Number of columns
 */
void shortRemoveRow(short *M, int row, int m, int n) {
	for (int i = row; i < m - 1; i++)
		for (int j = 0; j < n; j++)
			M[i + j * m] = M[i + 1 + j * m];
}

void shortLongRemoveRow(short *M, long row, long m, long n) {
	for (long i = row; i < m - 1; i++)
		for (long j = 0; j < n; j++)
			M[i + j * m] = M[i + 1 + j * m];
}
/**
 Remove columns from a 2D short integer array according to an boolean indicator array

 @param M    Array from which the columns are deleted
 @param keep indicator array
 @param m    number of columns in M
 @param n    number of rown in M
 */
void shortRemoveColumns(short *M, bool *keep, int m, int n) {
	int shift = 0;
	for (int i = 0; i < n; i++)
		if (!keep[i])
			shortRemoveCol(M, i - shift++, m, n);
}

/**
 Remove columns from a 2D integer array according to an boolean indicator array

 @param M    Array from which the columns are deleted
 @param keep indicator array
 @param m    number of columns in M
 @param n    number of rown in M
 */
void intRemoveColumns(int *M, bool *keep, int m, int n) {
	int shift = 0;
	for (int i = 0; i < n; i++)
		if (!keep[i])
			intRemoveCol(M, i - shift++, m, n);
}

/**
 Remove column by moving columns with higher index leftwards by one element.

 @param M   Modified matrix
 @param col Column to be removed
 @param m   Number of rows
 @param n   Number of columns
 */
void shortRemoveCol(short *M, int col, int m, int n) {
	for (int j = col; j < n - 1; j++)
		for (int i = 0; i < m; i++)
			M[i + j * m] = M[i + (j + 1) * m];
}

short *shortChooseColumns(short *M, int *cols, int n_new, int m, int n) {
	short *out;
	out = malloc(sizeof(short) * m * n_new);
	for (int i = 0; i < n_new; i++)
		for (int j = 0; j < m; j++)
			out[j + i * m] = M[j + cols[i] * m];
	return out;
}

short *flipArray(short *M, int m, int n) {
	short *out;
	out = malloc(sizeof(short) * m * n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			out[j + i * m] = M[j + (n - i - 1) * m];
	return out;
}

short *shortInsertColumn(short *M, int i_col, short *col, int m, int n) {

	short *M_new;
	M_new = malloc(sizeof(short) * m * (n + 1));

	/* copy columns before the insert */
	for (int i = 0; i < m; i++)
		for (int j = 0; j < i_col; j++)
			M_new[i + j * m] = M[i + j * m];

	/* insert */
	for (int i = 0; i < m; i++)
		M_new[i + (i_col) * m] = col[i];

	/* the remaining columns */
	for (int i = 0; i < m; i++)
		for (int j = i_col; j < n; j++)
			M_new[i + (j + 1) * m] = M[i + j * m];

	free(M);
	return M_new;
}

int *intInsertValue(int *a, int i_insert, int v, int n) {

	int *a_new;
	a_new = malloc(sizeof(int) * (n + 1));

	/* copy element before insert */
	for (int i = 0; i < i_insert; i++)
		a_new[i] = a[i];

	/* insert */
	a_new[i_insert] = v;

	/* elements after the insert */
	for (int i = i_insert; i < n; i++)
		a_new[i + 1] = a[i];

	free(a);
	return a_new;
}
/**
 Remove column by moving columns with higher index leftwards by one element.

 @param M   Modified matrix
 @param col Column to be removed
 @param m   Number of rows
 @param n   Number of columns
 */
void intRemoveCol(int *M, int col, int m, int n) {
	for (int j = col; j < n - 1; j++)
		for (int i = 0; i < m; i++)
			M[i + j * m] = M[i + (j + 1) * m];
}

/**
 Column sums of an 2D array of short integers

 @param M      Data array
 @param colsum Array to store the results
 @param m      Number of rows to add up
 @param n      Number of columns
 @param me     Number of rows allocated for M
 */
void shortColumnSum(short *M, short *colsum, int m, int n, int me) {
	for (int j = 0; j < n; j++)
		colsum[j] = 0;

	for (int j = 0; j < n; j++) // iterate columns
		for (int i = 0; i < m; i++) // iterate rows
			colsum[j] += M[i + j * me];
}

/**
 Sum of an array of short integers

 @param a      array
 @param n      number of elements
 */
int shortSum(const short *a, int n) {
	int s = 0;
	for (int j = 0; j < n; j++)
		s += (int) a[j];
	return s;
}

/**
 Sum of an array of integers

 @param a      array
 @param n      number of elements
 */
int intSum(const int *a, int n) {
	int s = 0;
	for (int j = 0; j < n; j++)
		s += (int) a[j];
	return s;
}

/**
 Row sums of an 2D array of short integers

 @param M      Data array
 @param rowsum Array to store the results
 @param m      Number of rows to add up
 @param n      Number of columns
 @param me     Number of rows allocated for M
 */
void shortRowSum(short *M, short *rowsum, int m, int n, int me) {
	for (int j = 0; j < n; j++) // iterate columns
		for (int i = 0; i < m; i++) // iterate rows
			rowsum[i] += M[i + j * me];
}

/**
 Single row sum of an 2D array of short integers

 @param M      Data array
 @param row    Row being summed
 @param n      Number of columns
 @param me     Number of rows allocated for M
 */
short shortSingleRowSum(short *M, int row, int n, int me) {
	short s = 0;
	for (int j = 0; j < n; j++) // iterate columns
		s += M[row + j * me];
	return s;
}

/**
 Print 2D array of short integers.

 @param M  Array
 @param m  Number of displayed rows
 @param n  Number of displayed columns
 @param me Number of actual rows allocated for the array
 */
void printShortArray(short *M, int m, int n, int me) {
	printf("\n");
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%i ", M[i + j * me]);
		printf("\n");
	}
	printf("\n");
}

void printLongArray(long *M, int m, int n, int me) {
	printf("\n");
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%ld ", M[i + j * me]);
		printf("\n");
	}
	printf("\n");
}

void printShortArrayW(short *M, int w, int m, int n, int me) {
	char formstr[10];
	sprintf(formstr, "%%%ii", w);

	printf("\n");
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			printf(formstr, M[i + j * me]);
		}
		printf("\n");
	}
	printf("\n");
}

void printUintArray(unsigned int *M, int m, int n, int me) {
	printf("\n");
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%u ", M[i + j * me]);
		printf("\n");
	}
	printf("\n");
}

/**
 Print 2D array of booleans.

 @param M  Array
 @param m  Number of displayed rows
 @param n  Number of displayed columns
 @param me Number of actual rows allocated for the array
 */
void printBoolArray(bool *M, int m, int n, int me) {
	printf("\n");
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%i ", (int) M[i + j * me]);
		printf("\n");
	}
	printf("\n");
}

/**
 Print 2D array of integers.

 @param M  Array
 @param m  Number of displayed rows
 @param n  Number of displayed columns
 @param me Number of actual rows allocated for the array
 */
void printIntArray(const int *M, int m, int n, int me) {
	printf("\n");
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%i ", M[i + j * me]);
		printf("\n");
	}
	printf("\n");
}

void printDoubleArray(const double *M, int m, int n, int me) {
	printf("\n");
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%8.10f ", M[i + j * me]);
		printf("\n");
	}
	printf("\n");
}

/**
 Fill 2D short integer array with zeros

 @param a  Array
 @param m  Number of filled rows
 @param n  Number of filled columns
 @param me Number of actual rows allocated for the array
 */
void shortZeros(short *a, int m, int n, int me) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			a[i + j * me] = 0;
}

void intZeros(int *a, int m, int n, int me) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			a[i + j * me] = 0;
}

void fillArray(int *a, int v, int m, int n, int me) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			a[i + j * me] = v;
}

void fillShortArray(short *a, short v, int m, int n, int me) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			a[i + j * me] = v;
}

void fillDoubleArray(double *a, double v, int m, int n, int me) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			a[i + j * me] = v;
}

/**
 Fill 2D short integer array with ones

 @param a  Array
 @param m  Number of filled rows
 @param n  Number of filled columns
 @param me Number of actual rows allocated for the array
 */
void shortOnes(short *a, int m, int n, int me) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			a[i + j * me] = 1;
}

/**
 Print array with column and row labels

 @param M      Array to be printed
 @param collab Column labels
 @param rowlab Row labels
 @param n      number of columns
 @param m      number of rows
 @param me     number of rows allocating the array
 */
void printWithLabels(short *M, int *collab, int *rowlab, int n, int m, int me) {
	printf("\n    ");
	for (int j = 0; j < n; j++)
		printf("%3i", collab[j]);
	printf("\n");

	for (int i = 0; i < m; i++) {
		printf("%3i:", rowlab[i]);
		for (int j = 0; j < n; j++)
			printf("%3i", M[i + j * me]);
		printf("\n");
	}
	printf("\n");
}
void printWithLabelsW(short *M, int *collab, int *rowlab, int w, int n, int m, int me) {
	char formstr[10];
	sprintf(formstr, "%%%ii", w);
	printf("\n ");
	for (int i = 0; i < w; i++)
		printf(" ");
	for (int j = 0; j < n; j++)
		printf(formstr, collab[j]);
	printf("\n");

	for (int i = 0; i < m; i++) {
		printf(formstr, rowlab[i]);
		printf(":");
		for (int j = 0; j < n; j++)
			printf(formstr, M[i + j * me]);
		printf("\n");
	}
	printf("\n");
}

/**
 Print array with column and row labels

 @param M      Array to be printed
 @param collab Column labels
 @param rowlab Row labels
 @param n      number of columns
 @param m      number of rows
 @param me     number of rows allocating the array
 */
void printIntWithLabels(int *M, int *collab, int *rowlab, int n, int m, int me) {
	printf("\n    ");
	for (int j = 0; j < n; j++)
		printf("%3i", collab[j]);
	printf("\n");

	for (int i = 0; i < m; i++) {
		printf("%3i:", rowlab[i]);
		for (int j = 0; j < n; j++)
			printf("%3i", M[i + j * me]);
		printf("\n");
	}
	printf("\n");
}

/**
 Print array with column and row labels

 @param M      Array to be printed
 @param collab Column labels
 @param rowlab Row labels
 @param n      number of columns
 @param m      number of rows
 @param me     number of rows allocating the array
 */
void printIntWithLabelsW(int *M, int *collab, int *rowlab, int w, int n, int m, int me) {

	char formstr[10];
	sprintf(formstr, "%%%ii", w);
	printf("\n ");
	for (int j = 0; j < w; j++)
		printf(" ");
	for (int j = 0; j < n; j++)
		printf(formstr, collab[j]);
	printf("\n");

	for (int i = 0; i < m; i++) {
		printf(formstr, rowlab[i]);
		printf(":");
		for (int j = 0; j < n; j++)
			printf(formstr, M[i + j * me]);
		printf("\n");
	}
	printf("\n");
}

/**
 Logical AND with short integers of two arrays of equal length.

 @param a Array 1
 @param b Array 2
 @param c Output
 @param n Common length of arrays a, b, and c
 */
void shortArrayAnd(short *a, short *b, short *c, int n) {
	for (int i = 0; i < n; i++)
		c[i] = a[i] == 1 && b[i] == 1 ? 1 : 0;
}

/**
 Max of a array of short integers.

 @param a   array
 @param len array length

 @return Returns the maximal value in the array.
 */
short shortMax(short *a, int len, int *argmax) {
	short max = SHRT_MIN;

	if (argmax != NULL) // return also the argmin
		for (int i = 0; i < len; i++) {
			if (a[i] > max) {
				max = a[i];
				*argmax = i;
			}
		}
	else
		// return just the maximum
		for (int i = 0; i < len; i++)
			if (a[i] > max)
				max = a[i];

	return max;
}

short shortMin(short *a, int len, int *argmin) {
	short min = SHRT_MAX;

	if (argmin != NULL) // return also the argmin
		for (int i = 0; i < len; i++) {
			if (a[i] < min) {
				min = a[i];
				*argmin = i;
			}
		}
	else
		// return just the maximum
		for (int i = 0; i < len; i++)
			if (a[i] < min)
				min = a[i];

	return min;
}

/**
 Trim rows at the end of matrix of short integers.

 @param A  Matrix to be trimmed
 @param B  Output
 @param m  Rows in output
 @param n  common number of columns
 @param me Number of rows in A
 */
void shortTrimRows(short *A, short *B, int m, int n, int me) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			A[i + j * me] = B[i + j * me];
}

/**
 Trim rows at the end of matrix of integers.

 @param A  Matrix to be trimmed
 @param B  Output
 @param m  Rows in output
 @param n  common number of columns
 @param me Number of rows in A
 */
void intTrimRows(int *A, int *B, int m, int n, int me) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			A[i + j * me] = B[i + j * me];
}

int any(int v, int *a, int len) {
	int out = 0;
	for (int i = 0; i < len; i++) {
		out = a[i] == v ? 1 : 0;
		if (out == 1)
			break;
	}
	return out;
}

void shortSimpleSeq(short *a, short len) {
	for (int i = 0; i < len; i++)
		a[i] = (short) i;
}

void intSimpleSeq(int *a, int len, int offset) {
	for (int i = 0; i < len; i++)
		a[i] = i + offset;
}

void intInsert(int *a, int v, int ind, int len) {
	for (int i = len - 1; i > ind; i--)
		a[i] = a[i - 1];
	a[ind] = v;
}

short intIncluded(int v, int *a, long len) {
	for (long i = 0; i < len; i++)
		if (v == a[i])
			return 1;
	return 0;
}

long long nchoosek(int n, int k) {
	long long ans = 1;
	k = k > n - k ? n - k : k;
	int j = 1;
	for (; j <= k; j++, n--) {
		if (n % j == 0) {
			ans *= n / j;
		} else if (ans % j == 0) {
			ans = ans / j * n;
		} else {
			ans = (ans * n) / j;
		}
	}
	return ans;
}

double findMin(struct DoubleVector vec) {
	double min = DBL_MAX;
	for (int i = 0; i < vec.length; i++)
		min = vec.v[i] < min ? vec.v[i] : min;
	return min;
}

double findMax(struct DoubleVector vec) {
	double max = DBL_MIN;
	for (int i = 0; i < vec.length; i++)
		max = vec.v[i] > max ? vec.v[i] : max;
	return max;
}

short firstCommon(struct ShortVector vec1, struct ShortVector vec2) {

	for (int i = 0; i < vec1.length; i++)
		for (int j = 0; j < vec2.length; j++)
			if (vec1.v[i] == vec2.v[j])
				return vec1.v[i];
	return -1;
}

