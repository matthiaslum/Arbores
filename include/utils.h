/*
 *  utils.h
 *  arbores
 *
 *  Created by Kari Heine on 30.9.2016.
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
#ifndef UTILS_H_
#define UTILS_H_

#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <stdbool.h>
#include "datastructures.h"

double findMax(struct DoubleVector vec);
double findMin(struct DoubleVector vec);
int intSum(const int *a, int n);
int shortSum(const short *a, int n);
int *intInsertValue(int *a, int i_insert, int v, int n);
long long nchoosek(int n, int k);
short firstCommon(struct ShortVector vec1, struct ShortVector vec2);
short intIncluded(int v, int *a, long len);
short shortMax(short *a, int len, int *argmax);
short shortSingleRowSum(short *M, int row, int n, int me);
short *flipArray(short *M, int m, int n);
short *shortInsertColumn(short *M, int i_col, short *col, int m, int n);
short *shortChooseColumns(short *M, int *cols, int n_new, int m, int n);
void doubleArrayCopy(double *dst, double *src, int mdst, int n, int msrc);
void fillArray(int *a, int v, int m, int n, int me);
void fillDoubleArray(double *a, double v, int m, int n, int me);
void fillShortArray(short *a, short v, int m, int n, int me);
void intArrayCopy(int *dst, int *src, int mdst, int n, int msrc);
void intInsert(int *a, int v, int ind, int len);
void intRemoveColumns(int *M, bool *keep, int m, int n);
void intRemoveCol(int *M, int col, int m, int n);
void intSimpleSeq(int *a, int len, int offset);
void shortArrayCopy(short *dst, short *src, int mdst, int n, int msrc);
void shortOnes(short *a, int m, int n, int me);
void shortRemoveColumns(short *M, bool *keep, int m, int n);
void shortRemoveRow(short *M, int row, int m, int n);
void shortRemoveCol(short *M, int col, int m, int n);
void shortRowSum(short *M, short *rowsum, int m, int n, int me);
void shortColumnSum(short *M, short *colsum, int m, int n, int me);
void shortZeros(short *a, int m, int n, int me);
void printDoubleArray(const double *M, int m, int n, int me);
void printIntArray(const int *M, int m, int n, int me);
void printIntWithLabels(int *M, int *collab, int *rowlab, int n, int m, int me);
void printIntWithLabelsW(int *M, int *collab, int *rowlab, int w, int n, int m, int me);
void printLongArray(long *M, int m, int n, int me);
void printShortArray(short *M, int m, int n, int me);
void printWithLabels(short *M, int *collab, int *rowlab, int n, int m, int me);
void printWithLabelsW(short *M, int *collab, int *rowlab, int w, int n, int m, int me);

#endif /* UTILS_H_ */
