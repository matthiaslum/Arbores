/*
 *  shrub.h
 *  arbores
 *
 *  Created by Kari Heine on 30.9.2016.
 *
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
#ifndef shrub_h
#define shrub_h

#include <stdlib.h>
#include "datastructures.h"
#include "utils.h"
#include "sorting.h"
#include "constants.h"

bool equalRows(short *M, int i1, int i2, int n, int me);
struct ShrubData createShrubData(struct Data data);
struct ShrubData mutateAndCoalesce(struct ShrubData sd);
struct ShrubData shrubDataFromData(short *C, short *M, int *rl, int *cl, int m, int n,
		int nni, int *rec_cnts, int rec_index, int m_orig);
struct ShrubData recombine(struct ShrubData sd, int i_der);
void createMatchMatrix(short *M, short *mm, short *derived, int m, int n, int i_der,
		int me);
void deleteShrubData(struct ShrubData sd);
void lifespanExtension(short *mm, short *alive, short *lifespan, int m, int col);
void minARG(struct ShrubARG *ARG, int n_ARGs, int *min_recs, int *argmin);
void resetLifeSpan(short* mm, short *alive, short *lifespan, int m, int col);
struct ShrubARG shrub(struct Data data);
void shrub_recursion(struct ShrubData sd, struct ShrubARG *ARG, int *n_ARG);
void updateGraph(short *C, int *R, int *rowlab, int m, int n_R, int der_lab,
		int *next_node_index, int *n_recs, int *ri);
void updateRecombinations(int *R, int n_R, int v1, int v2, int v3, int v4);

#endif /* shrub_h */
