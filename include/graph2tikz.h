/*
 * graph2tikz.h
 *
 *  Created on: 17.10.2016
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

#ifndef GRAPH2TIKZ_H_
#define GRAPH2TIKZ_H_

#include <stdio.h>

void writeTikzTexFileForTreePath(char *filename, struct Tree *t, short **o,
		enum Node_Colors **colors, short *new_nodes, short **glb_idx, int *sites,
		double *rec_times, short path_len);
void writeContent(FILE *file, struct Tree *t, short **o, enum Node_Colors **colors,
		short* new_nodes, short **glb_idx, short path_len, char *rel_path, int *sites,
		double *rec_times);
void writeTreeParameters(struct Tree t, short *o, enum Node_Colors *colors,
		short *glb_idx, short new_node, FILE *file, short nb, short n_l, int index,
		double maxt, long site, double regraft_time, double rec_time);
void appendChainTikzFile(struct Smc path, struct Data data, short first);
void closeChainTikzFile();
void initChainTikzFile(struct Smc path, struct Data data);
void printLatexPreamble(FILE *file);

#endif /* GRAPH2TIKZ_H_ */
