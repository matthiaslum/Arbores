/*
 * data.c
 *
 *  Created on: 14.11.2016
 *      Author: heine
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "data.h"
#include "datastructures.h"
#include "utils.h"

// NB: Site must be provided as a site label starting from 1, not 0
short getDataColumn(struct Data data, long site) {

	for (short i = 0; i < data.n_sites; i++)
		if (data.segregating_sites[i] == site)
			return i;

	return -1;
}

void printData(struct Data *data) {

	int *rowlab;
	int max_site = INT_MIN, w;

	rowlab = malloc(sizeof(int) * data->n_seq);
	for (int i = 0; i < data->n_seq; i++) {
		rowlab[i] = i;
		max_site =
				data->segregating_sites[i] > max_site ? data->segregating_sites[i] : max_site;
	}
	w = (int) floor(log10(max_site)) + 2;
	printf("Format width %i\n", w);

	/* Print the data matrix */
	printf("Data: %s\n", data->name);
	printf("Found:\n\t%i sequences\n\t%i sites\n\n", data->n_seq, data->n_sites);
	printf("DNA polymorphisms data:\n");
	printWithLabelsW(data->M, data->segregating_sites, rowlab, w, data->n_sites,
			data->n_seq, data->n_seq);
	free(rowlab);
}


