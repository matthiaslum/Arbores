/*
 * commandlineoutput.c
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
#include <stdio.h>
#include "commandlineoutput.h"

void printGeneralInfo(long it, short s, short n_rec, double loglik, double logpost,
		short overwrt) {

	if (overwrt == 1)
		printf("\033[A\033[A\033[A\033[A\033[A\033[A\033[A\033[A\033[A\033[A\033[A");

	printf("Iteration     : %ld\n", it);
	printf("Segment       : %hd\n", s);
	printf("Recombinations: %hd\n", n_rec);
	printf("Log likelihood: %f\n", loglik);
	printf("Log posterior : %f\n\n", logpost);
	printf("*\n");
	printf("*\n");
	printf("*\n");
	printf("*\n");
	printf("*\n");

}

void printCustomMessage(char **message, short lines) {

	// over write old message
	printf("\033[A\033[A\033[A\033[A\033[A");
	printf(
			"                                                                                  \n");
	printf(
			"                                                                                  \n");
	printf(
			"                                                                                  \n");
	printf(
			"                                                                                  \n");
	printf(
			"                                                                                  \n");
	printf("\033[A\033[A\033[A\033[A\033[A");
	for (int i = 0; i < lines; i++)
		printf("%s\n", message[i]);
	for (int i = 0; i < 5 - lines; i++)
		printf("\n");
}

