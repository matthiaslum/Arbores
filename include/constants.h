/*
 * constants.h
 *
 *  Created on: 21.10.2016
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

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

static const int MAX_SEGREGATING_SITES = 100;
static const int MAX_SEQUENCES = 10;
static const int MAX_READ_NODES = 100;
static const int MAX_NODES = 100;
static const int MAX_RECOMBINATIONS = 100;
static const int MAX_REC_ITERATIONS = 10;
static const int MAX_RECREATIONS = 10000;
static const int MAX_ARGS = 1e6;
static const int BR_LEN = 5; // bridge length
static const long MAX_PATHS = 2e6;
static const long MAX_DOMAIN_SIZE = 2e5;
static const long MAX_ADJ_SET_SIZE = 1e5;
static const long MAX_NEW_SITES = 10;
static const long MAX_SAMPLING_SET_SIZE = 1e6;
static const short MAX_STEPS = 2; // max number of recombinations between segregating sites
static const short MAX_PATH_LENGTH = 100;
static const double TOL = 1e-10;
static const double RELTOL = 1e-10;
static const double EPS = 1e-7;
static const char LIKELIHOOD_FILE[] = "/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/like.txt";
static const char DIAGNOSTICS_FILE[] = "/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/results/diagnostics.txt";
static const char CHAIN_FILE[] = "/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/results/chain.tex";


#endif /* CONSTANTS_H_ */
