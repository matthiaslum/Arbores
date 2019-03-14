/*
 *  fileaccess.h
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

#ifndef fileaccess_h
#define fileaccess_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "datastructures.h"
#include "treeutils.h"
#include "MCMCutils.h"
#include "constants.h"

char separator();
struct Data readDataFile(char *filename);
struct Data readData(char *filename);
struct Smc readInitialisationPathFile(struct Data data);
struct Smc readInitialisationFileRowFormat(struct Data data);
void createInitFilePath(char *initfile);
void createResultFullPahts(char *result_folder);
void parseGraphMatrixEntry(short *C, double *times, const char *sep, short row);
void parseOperationEntry(short *o, const char *sep);
void removeChainFile();
void removeDiagnosticsFile();
void removeMRCAFile();
void writeDiagnosticsFile(struct MCMCDiagnostics dgn);
void writeMrcaToFile(struct MRCA mrca);
void writePath(struct Smc path, FILE *file);
void writePathToChainFile(struct Smc path);
void writePathToInitialisationFile(struct Smc path);
void writeStateToFile(struct MCMCSummary *chain, struct Stats *stats, int len);
void writeTikzFileForMapPath(struct Smc map_path);
void writeTreePathToFile(struct Tree *p, short **ops, short path_len, int *op_sites);

#endif /* fileaccess_h */
