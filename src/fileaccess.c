/*
 * fileaccess.c
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
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "fileaccess.h"
#include "graph2tikz.h"
#include "pathutils.h"
#include "timeadjustment.h"
#include "utils.h"

char diagnostics_full_path[100];
char chain_file[100];
char mrca_file[100];
char initialisation_file[100];
char map_file[100];

char separator() {
#ifdef _WIN32
	return '\\';
#else
	return '/';
#endif
}

struct Data readData(char *filename) {

	struct Data data;
	data = readDataFile(filename);

	if (data.M != NULL)
		printf("Data successfully read.\n");
	else {
		printf("File read failed. Check file name.\n");
		assert(false);
	}

	if(data.n_seq >= 10) {
		printf("WARNING: Depending on the data, the algorithm may be unstable for data containing more than 9 sequences\n");
		printf("Press any key to proceed (or Ctrl+C to quit)\n");
		getchar();
	}

	return data;

}

struct Data readDataFile(char *filename) {

	FILE *file;
	struct Data data;
	char line[10000];
	const char sep[2] = " ";
	char *token;
	short *M;
	int *segr_sites;
	int n_segr_sites = 0, n_seq = 0, row = 0, index = 0, col = 0;
	size_t len;

	// TODO Make a note of the empty spaces at the end of lines

	file = fopen(filename, "r");
	if (file == NULL) {
		data.M = NULL;
		data.n_seq = -1;
		data.n_sites = -1;
		data.name = NULL;
		data.segregating_sites = NULL;
		return data;
	}

	/* Initialise the data structure */
	M = malloc(sizeof(short) * MAX_SEGREGATING_SITES * MAX_SEQUENCES);
	segr_sites = malloc(sizeof(int) * MAX_SEGREGATING_SITES);

	/* Read the file */
	while (fgets(line, sizeof line, file) != NULL) {
		/* Read the first line */
		if (row == 0) {
			token = strtok(line, sep);
			while (token != NULL) {
				if (strcmp(token, sep) != 0)
					segr_sites[index++] = (short) atoi(token);
				token = strtok(NULL, sep);
			}
			index = 0;
		}
		/* Read the actual data lines */
		else {
			token = strtok(line, sep);
			col = 0;
			while (token != NULL) {
				if (strcmp(token, sep) != 0)
					M[row - 1 + MAX_SEQUENCES * col++] = (short) atoi(token);
				token = strtok(NULL, sep);
			}
		}
		row++;
	}

	n_seq = row - 1;
	n_segr_sites = col;

	/* Create the data structure */

	data.M = malloc(sizeof(short) * n_seq * n_segr_sites);
	for (int i = 0; i < n_seq; i++)
		for (int j = 0; j < n_segr_sites; j++)
			data.M[i + j * n_seq] = M[i + j * MAX_SEQUENCES];
	free(M); // free the temporary M matrix

	for (int i = 0; i < n_segr_sites; i++)
		if (data.M[i * n_seq] != 0) {
			printf("The first row of the data must contain zeros only.\n");
			assert(false);
		}

	data.segregating_sites = malloc(sizeof(int) * n_segr_sites);
	for (int j = 0; j < n_segr_sites; j++)
		data.segregating_sites[j] = segr_sites[j];
	free(segr_sites); // free the temporary M matrix

	data.n_seq = n_seq;
	data.n_sites = n_segr_sites;

	len = strlen(filename);
	data.name = malloc(sizeof(char) * len);
	strncpy(data.name, filename, len);

	fclose(file);

	return data;
}

void writeTreePathToFile(struct Tree *p, short **ops, short path_len, int *op_sites) {

	FILE *file;
	short n_b, n_l;
	n_l = (p[0].n_nodes + 1) / 2;
	n_b = n_l - 1;

	file = fopen("initialisation.txt", "w");
	for (int i = 0; i < path_len; i++) {
		for (int j = 0; j < n_b; j++)
			fprintf(file, "%i %hd %hd %10.4f\n", j + n_l, p[i].C[j], p[i].C[j + n_b],
					p[i].times[j]);
		fprintf(file, "SITE %i\n", op_sites[i]);
		if (i < path_len - 1)
			fprintf(file, "SPR %hd %hd\n", ops[i][0], ops[i][1]);
	}
	fclose(file);
}

void writeStateToFile(struct MCMCSummary *chain, struct Stats *stats, int len) {

	FILE *file;
	short n_n, n_l, n_b, nops, path_len;

	file = fopen("output.txt", "a");

	if (len <= 0) {
		printf("Nothing to print. Array is empty.\n");
		return;
	}

	// these are constant over all trees and iterations
	n_n = chain[0].path.tree_path[0].n_nodes;
	n_l = (n_n + 1) / 2;
	n_b = n_l - 1;

	for (int k = 0; k < len; k++) {

		path_len = chain[k].path.path_len;
		nops = path_len - 1;

		fprintf(file, "%hd,%hd,", n_l, nops);

		/* write graph matrices */
		for (int j = 0; j < path_len; j++) {
			for (int i = 0; i < n_b; i++)
				fprintf(file, "%hd,", chain[k].path.tree_path[j].C[i]);
			for (int i = 0; i < n_b; i++)
				fprintf(file, "%hd,", chain[k].path.tree_path[j].C[i + n_b]);
		}

		/* write times */
		for (int j = 0; j < path_len; j++)
			for (int i = 0; i < n_b; i++)
				fprintf(file, "%e,", chain[k].path.tree_path[j].times[i]);

		/* write operations */
		for (int j = 0; j < nops; j++)
			fprintf(file, "%hd,%hd,", chain[k].path.opers[j][0], chain[k].path.opers[j][1]);

		/* write recombination sites */
		for (int j = 0; j < nops; j++)
			fprintf(file, "%i,", chain[k].path.sites[j]);

		/* write recombination times */
		for (int j = 0; j < nops; j++)
			fprintf(file, "%e", chain[k].path.rec_times[j]);
		fprintf(file, "\n");
	}

	fclose(file);
}

struct Smc readInitialisationFileRowFormat(struct Data data) {

	FILE *file;
	struct Smc path;
	char line[10000];
	char *token;
	const char sep[2] = " ";
	short nl, nb, cnt = 0, n_global = 0;
	struct ShortVector selector;

	file = fopen(initialisation_file, "r");
	if (!file) {
		printf("No initialisation found.\n");
		path.path_len = 0;
		assert(false);
		return path;
	}

	/* All data is on one line */
	fgets(line, sizeof(line), file);

	/* Read the number of leaves */
	token = strtok(line, sep);
	nl = (short) atoi(token);
	nb = nl - 1;
	if (nl <= 0) {
		printf("Invalid initialisation file.\n");
		path.path_len = 0;
		return path;
	}

	path = createPath(MAX_PATH_LENGTH);

	while (true) {
		path.tree_path[cnt] = createTree(nl);

		/* graph matrix */
		for (int i = 0; i < nb; i++) {
			token = strtok(NULL, sep);
			path.tree_path[cnt].C[i] = (short) atoi(token);
		}
		for (int i = 0; i < nb; i++) {
			token = strtok(NULL, sep);
			path.tree_path[cnt].C[i + nb] = (short) atoi(token);
		}

		/* times */
		for (int i = 0; i < nb; i++) {
			token = strtok(NULL, sep);
			path.tree_path[cnt].times[i] = (double) atof(token);
		}

		/* tree site */
		token = strtok(NULL, sep);
		path.sites[cnt] = (int) atoi(token);

		token = strtok(NULL, sep);
		if (atoi(token) == -1) {
			cnt++;
			// last tree with identity operation --> terminate parsing loop
			break;
		}

		/* operation */
		path.opers[cnt] = malloc(sizeof(short) * 2);
		path.opers[cnt][0] = (short) atoi(token);
		token = strtok(NULL, sep);
		path.opers[cnt][1] = (short) atoi(token);
		token = strtok(NULL, sep);

		/* recombination time */
		path.rec_times[cnt] = (double) atof(token);

		token = strtok(NULL, sep); // number of leaves in the next tree, which can be ignored

		cnt++;
	}

	/* trim preallocated arrays to right size */
	path.tree_path = realloc(path.tree_path, sizeof(struct Tree) * cnt);
	path.opers = realloc(path.opers, sizeof(short *) * (cnt - 1));
	path.sites = realloc(path.sites, sizeof(int) * cnt);
	path.rec_times = realloc(path.rec_times, sizeof(double) * (cnt - 1));
	free(path.is_free);
	path.is_free = malloc(sizeof(short) * cnt);
	fillShortArray(path.is_free, -1, 1, cnt, 1);

	/* global indexing */
	free(path.global_index);
	path.global_index = malloc(sizeof(short *) * cnt);
	for (int i = 0; i < cnt; i++)
		assignGlobalIndices(path.global_index, &n_global, i, path.tree_path);

	path.path_len = cnt;

	selector = createTreeSelector(path, data);
	path.tree_selector = selector.v;
	path.selector_length = (int) selector.length;

	return path;

}

struct Smc readInitialisationPathFile(struct Data data) {

	const int MAX_PATH_LEN = 100;

	struct ShortVector selector;
	struct Smc path;
	struct Tree *p, *tree_path;
	short **ops, **opers;
	double *times;
	short *C;
	short n_b = -1, row = 0, n_trees = 0;
	const char sep[2] = " ";
	char line[2 * MAX_SEGREGATING_SITES];
	char *token;
	int *op_sites, *sites;
	FILE *file;
	file = fopen("initialisation.txt", "r");

	if (!file) {
		printf("No precalculated initialisation found.\n");
		path.path_len = 0;
		return path;
	}

	p = malloc(sizeof(struct Tree) * MAX_PATH_LEN);
	ops = malloc(sizeof(short *) * MAX_PATH_LEN);
	op_sites = malloc(sizeof(int) * MAX_PATH_LEN);
	C = malloc(sizeof(short) * 2 * MAX_READ_NODES);
	times = malloc(sizeof(double) * MAX_READ_NODES);

	while (fgets(line, sizeof line, file) != NULL) {
		token = strtok(line, sep);
		if (strcmp(token, "SPR") == 0) {
			n_b = row;
			p[n_trees] = treeFromData(C, times, 2 * (n_b + 1) - 1, MAX_READ_NODES);
			ops[n_trees] = malloc(sizeof(short) * 2);
			parseOperationEntry(ops[n_trees], sep);
			row = 0;
			n_trees++;
		} else if (strcmp(token, "SITE") == 0) {
			token = strtok(NULL, sep);
			op_sites[n_trees] = (int) atoi(token);
		} else {
			parseGraphMatrixEntry(C, times, sep, row);
			row++;
		}
	}
	fclose(file);
	if (n_b < 0) {
		free(p);
		free(ops);
		free(op_sites);
		free(C);
		free(times);
		path.path_len = 0;
		return path;
	}

	/* path should have one less operations than trees */
	p[n_trees] = treeFromData(C, times, 2 * (n_b + 1) - 1, MAX_READ_NODES);
	n_trees++;

	free(times);
	free(C);

	/* Trim the path data to save memory */
	tree_path = malloc(sizeof(struct Tree) * n_trees);
	for (int i = 0; i < n_trees; i++) {
		tree_path[i] = duplicateTree(p[i]);
		deleteTree(p[i]);
	}
	opers = malloc(sizeof(short*) * n_trees);
	for (int i = 0; i < n_trees - 1; i++) {
		opers[i] = malloc(sizeof(short) * 2);
		opers[i][0] = ops[i][0];
		opers[i][1] = ops[i][1];
		free(ops[i]);
	}
	sites = malloc(sizeof(int) * n_trees);
	for (int i = 0; i < n_trees; i++)
		sites[i] = op_sites[i];

	path.opers = opers;
	path.path_len = n_trees;
	path.rec_times = malloc(sizeof(double) * n_trees);
	path.sites = sites;
	path.tree_path = tree_path;
	selector = createTreeSelector(path, data);
	path.tree_selector = selector.v;
	path.selector_length = (int) selector.length;

	free(op_sites);
	free(p);
	free(ops);

	return path;
}

void parseGraphMatrixEntry(short *C, double *times, const char *sep, short row) {
	char *token;
	token = strtok(NULL, sep);
	C[row] = (short) atoi(token);
	token = strtok(NULL, sep);
	C[row + MAX_READ_NODES] = (short) atoi(token);
	token = strtok(NULL, sep);
	times[row] = atof(token);
}

void parseOperationEntry(short *o, const char *sep) {
	char *token;
	token = strtok(NULL, sep);
	o[0] = (short) atoi(token);
	token = strtok(NULL, sep);
	o[1] = (short) atoi(token);
}

void makeGraphFile(int *C, int n_nodes, int n_leaves, int *parent_cnts) {
	FILE *file;

	file = fopen("graph.txt", "w");
	fprintf(file, "digraph mygraph{\n\t{rank=same;");
	/* Ensure that leaves are on the same row in the picture */
	for (int i = 0; i < n_leaves; i++)
		fprintf(file, "%i ", i);
	fprintf(file, "}\n");

	for (int i = 0; i < n_nodes; i++) {
		if (i < n_leaves) {
			// Leaf node formatting
			fprintf(file, "%i [label=\"%i\",", i, i);
			fprintf(file, "fillcolor=black, fontcolor=white, shape=circle, style=filled];\n");
		} else if (parent_cnts[i] <= 1) {
			// Normal coalescence node formatting
			fprintf(file, "%i [label=\"%i\",", i, i);
			fprintf(file, "fillcolor=white, fontcolor=black, shape=circle, style=filled];\n");
		} else {
			// Recombination node formatting
			fprintf(file, "%i [label=\"%i\",", i, i);
			fprintf(file,
					"fillcolor=white, fontcolor=black, shape=rectangle, style=filled];\n");
		}
	}
	for (int i = n_leaves; i < n_nodes; i++) {
		fprintf(file, "%i -> %i; \n", i, C[i]);
		fprintf(file, "%i -> %i; \n", i, C[i + n_nodes]);
	}
	fprintf(file, "}\n");
	fclose(file);
}

void createResultFullPahts(char *result_folder) {

	char sep[2];
	sep[0] = separator();
	sep[1] = '\0';

	diagnostics_full_path[0] = '.';
	diagnostics_full_path[1] = separator();
	diagnostics_full_path[2] = '\0';

	strcat(diagnostics_full_path, result_folder);
	strcat(diagnostics_full_path, sep);
	strcat(diagnostics_full_path, "diagnostics.txt");

	mrca_file[0] = '.';
	mrca_file[1] = separator();
	mrca_file[2] = '\0';

	strcat(mrca_file, result_folder);
	strcat(mrca_file, sep);
	strcat(mrca_file, "mrca.txt");

	chain_file[0] = '.';
	chain_file[1] = separator();
	chain_file[2] = '\0';

	strcat(chain_file, result_folder);
	strcat(chain_file, sep);
	strcat(chain_file, "chain.txt");

	map_file[0] = '.';
	map_file[1] = separator();
	map_file[2] = '\0';

	strcat(map_file, result_folder);
	strcat(map_file, sep);
	strcat(map_file, "map.tex");

}

void createInitFilePath(char *initfile) {

	char sep[2];
	sep[0] = separator();
	sep[1] = '\0';

	initialisation_file[0] = '.';
	initialisation_file[1] = separator();
	initialisation_file[2] = '\0';

	strcat(initialisation_file, initfile);

}

void removeChainFile() {
	remove(chain_file);
}

void removeMRCAFile() {
	remove(mrca_file);
}

void removeDiagnosticsFile() {
	remove(diagnostics_full_path);
}

void writePath(struct Smc path, FILE *file) {

	short nl = (path.tree_path[0].n_nodes + 1) / 2, nb = nl - 1;

	for (int i = 0; i < path.path_len; i++) {
		fprintf(file, "%hd ", nl);
		for (int j = 0; j < nb; j++)
			fprintf(file, "%hd ", path.tree_path[i].C[j]);
		for (int j = 0; j < nb; j++)
			fprintf(file, "%hd ", path.tree_path[i].C[j + nb]);
		for (int j = 0; j < nb; j++)
			fprintf(file, "%15.10e ", path.tree_path[i].times[j]);
		fprintf(file, "%i ", path.sites[i]);
		if (i == path.path_len - 1) {
			fprintf(file, "%i %i ", -1, -1);
			fprintf(file, "%15.10e ", (double) -1);
		} else {
			fprintf(file, "%hd %hd ", path.opers[i][0], path.opers[i][1]);
			fprintf(file, "%15.10e ", path.rec_times[i]);
		}
	}
	fprintf(file, "\n");
}

void writePathToChainFile(struct Smc path) {

	FILE *file;
	file = fopen(chain_file, "a");
	writePath(path, file);
	fclose(file);
}

void writePathToInitialisationFile(struct Smc path) {

	FILE *file;

	remove(initialisation_file);

	file = fopen(initialisation_file, "w");
	writePath(path, file);
	fclose(file);

}

void writeMrcaToFile(struct MRCA mrca) {
	FILE *file;
	file = fopen(mrca_file, "a");

	for (int i = 0; i < mrca.nl - 1; i++)
		for (int j = 0; j < mrca.nl - 1 - i; j++)
			fprintf(file, "%15.10e ", mrca.times[i][j]);
	fprintf(file, "\n");
	fclose(file);
}

void writeTikzFileForMapPath(struct Smc map_path) {

	writeTikzTexFileForTreePath(map_file, map_path.tree_path, map_path.opers, NULL, NULL,
			map_path.global_index, map_path.sites, map_path.rec_times, map_path.path_len);

}

void writeDiagnosticsFile(struct MCMCDiagnostics dgn) {

    FILE *file = fopen(diagnostics_full_path, "a");
    fprintf(file, "%15.10e ", dgn.alpha); // 1
    fprintf(file, "%15.10e ", dgn.current_free_time_density);
    fprintf(file, "%15.10e ", dgn.proposed_free_time_density);
    fprintf(file, "%15.10e ", dgn.current_log_likelihood);
    fprintf(file, "%15.10e ", dgn.proposed_log_likelihood); // 5
    fprintf(file, "%15.10e ", dgn.current_log_prior);
    fprintf(file, "%15.10e ", dgn.proposed_log_prior);
    fprintf(file, "%15hd ", dgn.proposed_number_of_free_times);
    fprintf(file, "%15hd ", dgn.current_number_of_free_times);
    fprintf(file, "%15hd ", dgn.proposed_number_of_recombinations); //10
    fprintf(file, "%15hd ", dgn.current_number_of_recombinations);
    fprintf(file, "%15hd ", dgn.accept_indicator);
    fprintf(file, "%15hd ", dgn.jitter_step);
    fprintf(file, "%15.10e ", dgn.log_likelihood); // 14
    fprintf(file, "%15.10e ", dgn.log_prior); // 15
    fprintf(file, "%15.10e ", dgn.log_posterior); // 16
    fprintf(file, "%15.10e ", dgn.cardinality_ratio); // 17
    fprintf(file, "%15.10e ", dgn.proposed_recombination_density); // 18
    fprintf(file, "%15.10e ", dgn.current_recombination_density); // 19
    fprintf(file, "%15hd ", dgn.irreducibility); // 20
    fprintf(file, "%15.10e ", dgn.u); //21
    fprintf(file, "%s ", dgn.indicators); //22
    fprintf(file, "\n");
    fclose(file);
}