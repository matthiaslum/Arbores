/*
 * graph2tikz.c
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
#include <float.h>
#include <stdio.h>
#include "constants.h"
#include "data.h"
#include "datastructures.h"
#include "graph2tikz.h"
#include "pathutils.h"
#include "treeutils.h"
#include "utils.h"

void initChainTikzFile(struct Smc path, struct Data data) {

	FILE *file;

	file = fopen(CHAIN_FILE, "w");
	printLatexPreamble(file);
	fprintf(file, "\\begin{tabular}{l}\n");
	fclose(file);

	appendChainTikzFile(path, data, 1);

}

void appendChainTikzFile(struct Smc path, struct Data data, short first) {

	FILE *file;
	enum Node_Colors **colors;
	short col;

	file = fopen(CHAIN_FILE, "a");

	if (first == 0)
		fprintf(file, "\n\\\\\n\n");

	fprintf(file, "\\begin{tikzpicture}\n");
	fprintf(file, "\\def\\hstep{.5pt}\n\\def\\vstep{.25pt}\n");

	colors = malloc(sizeof(enum Node_Colors*) * path.path_len);
	for (int i = 0; i < path.path_len; i++) {
		col = getDataColumn(data, path.sites[i]);
		if (col >= 0) {
			colors[i] = colorTree(path.tree_path[i], data.M, col);
		} else {
			colors[i] = malloc(sizeof(enum Node_Colors) * (2 * data.n_seq - 1));
			for (int j = 0; j < 2 * data.n_seq - 1; j++)
				colors[i][j] = white;
		}
	}

	writeContent(file, path.tree_path, path.opers, colors, NULL, path.global_index,
			path.path_len, "../", path.sites, path.rec_times);

	fprintf(file, "\\end{tikzpicture}\n");

	for (int i = 0; i < path.path_len; i++)
		free(colors[i]);
	free(colors);
	fclose(file);
}

void closeChainTikzFile() {
	FILE *file;
	file = fopen(CHAIN_FILE, "a");
	fprintf(file, "\\end{tabular}\n");
	fprintf(file, "\\end{document}\n");
	fclose(file);
}

void writeTikzTexFileForTreePathSet(struct Smc *path, long n_paths, struct Data data) {

	FILE *file;
	enum Node_Colors **colors;
	short col;

	if (n_paths < 1)
		return;

	file = fopen("/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/tikzgraph.tex",
			"w");

	printLatexPreamble(file);

	fprintf(file, "\\begin{tabular}{l}\n");

	for (int i = 0; i < n_paths; i++) {

		if (path[i].path_len < 1)
			continue;

		colors = malloc(sizeof(enum Node_Colors *) * path[i].path_len);
		for (int j = 0; j < path[i].path_len; j++) {
			col = getDataColumn(data, path[i].sites[j]);
			printf("%i (%hd), ", path[i].sites[j], col);
			colors[j] = colorTree(path[i].tree_path[j], data.M, col);
		}
		printf("\n");

		fprintf(file, "\\begin{minipage}{\\textwidth}\n");
		fprintf(file, "\\begin{tikzpicture}\n");
		fprintf(file, "\\def\\hstep{.5pt}\n\\def\\vstep{.25pt}\n");

		writeContent(file, path[i].tree_path, path[i].opers, colors, NULL,
				path[i].global_index, path[i].path_len, "./", path[i].sites, NULL);

		for (int j = 0; j < path[i].path_len; j++)
			free(colors[j]);
		free(colors);

		if (i == n_paths - 1)
			fprintf(file, "\\end{tikzpicture}\n\\end{minipage}\n");
		else
			fprintf(file, "\\end{tikzpicture}\n\\end{minipage}\\\\\n");

	}

	fprintf(file, "\\end{tabular}\n");

	fprintf(file, "\\end{document}\n");
	fclose(file);
}

void printLatexPreamble(FILE *file) {

	fprintf(file, "\\documentclass{standalone}\n\n");
	fprintf(file,
			"\\usepackage{tikz}\n\\usetikzlibrary{shapes,arrows,positioning,shapes.geometric}\n");
	fprintf(file, "\\usepackage{ifthen}\n");
	fprintf(file, "\n");
	fprintf(file, "\\setlength{\\textwidth}{40cm}\n");
	fprintf(file, "\n");
	fprintf(file, "\\begin{document}\n");

}

void writeTikzTexFileForTreePath(char *filename, struct Tree *t, short **o,
		enum Node_Colors **colors, short *new_nodes, short **glb_idx, int *sites,
		double *rec_times, short path_len) {

	FILE *file;
	short sites_allocated = 0;

	if (filename == NULL)
		file = fopen("/Users/heine/Dropbox/Work/workspaceMars/Arbores/Debug/tikzgraph.tex",
				"w");
	else
		file = fopen(filename, "w");

	if (sites == NULL) {
		sites = malloc(sizeof(int) * path_len);
		for (int i = 0; i < path_len; i++)
			sites[i] = 0;
		sites_allocated = 1;
	}
	printLatexPreamble(file);

// Begin the tikz-picture itself
	fprintf(file, "\\begin{tikzpicture}\n");
	fprintf(file, "\\def\\hstep{.5pt}\n\\def\\vstep{.25pt}\n");
	writeContent(file, t, o, colors, new_nodes, glb_idx, path_len, "./", sites, rec_times);
	fprintf(file, "\\end{tikzpicture}\n");
	fprintf(file, "\\end{document}\n");

	fclose(file);

	if (sites_allocated == 1)
		free(sites);
}

void writeContent(FILE *file, struct Tree *t, short **o, enum Node_Colors **colors,
		short* new_nodes, short **glb_idx, short path_len, char *rel_path, int *sites,
		double *rec_times) {

	short nn = t[0].n_nodes;
	short n_l = (nn + 1) / 2, n_b = n_l - 1, glb_idx_allocated = 0, rec_times_allocated = 0;
	double maxt = -DBL_MAX;
	double *regraft_times;

	if (glb_idx == NULL) {
		glb_idx = malloc(sizeof(short *) * path_len);
		glb_idx_allocated = 1;
		for (int i = 0; i < path_len; i++)
			glb_idx[i] = NULL;
	}

	if (o != NULL)
		regraft_times = calculateRegraftTimes(t, path_len);
	else {
		regraft_times = malloc(sizeof(double) * path_len);
		fillDoubleArray(regraft_times, -1.0, 1, path_len, 1);
	}

	if (rec_times == NULL) {
		rec_times_allocated = 1;
		rec_times = malloc(sizeof(double) * (path_len - 1));
		fillDoubleArray(rec_times, -1.0, 1, path_len - 1, 1);
	}

	for (int i = 0; i < path_len; i++)
		for (int j = 0; j < n_b; j++)
			maxt = t[i].times[j] > maxt ? t[i].times[j] : maxt;

	if (colors != NULL) {
		for (int i = 0; i < path_len; i++) {
			nn = (new_nodes == NULL) ? -1 : new_nodes[i];
			if (i < path_len - 1 && o != NULL)
				writeTreeParameters(t[i], o[i], colors[i], glb_idx[i], nn, file, n_b, n_l, i,
						maxt, (long) sites[i], regraft_times[i], rec_times[i]);
			else
				writeTreeParameters(t[i], NULL, colors[i], glb_idx[i], nn, file, n_b, n_l, i,
						maxt, (long) sites[i], regraft_times[i], -1.0);
			fprintf(file, "\\input{%sphytree.tex}\n\n", rel_path);
		}
	} else {
		for (int i = 0; i < path_len; i++) {
			nn = (new_nodes == NULL) ? -1 : new_nodes[i];
			if (i < path_len - 1 && o != NULL)
				writeTreeParameters(t[i], o[i], NULL, glb_idx[i], nn, file, n_b, n_l, i, maxt,
						(long) sites[i], regraft_times[i], rec_times[i]);
			else
				writeTreeParameters(t[i], NULL, NULL, glb_idx[i], nn, file, n_b, n_l, i, maxt,
						(long) sites[i], regraft_times[i], -1.0);
			fprintf(file, "\\input{%sphytree.tex}\n\n", rel_path);
		}
	}
	if (glb_idx_allocated == 1)
		free(glb_idx);
	free(regraft_times);
	if (rec_times_allocated == 1)
		free(rec_times);
}

void writeTreeParameters(struct Tree t, short *o, enum Node_Colors *colors,
		short *glb_idx, short new_node, FILE *file, short nb, short n_l, int index,
		double maxt, long site, double regraft_time, double rec_time) {

	short pr_parent = -1, rg_parent = -1;
	double *x;
	double tree_dist = (double) ((n_l + 1) * index), minx = DBL_MAX, maxx = -DBL_MAX;

// Write first columns of the canonical matrix
	fprintf(file, "\\def\\pone{{");
	for (int i = 0; i < nb - 1; i++)
		fprintf(file, "%hd,", t.C[i]);
	fprintf(file, "%hd}}\n", t.C[nb - 1]);
// Write second column of the canonical matrix
	fprintf(file, "\\def\\ptwo{{");
	for (int i = 0; i < nb - 1; i++)
		fprintf(file, "%hd,", t.C[i + nb]);
	fprintf(file, "%hd}}\n", t.C[2 * nb - 1]);

	fprintf(file, "\\def\\site{%ld}\n", site);

// Write node x-coordinate data
	x = malloc(sizeof(double) * t.n_nodes);
	fprintf(file, "\\def\\X{{");
	leafXCoordinates(t, x);

	for (int i = 0; i < n_l; i++) {
		x[i] += tree_dist;
		fprintf(file, "%8.4f,", x[i]);
	}
	for (int i = 0; i < nb - 1; i++) {
		if (t.C[i] >= 0 && t.C[i + nb] >= 0)
			x[i + n_l] = (x[t.C[i]] + x[t.C[i + nb]]) / (double) 2;
		else if (t.C[i] >= 0)
			x[i + n_l] = x[t.C[i]];
		else
			x[i + n_l] = x[t.C[i + nb]];
		fprintf(file, "%8.4f,", x[i + n_l]);
	}
	if (t.C[nb - 1] >= 0 && t.C[nb - 1 + nb] >= 0)
		x[nb - 1 + n_l] = (x[t.C[nb - 1]] + x[t.C[nb - 1 + nb]]) / (double) 2;
	else if (t.C[nb - 1] >= 0)
		x[nb - 1 + n_l] = x[t.C[nb - 1]];
	else
		x[nb - 1 + n_l] = x[t.C[nb - 1 + nb]];
	fprintf(file, "%8.4f}}\n", x[nb - 1 + n_l]);

// minimum and maximum of node x-coordinates for picture formatting purposes
	for (int i = 0; i < n_l; i++) {
		minx = x[i] < minx ? x[i] : minx;
		maxx = x[i] > maxx ? x[i] : maxx;
	}

	fprintf(file, "\\def\\minx{%8.4f}\n\\def\\maxx{%8.4f}\n", minx, maxx);
	free(x);

// write node y-coordinate data
	fprintf(file, "\\def\\T{{");
	for (int i = 0; i < n_l; i++)
		fprintf(file, "%8.4f,", (double) 0);
	for (int i = 0; i < nb - 1; i++)
		fprintf(file, "%8.4f,", t.times[i]);
	fprintf(file, "%8.4f}}\n", t.times[nb - 1]);

// write normalised node y-coordinate data
	fprintf(file, "\\def\\Tnorm{{");
	for (int i = 0; i < n_l; i++)
		fprintf(file, "%8.4f,", (double) 0);
	for (int i = 0; i < nb - 1; i++)
		fprintf(file, "%8.4f,", t.times[i] / maxt * 20);
	fprintf(file, "%8.4f}}\n", t.times[nb - 1] / maxt * 20);

// node labels
	fprintf(file, "\\def\\LAB{{");
	for (int i = 0; i < t.n_nodes - 1; i++)
		fprintf(file, "%i,", i);
	fprintf(file, "%i}}\n", t.n_nodes - 1);

// Node colours
	if (colors != NULL) {
		fprintf(file, "\\def\\colors{{");
		for (int i = 0; i < t.n_nodes - 1; i++)
			fprintf(file, "%u,", colors[i]);
		fprintf(file, "%u}}\n", colors[t.n_nodes - 1]);
		fprintf(file, "\\def\\coloring{1}\n");
	} else {
		fprintf(file, "\\def\\coloring{-1}\n");
	}

// Global indices
	if (glb_idx != NULL) {
		fprintf(file, "\\def\\globalindex{{");
		for (int i = 0; i < nb - 1; i++)
			fprintf(file, "%hd,", glb_idx[i]);
		fprintf(file, "%hd}}\n", glb_idx[nb - 1]);
		fprintf(file, "\\def\\globalindiceson{1}\n");
	} else {
		fprintf(file, "\\def\\globalindiceson{0}\n");
	}

// regraft time
	fprintf(file, "\\def\\regrafttime{%8.4f}\n\\def\\regrafttimenorm{%8.4f}\n",
			regraft_time, regraft_time / maxt * 20);

	// recombination time
	fprintf(file, "\\def\\recombinationtime{%8.4f}\n\\def\\recombinationtimenorm{%8.4f}\n",
			rec_time, rec_time / maxt * 20);

// node numbers
	fprintf(file, "\\def\\NL{%hd}\n\\def\\NN{%hd}\n", n_l, t.n_nodes);

// prune and regraft nodes and their parents
	if (o == NULL) {
		fprintf(file, "\\def\\prune{%hd}\n\\def\\regraft{%hd}\n\n", (short) -1, (short) -1);
		fprintf(file, "\\def\\prpar{%hd}\n\\def\\rgpar{%hd}\n\n", (short) -1, (short) -1);
	} else {
		fprintf(file, "\\def\\prune{%hd}\n\\def\\regraft{%hd}\n\n", o[0], o[1]);
		if (findParentRows(t.C, nb, o[1], &rg_parent) == 0)
			fprintf(file, "\\def\\prpar{%hd}\n\\def\\rgpar{%hd}\n\n", (short) (pr_parent + n_l),
					(short) (-1));
		else
			fprintf(file, "\\def\\prpar{%hd}\n\\def\\rgpar{%hd}\n\n", (short) (pr_parent + n_l),
					(short) (rg_parent + n_l));
	}

// new node
	fprintf(file, "\\def\\newnode{%hd}\n\n", new_node);

}

