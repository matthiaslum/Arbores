/*
 * exhaustiveSearch.c
 *
 *  Created on: 16.11.2016
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
#include <assert.h>
#include <math.h>
#include "constants.h"
#include "exhaustiveSearch.h"
#include "utils.h"

int termination_counter = 0;

struct AdjacencySet exhaustiveSearch(struct Tree t, struct Conditioning cond, short col,
		short forced, short max_len) {

	struct AdjacencySet aSet;
	struct SearchPath path;
	long n;

	path = createSearchPath(t, max_len);

	aSet.trees = malloc(sizeof(struct Tree) * MAX_ADJ_SET_SIZE * max_len);
	aSet.opers = malloc(sizeof(short) * 3 * MAX_ADJ_SET_SIZE * (max_len - 1));
	aSet.n_entries = malloc(sizeof(long));
	aSet.length = max_len;
	aSet.n_entries[0] = 0;

	assert(checkTree(t) == 1);

	step(merge, path, aSet, cond, col, max_len, forced);
	step(prune, path, aSet, cond, col, max_len, forced);

	for (int i = 0; i < max_len; i++)
		deleteTree(path.trees[i]);
	free(path.trees);
	free(path.opers);

	// Release the extra memory allocated above
	n = *aSet.n_entries;
	if (n == 0) {
		free(aSet.trees);
		free(aSet.opers);
	} else {
		for (long i = 0; i < n * max_len; i++)
			aSet.trees[i] = aSet.trees[i % n + i / n * MAX_ADJ_SET_SIZE];
		aSet.trees = (struct Tree*) realloc(aSet.trees, sizeof(struct Tree) * n * max_len);
		for (long i = 0; i < n * 3 * (max_len - 1); i++)
			aSet.opers[i] = aSet.opers[i % n + i / n * MAX_ADJ_SET_SIZE];
		aSet.opers = (short*) realloc(aSet.opers, sizeof(short) * n * 3 * (max_len - 1));
	}
	return aSet;
}

/* This function creates the whole tree path in advance, but the trees are only
 * valid up to the tree path.length-1.
 */
struct SearchPath createSearchPath(struct Tree t, short length) {

	struct SearchPath path;

	path.trees = malloc(sizeof(struct Tree) * length);
	path.opers = malloc(sizeof(short) * (length - 1) * 3);
	for (int i = 0; i < length; i++)
		path.trees[i] = createCopy(t);
	path.length = 1;

	return path;
}

void deleteSearchPath(struct SearchPath path) {

	for (int i = 0; i < path.length; i++)
		deleteTree(path.trees[i]);
	free(path.trees);
	free(path.opers);

}

void step(enum Step_Type st, struct SearchPath path, struct AdjacencySet a,
		struct Conditioning cond, short col, short max_len, short extra) {

//	struct Tree *visu;
	struct Tree last_tree, new_tree, tmp_tree;
	struct OperationSet opSet;
	struct RegraftTimes regraft_times;
	struct SearchPath new_path;
	long ops_max = MAX_ADJ_SET_SIZE, n;
	short new_node, n_l;
	short *o;

//	visu = malloc(sizeof(struct Tree) * 3);

	last_tree = path.trees[path.length - 1];
	n_l = (last_tree.n_nodes + 1) / 2;
	n = *a.n_entries;
//	for (int i = 0; i < path.length; i++)
//		assert(checkTree(path.trees[i]) == 1);

	if (st == terminate) {
		if (checkFeasibility(last_tree, cond.M, (int) col) == 1) {
			for (int i = 0; i < path.length - 1; i++) {
//				assert(checkTree(path.trees[i]) == 1);
				a.trees[n + i * ops_max] = createCopy(path.trees[i]);
				a.opers[n + (3 * i + 0) * ops_max] = path.opers[i];
				a.opers[n + (3 * i + 1) * ops_max] = path.opers[i + max_len - 1];
				a.opers[n + (3 * i + 2) * ops_max] = path.opers[i + 2 * (max_len - 1)];
			}
			a.trees[n + (path.length - 1) * ops_max] = createCopy(path.trees[path.length - 1]); // last tree
//			assert(checkTree(path.trees[path.length - 1]) == 1);
			a.n_entries[0]++;
		}
	} else {

		// create the sets of extending operations
		opSet =
				st == merge ?
						mergeSearchSpace(path, cond, col, extra) : pruneSearchSpace(path, cond, col);
//		printf("%s, %ld\n", st == merge ? "Merge set size " : "Prune set size ", opSet.n_ops);

		// create new search path
		new_path.trees = path.trees;
		new_path.opers = path.opers;
		new_path.length = path.length + 1;

		for (int i = 0; i < opSet.n_ops; i++) {

			// add the latest operation to the path
			new_path.opers[new_path.length - 2] = opSet.opers[i];
			new_path.opers[new_path.length - 2 + max_len - 1] = opSet.opers[i + ops_max];
			new_path.opers[new_path.length - 2 + 2 * (max_len - 1)] = -1;

			new_tree = new_path.trees[new_path.length - 1];

			if (opSet.opers[i] == -1 && opSet.opers[i + ops_max] == -1) {
				// for no-op, extend the path with a copy of the last tree and
				// take the recursion step

				copyTree(&new_tree, last_tree);

				if (new_path.length == max_len) {
					step(terminate, new_path, a, cond, col, max_len, extra);
				} else {
					step(merge, new_path, a, cond, col, max_len, extra);
					step(prune, new_path, a, cond, col, max_len, extra);
				}

			} else {
				// non no-op case

				// create a temporary tree for the SPR outcome
				tmp_tree = createCopy(last_tree);
//				writeTikzTexFileForTreePath(&last_tree, NULL, NULL, NULL, 1);
				// execute the SPR operation
				o = pickOp(opSet.opers, (long) i, ops_max);
				new_node = SPR(tmp_tree, o);
				timeSorting(tmp_tree, NULL);
//				assert(checkTree(tmp_tree) == 1);
//				printf("Tree before operation %hd %hd \n", o[0], o[1]);
//				printTree(last_tree);

				// find the possible regrafting times
				regraft_times = getRegraftTimeRange(last_tree, o);
//				printf("%i %hd %hd %i\n",i,o[0],o[1],regraft_times.n_times);
				free(o);

				// iterate over the regrafting times and deepen the recursion
				for (int j = 0; j < regraft_times.n_times; j++) {

					// adjust new node time
					tmp_tree.times[new_node - n_l] = regraft_times.times[j];
					copyTree(&new_tree, tmp_tree);
					// ensure time ordering
					timeSorting(new_tree, NULL);
//					assert(checkTree(new_tree)==1);
//					if (checkTree(new_tree) != 1) {
//						writeTikzTexFileForTreePath(NULL, &last_tree, NULL, NULL, NULL, NULL, NULL,
//								1);
//						printf("Before sorting\n");
//						printTree(tmp_tree);
//						printf("Failing tree\n");
//						printTree(new_tree);
//						printDoubleArray(regraft_times.times, 1, regraft_times.n_times, 1);
//						assert(false);
//					}

					new_path.opers[new_path.length - 2 + 2 * (max_len - 1)] = j;

					if (new_path.length == max_len) {
						step(terminate, new_path, a, cond, col, max_len, extra);
					} else {
						step(merge, new_path, a, cond, col, max_len, extra);
						step(prune, new_path, a, cond, col, max_len, extra);
					}
				}
				free(regraft_times.times);
				free(tmp_tree.C);
				free(tmp_tree.times);
			}
		}
		free(opSet.opers);
	}
}

short isNoop(short *o) {
	return (o[0] == -1 && o[1] == -1) ? 1 : 0;
}

/* Return a 1-by-2 array of short integers to represent an SPR operation.
 * short *ops 	- n_ops-by-2 array of operations
 * short row		- the row from which the operation is picked
 * short n_ops  - number of rows in ops
 */
short *pickOp(short *ops, long row, long n_ops) {
	short *o = malloc(sizeof(short) * 2);
	o[0] = ops[row];
	o[1] = ops[row + n_ops];
	return o;
}

struct RegraftTimes getRegraftTimeRange(struct Tree t, short *op) {

	struct RegraftTimes r_times;
	short n_l, n_b, regraft_parent;
	double min_time, max_time, pr_node_time, rg_node_time, next_time, curr_time, avg_time;
	double *sorted_times;

	n_l = (t.n_nodes + 1) / 2;
	n_b = n_l - 1;

	/* regrafting to root */
	if (op[1] == t.n_nodes - 1) {
		r_times.n_times = 1;
		r_times.times = malloc(sizeof(double));
		r_times.times[0] = t.times[op[1] - n_l] + 1;
		return r_times;
	}

	/* identity operation */
	if (op[0] == -1 && op[1] == -1) {
		r_times.n_times = 1;
		r_times.times = NULL;
		return r_times;
	}

	r_times.n_times = 0;
	r_times.times = malloc(sizeof(double) * (n_b + 1));

	/* general case */
	sorted_times = malloc(sizeof(double) * (n_b + 1));
	sorted_times[0] = 0;
	for (int i = 0; i < n_b; i++)
		sorted_times[i + 1] = t.times[i];

	pr_node_time = op[0] < n_l ? 0 : t.times[op[0] - n_l];
	rg_node_time = op[1] < n_l ? 0 : t.times[op[1] - n_l];
	regraft_parent = findParent(t.C, n_b, n_l, op[1]);

	// bounds for regrafting time
	min_time = pr_node_time > rg_node_time ? pr_node_time : rg_node_time;
	max_time = t.times[regraft_parent - n_l];
	for (int j = 0; j < n_b; j++) {
		curr_time = sorted_times[j];
		next_time = sorted_times[j + 1];
		avg_time = (curr_time + next_time) / (double) 2;
		//		if (min_time - curr_time < TOL)  // in the range
		if (avg_time > min_time && avg_time < max_time) 			// in the range
			r_times.times[r_times.n_times++] = avg_time;
		//		if (next_time - max_time > TOL)
		if (avg_time > max_time)
			break;
	}
	free(sorted_times);
	return r_times;
}

struct OperationSet mergeSearchSpace(struct SearchPath path, struct Conditioning cond,
		short col, short extra_rec) {

	struct OperationSet opSet;
	struct Tree t;
	struct Roots r;
	struct Descendants descs;
	enum Node_Colors *colors, **clrs;
	long ops_max = MAX_ADJ_SET_SIZE;
	short *complement, *subset;
	short nl, nb, pr_root, rg_root, pr_node, rg_node, pr_prnt, rg_prnt, black_leaf_count,
			n_comp, n_sub;
	double pr_node_time;

	clrs = malloc(sizeof(enum Node_Colors *));

	t = path.trees[path.length - 1];
	nl = (t.n_nodes + 1) / 2;
	nb = nl - 1;
	colors = colorTree(t, cond.M, col);
	r = findSingleColorSubtreeRoots(t, colors, black);

	opSet.opers = malloc(sizeof(short) * 2 * ops_max);
	opSet.n_ops = 0;

	descs.nodes = malloc(sizeof(short) * t.n_nodes);

// Allow identity operation
	opSet.opers[0] = -1;
	opSet.opers[0 + ops_max] = -1;
	opSet.n_ops++;

	if (r.n_roots > 1) {

		for (int i = 0; i < r.n_roots; i++) {
			for (int j = i + 1; j < r.n_roots; j++) {

				pr_root = r.roots[i];
				rg_root = r.roots[j];

				for (int ii = 0; ii < 2; ii++) {
					// a trick to repeat operation with prune and regraft swapped

					descs.nodes[0] = rg_root;
					descs.n_desc = 1;
					descs = descendants(t, rg_root, descs);
					/* iterate over the descendants of the rg_root */
					for (int k = 0; k < descs.n_desc; k++) {

						rg_prnt = findParent(t.C, nb, nl, descs.nodes[k]);
						pr_node_time = pr_root < nl ? 0.0 : t.times[pr_root - nl];

						/* Regrafting is possible only if regraft parent is older that prune node. */
						if (rg_prnt == -1 || t.times[rg_prnt - nl] > pr_node_time) {

							opSet.opers[opSet.n_ops] = pr_root;
							opSet.opers[opSet.n_ops + ops_max] = descs.nodes[k];
							opSet.n_ops++;
						}
					}
					// Swap prune and regraft nodes and repeat
					pr_root = r.roots[j];
					rg_root = r.roots[i];
				}
			}
		}
	} else if (extra_rec == 1 && r.n_roots == 1) {
		// This branch is entered only if both conditions are satisfied:
		//  1) There is only one black subtree (i.e. tree is compatible)
		//  2) Extra recombinations are allowed

		black_leaf_count = 0;
		for (int ii = 0; ii < nl; ii++)
			black_leaf_count = black_leaf_count + (colors[ii] == black ? 1 : 0);
		clrs[0] = colors;
		assert(r.n_roots == 1);

		// descendants of the ONLY black subtree
		descs.nodes[0] = r.roots[0];
		descs.n_desc = 1;
		descs = descendants(t, r.roots[0], descs);

		if (black_leaf_count == 1) {
			// the black subtree is a single leaf
			n_comp = t.n_nodes - 1;
			complement = subtreeComplement(descs, t.n_nodes);
			complement[n_comp - 1] = r.roots[0];
			clrs[0] = colors;
			printf("Data containing segregating sites with only one mutation are not supported. Exiting...\n");
			assert(false);
		} else {
			// If there is only one black subtree, then the nodes in the
			// complement are either white or gray. Prune-search already
			// takes care of the pruning of white nodes and gluing back
			// to white (non-descendant) or gray node. Therefore the only
			// thing that needs to be done is pruning gray node and gluing
			// it back anywhere (non-descendant).
			n_comp = t.n_nodes - descs.n_desc;
			complement = subtreeComplement(descs, t.n_nodes);
		}
		n_sub = n_comp;
		subset = complement;

		for (int i_set = 0; i_set < 2; i_set++) {
			// This loop takes care of pruning and regrafting
			// 1) inside a black subtree (second iteration)
			// 2) inside the complement of a black subtree (first iteration)

			// Iterate over all pairs of different nodes outside the black subtree
			for (int i = 0; i < n_sub; i++) {
				for (int j = i + 1; j < n_sub; j++) {

					pr_node = subset[i];
					rg_node = subset[j];

					for (int i_swap = 0; i_swap < 2; i_swap++) {

						pr_prnt = findParent(t.C, nb, nl, pr_node);
						rg_prnt = findParent(t.C, nb, nl, rg_node);
						// Regraft is possible only if the rg parent is older that prune node OR
						// if regrafting to the (global) root.
						if (rg_prnt > pr_node || rg_prnt == -1) {

							// In addition you are not allowed to:
							//  1) Regraft to your parent
							//  2) Regraft to your sibling
							if ((rg_node != pr_prnt) && (rg_prnt != pr_prnt) && pr_node != rg_prnt) {

								opSet.opers[opSet.n_ops] = pr_node;
								opSet.opers[opSet.n_ops + ops_max] = rg_node;
								opSet.n_ops++;
								assert(pr_node < t.n_nodes && rg_node <= t.n_nodes);
							}
						}
						// Swap the prune and regraft nodes
						pr_node = subset[j];
						rg_node = subset[i];
					}
				}
			}
			// Switch to considering the inside of the black subtree
			n_sub = descs.n_desc;
			subset = descs.nodes;
		}
		free(complement);
	}
	free(clrs);
	free(colors);
	free(r.roots);
	free(descs.nodes);

	return opSet;
}

struct OperationSet pruneSearchSpace(struct SearchPath path, struct Conditioning cond,
		short col) {

	struct OperationSet opSet;
	struct Tree t, t_cpy;
	struct Roots r;
	struct Removed_Nodes removed;
	short n_l, n_b, rg_parent, prune_parent;
	enum Node_Colors *colors1, *colors2;
	short not_black, old_enough, black_root, global_root;
	long ops_max = MAX_ADJ_SET_SIZE;

	t = path.trees[path.length - 1];
	n_l = (t.n_nodes + 1) / 2;
	n_b = n_l - 1;

	t_cpy = createCopy(t);
	colors1 = colorTree(t_cpy, cond.M, col);

	opSet.opers = malloc(sizeof(short) * 2 * ops_max);
	opSet.n_ops = 0;

	removed.nodes = malloc(sizeof(short) * t.n_nodes);

	/* Iterate over all nodes to see if node is the root of a white subtree */
	for (int i = 0; i < t.n_nodes; i++) {

		if (colors1[i] == white) {

			prune_parent = findParent(t.C, n_b, n_l, i);

			copyTree(&t_cpy, t);

			// Mark white root for removal
			removed.nodes[0] = i;
			removed.n_nodes = 1;
			// Recursively remove the white subtree,
			// starting from the root and removing children
			removed = removeRecursion(t_cpy, removed);

			colors2 = colorTree(t_cpy, cond.M, col);
			r = findSingleColorSubtreeRoots(t_cpy, colors2, black);
			/* Possible nodes for regrafting are:
			 1) Any white nodes
			 2) Any gray nodes
			 3) Any black subtree root node
			 */
			for (int j = 0; j < t.n_nodes; j++) { // This loop makes things O(N^2)
				if (j != prune_parent) {
					// don't regraft to parent
					rg_parent = findParent(t_cpy.C, n_b, n_l, j);

					if (rg_parent != prune_parent) {
						// don't regraft to sibling
						global_root = (j == t.n_nodes - 1) ? 1 : 0; // try regrafting to root

						if (rg_parent != -1) {
							// regraft node has a parent
							//
							// NB: This branch takes care of the case when attempting
							// to regraft to a descendant of the pruned white subtree:
							// after pruning the white subtree, all the nodes in the
							// pruned subtree are parentless, therefore this branch
							// is not entered.

							not_black = (colors2[j] != black) ? 1 : 0;
							old_enough = (rg_parent >= i) ? 1 : 0;
							black_root = isroot(t_cpy, j, r);
							if ((not_black == 1 && old_enough == 1 && j != i)
									|| (black_root == 1 && j != i && old_enough == 1)) {
								opSet.opers[opSet.n_ops] = i;
								opSet.opers[opSet.n_ops + ops_max] = j;
								opSet.n_ops++;
								if (i > t.n_nodes || j > t.n_nodes)
									printf("ERROR: Invalid operation generated\n");

							}

						} else if (global_root == 1 && i != j) {
							// regraft does not have a parent
							// Regrafting to global root provided that global root was not the original
							// prune node.
							// NB: The case when attempting to regraft to descendants of the pruned subtree
							// will not lead to this branch either. This branch is entered only when attempting
							// to regraft to the global root, which never can be a descendant of anything.
							opSet.opers[opSet.n_ops] = i;
							opSet.opers[opSet.n_ops + ops_max] = j;
							opSet.n_ops++;
							if (i > t.n_nodes || j > t.n_nodes)
								printf("ERROR: Invalid operation generated\n");
						}
					}
				}
			}
			free(r.roots);
			free(colors2);
		}
	}
	free(removed.nodes);
	free(t_cpy.C);
	free(t_cpy.times);
	free(colors1);
	return opSet;
}

short isroot(struct Tree t, short node, struct Roots r) {
	for (int i = 0; i < r.n_roots; i++)
		if (r.roots[i] == node)
			return 1;
	return 0;
}

struct Roots findSingleColorSubtreeRoots(struct Tree t, enum Node_Colors *colors,
		enum Node_Colors clr) {

	struct Roots r;
	short n_l, n_b, parent;
	n_l = (t.n_nodes + 1) / 2;
	n_b = n_l - 1;
	r.roots = malloc(sizeof(short) * t.n_nodes);
	r.n_roots = 0;

	/* Iterate over all nodes */
	for (int i = 0; i < t.n_nodes; i++) {
		// if current node is of chosen color AND its parent is different
		// color, then node is a root
		parent = findParent(t.C, n_b, n_l, i);
		if (parent == -1) {
			if (colors[i] == clr)
				r.roots[r.n_roots++] = (short) i;
		} else {
			if (colors[i] == clr && colors[parent] != clr)
				r.roots[r.n_roots++] = (short) i;
		}
	}
	return r;
}

struct Descendants descendants(struct Tree t, short node, struct Descendants descs) {

	short n_l = (t.n_nodes + 1) / 2, n_b = n_l - 1;

	/* terminate recursion */
	if (node < n_l)
		return descs;

	/* left child */
	if (t.C[node - n_l] >= 0) {
		descs.nodes[descs.n_desc++] = t.C[node - n_l];
		descs = descendants(t, t.C[node - n_l], descs);
	}
	/* right child */
	if (t.C[node - n_l + n_b] >= 0) {
		descs.nodes[descs.n_desc++] = t.C[node - n_l + n_b];
		descs = descendants(t, t.C[node - n_l + n_b], descs);
	}

	return descs;
}

short *subtreeComplement(struct Descendants desc, short n_n) {

	short ptr = 0;
	short *selector, *complement;

	complement = malloc(sizeof(short) * n_n);
	selector = malloc(sizeof(short) * n_n);

	for (int i = 0; i < n_n; i++)
		selector[i] = 1;

	for (int i = 0; i < desc.n_desc; i++)
		selector[desc.nodes[i]] = 0;

	for (int i = 0; i < n_n; i++)
		if (selector[i] == 1)
			complement[ptr++] = i;

	free(selector);
	return complement;
}

struct Removed_Nodes removeRecursion(struct Tree t, struct Removed_Nodes remove) {

	short *removed_children;
	short n_l, n_b, n_rm_children = 0, any = 0;

	n_l = (t.n_nodes + 1) / 2;
	n_b = n_l - 1;

// array for the children of all removed nodes
	removed_children = malloc(sizeof(short) * t.n_nodes);

// iterate over removed nodes and remove their children
	for (int i = 0; i < remove.n_nodes; i++) {
		if (remove.nodes[i] >= n_l) { // leaves don't have children
			if (t.C[remove.nodes[i] - n_l] != -1) { // removed node has a left child
				any = 1;
				removed_children[n_rm_children++] = t.C[remove.nodes[i] - n_l];
				t.C[remove.nodes[i] - n_l] = -1;
			}
			if (t.C[remove.nodes[i] - n_l + n_b] != -1) {
				any = 1;
				removed_children[n_rm_children++] = t.C[remove.nodes[i] - n_l + n_b];
				t.C[remove.nodes[i] - n_l + n_b] = -1;
			}
		}
	}

// Remove edges to parents
	for (int i = 0; i < n_b; i++)
		for (int j = 0; j < remove.n_nodes; j++)
			if (t.C[i] == remove.nodes[j])
				t.C[i] = -1;
			else if (t.C[i + n_b] == remove.nodes[j])
				t.C[i + n_b] = -1;

	if (any == 0) {
		free(removed_children);
		return remove;
	}

	for (int i = 0; i < n_rm_children; i++)
		remove.nodes[remove.n_nodes++] = removed_children[i];

	free(removed_children);

	remove = removeRecursion(t, remove);
	return remove;
}

/* This function if used solely for debugging purposes */
short checkOperationValidity(struct Tree t, short *o) {

	struct Descendants descs;

	if (o[0] == -1 && o[1] == -1)
		return 0;

	descs.nodes = malloc(sizeof(short) * t.n_nodes);
	descs.n_desc = 0;

// find the descendants of the pruning node
	descs.nodes[descs.n_desc++] = o[0];
	descs = descendants(t, o[0], descs);

	for (int i = 0; i < descs.n_desc; i++)
		if (o[1] == descs.nodes[i])
			return -1;
	return 0;
}

