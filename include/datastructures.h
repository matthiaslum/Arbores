/*
 *  datastructures.h
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

#ifndef datastructures_h
#define datastructures_h

#include <stdlib.h>

enum Node_Colors {
	black, gray, white, undef
};

struct Data {
	char *name;
	short *M;
	int *segregating_sites;
	int n_seq;
	int n_sites;
};

struct Descendants {
	short *nodes;
	short n_desc;
};

struct ShrubARG {
	short *C;
	int n_rec;
	short n_nodes;
	short n_seq;
};

struct Tree {
	short *C;
	double *times;
	short n_nodes;
};

struct Tree_array_version{
    short C[15];
    double times[15];
    short n_nodes;
};

struct Smc {
	struct Tree *tree_path;
	short **opers; //recombination prune n regraft nodes, equal to number of recombinations
	int *sites; //recombination sites
	double *rec_times; //recombintion times of the Opers
	short path_len; //total number of trees
	short *tree_selector;	 // this is a full sequence length array of short    //can be obtained from the function
												 // integers taking values in {0,...,path_len-1}
												 // so that the tree at i'th site of the full sequence
												 // is tree_path[tree_selector[i]];
	short **global_index;  // global indexing is an indexing of the branch nodes
												 // such that each node in the whole path with unique time
												 // has also a unique index
	int selector_length;
	short *is_free; 			 // array of length equal length as the path. 1 indicates
												 // a free time 0 fixed.
};

struct Smc_array_version{
    struct Tree_array_version tree_path[8];
    short opers[7][2];
    int sites[8];
    double rec_times[7];
    short path_len;
    short is_free[8];
};



struct SimplePath {
	struct Tree *trees;
	short **opers;
	long *sites;
	int length;
};

struct SamplingSet {
	struct Smc *paths;
	long *n_paths;
	short *set_selector; 			// this is *n_paths element array whose interpretation is that
														// the path paths[i] is sampled from a set whose cardinality is
														// set_cardinalities[set_selector[i]] ...
	long *set_cardinalities; // this is the collection of all possible set cardinalities. It
													 // is n_sets elements long (n_sets is small, i.e. 2)
	short n_sets;
};

struct ForcingConfig {
	short ** indicators;
	short n_configs;
	short length;
};

struct ShrubData {
	short *C;
	short *M;
	int *rowlab;
	int *collab;
	int m;
	int n;
	int nni; // next node index, i.e. node count
	int *rec_counts; // array for number or recombinations per SHRUB recombination step
	int rec_index; // number of SHRUB recombination steps
};

struct BridgePoint {
	int start;
	int end;
};

struct BridgePoints {
	struct BridgePoint *points; // initial and terminal sites of the bridging points
	int length; // number of segments
};

struct State {
	struct Smc path;
	long iteration;
};

struct Stats {
	double log_likelihood;
	short n_recomb;
	double alpha;
};

struct Conditioning {
	struct Tree tl;
	struct Tree tr;
	short *M; // data over the segment
	int *sites; // indices of the segregating sites in the segment
	short n_seq;
	int length;
	short backward;
	short twosided;
};

struct Conditioning_array_version{
    struct Tree_array_version tl;
    struct Tree_array_version tr;
    short M[40];
    int sites[10];
    short n_seq;
    int length;
    short backward;
    short twosided;
};

struct IntVector {
	int *values;
	int length;
};

struct LongVector {
	long *values;
	long length;
};

struct ShortVector {
	short *v;
	long length;
};

struct DoubleVector {
	double *v;
	long length;
};

struct FindResult {
	int *values;
	int *indices;
	int length;
};

//struct PathSet {
//	double **F;
//	long **L;
//	short **O;
//	int *sites;
//	short *extra_indicator;
//};

struct CanonicalMatrix {
	short *data; // the data is a m-by-n matrix
	long n;
	long m; // number of rows, i.e. the number of canonised trees
};

struct UniqueRowsWithCounts {
	short *data;
	long n;
	long m;
	long *counts;
	long *inds;
	long *selector;
};

struct TimeParms {
	double *times;
	double *parms;
	int length;
};

struct AdjacencySet {
	struct Tree *trees;
	short *opers;
	long *n_entries;
	short length;
};

struct AdjacencySets {
	struct Tree ** tree_adj_set;
	short ** oper_adj_set;
	long *set_size;
	long n_sets;
	long n_active;
	short length;
	long *sites;
};

struct AdjSetUnion {
	struct Tree *trees;
	short *opers;
	long *links;
	long n_paths;
	short length;
};

struct LinkedSet {
	struct Tree *trees; // linked set contain n_trees UNIQUE trees
	short **oper; // for ith (out of n_trees) tree there is a set
								// of operations opers[i] (set_size[i]-by-2)
	long **links; // and links[i] is set_size[i]-by-1. The rationale
								// is that trees[i] is obtained by applying operation
								// (opers[i][j],opers[i][j+set_size[i]]) to the tree
								// domain[links[i][j]], for any j in 0...set_size[i]-1
								// NB: It may be that instead of single operation,
								// multiple operations are needed, in which case,
								// opers[i] is set_size[i]-by-2*n_op and links[i] is
								// set_size[i]-by-n_op, where n_op is the required number
								// operations
	long *set_size;
	long *sites;
	short length;
	long n_trees;
};

struct LikelihoodData {
	double log_likelihood;
	double *increments;
	double *total_branch_length;
	double *mutated_branch_lengths;
	long length;
};

struct LinkedSetArray {
	struct LinkedSet *sets;
	short length;
	short completed;
};

struct MRCA {
	double **times;
	short nl;
};

struct ShortRowSortedMtx {

	short *data;
	long m;
	long n;
	long *i_sort;

};

struct Recombinations {
	double log_p;
	short n_recombinations;
};

struct DoubleShortPair {
	double v;
	short i;
};

struct BranchLength {
	double total_length;
	double mutated_length;
	double *lenghts;
	short n;
};

struct Parameters {
	double mu;
	double rho;
	int n_eff; // effective population size
	short verb;
};

struct RegraftTimeData {
	double p_rg_time;
	short n_active;
};

struct SmcPriorData {
	double density;
	double *components;
	double *increments;
	int length;
	short number_of_recombinations;
};

struct JitterSychro {
	short accept;
	double new_time;
	double *rec_times;
	short length;
};

struct SegmentSamplerSynchro {
	short **C;
	double **times;
	short **opers;
	int *sites;
	short length;
	short nl;
};

struct InitTimeSynchro {
	double *times;
	short *indices;
	short length;
};

struct FreeTimeSynchro {
	double *times;
	short length;
};

struct RecTimeSynchro {
	double *times;
	short length;
};

struct AcceptanceSynchro {
	short accept;
};

struct MCMCDiagnostics {

    double proposed_log_likelihood;
    double current_log_likelihood;
    double proposed_log_prior;
    double current_log_prior;
    double proposed_free_time_density;
    double current_free_time_density;
    double proposed_recombination_density;
    double current_recombination_density;
    short proposed_number_of_free_times;
    short current_number_of_free_times;
    short proposed_number_of_recombinations;
    short current_number_of_recombinations;
    double alpha;
    short accept_indicator;
    short jitter_step;
    double log_posterior;
    double log_prior;
    double log_likelihood;
    double cardinality_ratio;
    short irreducibility;
    double u;
    char indicators[7];
};

struct NonIdentityOps {
	short **ops;
	double *rec_times;
	short nops;
};

struct MCMCSummary {
	struct Smc path;
	struct MCMCDiagnostics data;
	short full_scan;
};

struct segment_output{
    struct Smc new_segment;
    short accept_indicator;
};

struct arraySegmentOutput{
    struct Smc_array_version new_segment;
    short accept_indicator;
};

enum Boundary {
	up, down
};

struct AdjustedPath {
	struct Tree *trees;
	short **opers;
	short *free_times;
};

struct GlobalTimes {
	double *t;
	short *step;
	short n_times;
};

struct LimitData {
	short *inds;
	double *times;
	short *is_free;
	short *glb_inds;
	short step;
};

struct LimitSummary {
	struct LimitData *limits;
	struct GlobalTimes glb_times;
	short *is_new;
	short *conditioning_index;
	short *freedom_indicator;
	double *adjusted_times;
	short *new_node;
	short length;
};

struct SamplingSetData {
	struct SamplingSet set;
	long i_curr;
};

enum Step_Type {
	prune, merge, terminate
};

struct SearchPath {
	struct Tree *trees; // array of trees of length 'length'
	short * opers; // length-by-2 array of operations
	short length; // length of the path
};

struct OperationSet {
	short *opers;
	long n_ops;
};

struct RegraftTimes {
	double *times;
	int n_times;
};

struct Roots {
	short *roots;
	short n_roots;
};

struct Removed_Nodes {
	short *nodes;
	short n_nodes;
};

struct SearchData {
	struct AdjacencySet set;
	struct SearchPath path;
};

struct AcceptanceProbabilityData {
	double alpha;
	double *log_likelihood;
};

#endif /* datastructures_h */
