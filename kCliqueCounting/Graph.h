/*
 * Graph.h
 *
 *  Created on: 24 Sep 2023
 *      Author: ljchang
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"
#include "Pivotor_Matrix.h"

class Graph {
private:
	std::string dir; //input graph directory
	ui n; //number of nodes of the graph
	ept m; //number of edges of the graph

	ept *pstart; //offset of neighbors of nodes
	ui *edges; //adjacent ids of edges

	ui *degree;
	ui *rid;
	char *matrix; // adjacency matrix for subgraphs
	char *vis;

	double *dp;
	double *prob;
	ui *alias;
	ui *neighbors;
	ui *order_by_color;
	std::vector<ui> *heaps;
	std::vector<double> *ubs;

	double total_cnt; // total number of exact k-clique counts
	double total_ub; // total upper bound of k-clique counts
	double total_exp; // total expected number of k-cliques

	std::mt19937 gen;

	double density_factor = 0;

public:
	Graph(const char *_dir) ;
	~Graph() ;

	void read_graph_binary() ;
	void read_graph() ;
	void kCliqueCounting(ui K, std::string alg, double epsilon) ;

private:
	ui degen(ui degree_threshold, ui *peel_sequence, ui *core, ListLinearHeap *heap) ;
	void core_shrink_graph(ui core_threshold, ui &n, ept &m, ui *peel_sequence, ui *core, ept *&pstart, ui *&edges) ;
	void reorganize_adjacency_lists(ui n, ui *peel_sequence, ui *edges) ;

	ept get_number_of_edges(ui u, char *exists, ept *pstart, ept *pend, ui *edges) ;
	void get_subgraph_statistics(ept *pstart, ept *pend, ui *edges, ui K, char *vis) ;
	ept extract_subgraph(ui u, ui *out_mapping, ept *pend) ;
	bool matrix_degen(ui subgraph_n, ui *peel_sequence, ui &peel_sequence_n, ui threshold, ept &number_of_edges) ;
	ui matrix_coloring(ui subgraph_n, ui *peel_sequence, ui peel_sequence_n, ui *local_color) ;
	void reorder_peel_sequence_based_on_color(ui peel_sequence_n, ui *peel_sequence, ui *local_color, ui *order_by_color) ;
	double construct_DP_table_color_path(ui subgraph_n, ui peel_sequence_n, ui *peel_sequence, ui K) ;
	double construct_DP_table_color_path_and_alias(ui subgraph_n, ui peel_sequence_n, ui *peel_sequence, ui K, ui *lists) ;
	void print_matrix(ui n, char *matrix) ;
	void print_matrix(ui subgraph_n, char *matrix, ui ids_n, ui *ids) ;
	void print_array(const char *str, ui *array, ui array_n) ;
	void print_array(const char *str, ui *array, ui array_n, ui *out_mapping) ;
	void print_array(const char *str, char *array, ui array_n) ;
	void print_dp(ui peel_sequence_n, ui K) ;
	void assert_unique_and_bound(const char *str, ui ids_n, ui *ids, ui subgraph_n) ;

	void get_out_neighbors(char *adj, ui ids_n, ui *ids, ui &peel_sequence_n, ui *peel_sequence) ;
	ept count_number_of_edges(ui subgraph_n, ui ids_n, ui *ids) ;

	double build_alias(ui num, double *prob, ui *alias, ui *lists) ;
	void sample_a_batch(std::vector<ui> &nodes, double *ub, ul sample_size, ul &n_cliques,
		ul &n_samples, ui *peel_sequence, ui *out_mapping, ui *local_color, ept *pend, bool opt, ui K) ;
	void sample_a_batch_opt(std::vector<ui> &nodes, double *ub, ul sample_size, ul &n_cliques,
		ul &n_samples, ui *peel_sequence, ui *out_mapping, ui *local_color, ept *pend, ui K) ;
	ui sample_once(ui subgraph_n, ui peel_sequence_n, ui *peel_sequence, ui *cliques, ui K) ;
	ui sample_once_opt(ui subgraph_n, ui peel_sequence_n, ui *peel_sequence, ui *cliques, ui K) ;

	void sample_a_batch_for_one_subgraph(ui subgraph_n, ui peel_sequence_n, ui *order_by_color, ui *clique, ui l, ul &success, ul &total) ;
	void push_to_heap(std::vector<ui> &hh, ui peel_sequence_n, ui *peel_sequence, ui *out_mapping, ui success, ui trails, ui l) ;
	ept get_subgraph(ui ids_n, ui *ids, std::unordered_set<ui> *hash_tables, ui *out_mapping, ept *pend) ;
	void process_a_subgraph(ui leading_vertex, ui subgraph_n, ui peel_sequence_n, ui *peel_sequence, ui *out_mapping, ui l,
			Pivotor_Matrix *pm, ui *local_color, bool record_time, ul &total_sample_n, ul &total_sample_time,
			ui start_idx, std::vector<double> &density_thresholds, bool opt, ui *ids, double dd) ;
};

#endif /* GRAPH_H_ */
