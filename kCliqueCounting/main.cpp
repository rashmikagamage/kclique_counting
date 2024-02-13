#include "Graph.h"
#include "Utility.h"
#include "Timer.h"
#include "popl.hpp"

using namespace std;
using namespace popl;

void print_usage() {
	// printf("Example usage: ./kDefectiveClique -g path_to_graph -k 3 -o -b\n");
}

int main(int argc, char *argv[]) {
#ifndef NDEBUG
	printf("**** kCliqueCounting (Debug) build at %s %s ***\n", __TIME__, __DATE__);
	printf("!!! You may want to define NDEBUG in Utility.h to get better performance!\n");
#else
	printf("**** kCliqueCounting (Release) build at %s %s ***\n", __TIME__, __DATE__);
#endif

	bool binary_input = false;
	string alg;
	double epsilon;

	OptionParser op("Allowed options");
	auto help_option = op.add<Switch>("h", "help", "\'produce help message\'");
	auto graph_option = op.add<Value<string>>("g", "graph", "\'path to input graph file\'");
	auto alg_option = op.add<Value<string>>("a", "alg", "\'algorithm name\' (exact | exact1 | exact2 | TDKE23 | TKDE23_opt | arbitrary)", "exact", &alg);
	auto k_option = op.add<Value<int>>("k", "k", "\'the value of k for k-clique\'");
	auto epsilon_option = op.add<Value<double>>("e", "epsilon", "\'the value of epsilon\'", 0.001, &epsilon);
	op.add<Switch>("b", "binary", "\'read the input graph from binary files b_adj.bin and b_degree.bin\'", &binary_input);

	op.parse(argc, argv);

	if(help_option->is_set()||argc <= 1) {
		cout << op << endl;
		if(argc <= 1) {
			print_usage();
			return 0;
		}
	}

	if(!graph_option->is_set()) {
		printf("!!! The argument -g is required! Exit !!!\n");
		return 0;
	}
	if(!k_option->is_set()) {
		printf("!!! The argument -k is required! Exit !!!\n");
		return 0;
	}

	Graph *graph = new Graph(graph_option->value().c_str());
	if(binary_input) graph->read_graph_binary();
	else graph->read_graph();

#ifndef NDEBUG
	printf("\t*** Finished reading graph\n");
#endif

	// if(strcmp(alg.c_str(), "degen") == 0) graph->kDefectiveClique_degen(); //degeneracy-based k-plex
	// else print_usage();
	graph->kCliqueCounting(k_option->value(), alg, epsilon);

	delete graph;

	printf("\n");
	return 0;
}
