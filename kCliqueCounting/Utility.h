#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <random>
#include <tuple>
#include <queue>
#include <set>
#include <unordered_set>

#define NDEBUG
#include <cassert>

using ui = unsigned int; // vertex type
using ept = unsigned int; // edge pointer type; unsigned int can be used to process upto two billion undirected edges
using ul = unsigned long long;

#define pb push_back
#define mp make_pair

#define mmax(a,b) ((a)>(b)?(a):(b))
#define mmin(a,b) ((a)<(b)?(a):(b))

class Utility {
public:
	static FILE *open_file(const char *file_name, const char *mode) {
		FILE *f = fopen(file_name, mode);
		if(f == nullptr) {
			printf("Can not open file: %s\n", file_name);
			exit(1);
		}

		return f;
	}

	static std::string integer_to_string(long long number) {
		std::vector<ui> sequence;
		if(number == 0) sequence.pb(0);
		while(number > 0) {
			sequence.pb(number%1000);
			number /= 1000;
		}

		char buf[5];
		std::string res;
		for(ui i = sequence.size();i > 0;i --) {
			if(i == sequence.size()) sprintf(buf, "%u", sequence[i-1]);
			else sprintf(buf, ",%03u", sequence[i-1]);
			res += std::string(buf);
		}
		return res;
	}

	static void print_neighbors(ui u, const ui *pstart, const ui *pend, const ui *edges) {
		std::vector<ui> neighbors;
		for(ui i = pstart[u];i < pend[u];i ++) neighbors.push_back(edges[i]);
		//sort(neighbors.begin(), neighbors.end());
		printf("neighbors of %u:", u);
		for(ui i = 0;i < neighbors.size();i ++) printf(" %u", neighbors[i]);
		printf("\n");
	}

	static void print_array(const char *str, const ui *array, ui idx_start, ui idx_end, ui l) {
		for(ui i = 0;i < l;i ++) printf(" ");
		printf("%s:", str);
		for(ui i = idx_start;i < idx_end;i ++) printf(" %u", array[i]);
		printf("\n");
	}
};

#endif
