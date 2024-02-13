/*
 * Graph.cpp
 *
 *  Created on: 24 Sep 2023
 *      Author: ljchang
 */

#include "Graph.h"
#include "CTPrune.h"
#include "Pivotor_Matrix.h"

using namespace std;

Graph::Graph(const char *_dir) {
	dir = string(_dir);

	n = m = 0;

	pstart = nullptr;
	edges = nullptr;

	degree = nullptr;
	rid = nullptr;
	matrix = nullptr;
	vis = nullptr;
	dp = nullptr;
	prob = nullptr;
	alias = nullptr;
	neighbors = nullptr;
	order_by_color = nullptr;
	heaps = nullptr;
	ubs = nullptr;

	gen.seed(12345);

	total_cnt = total_ub = total_exp = 0;
}

Graph::~Graph() {
	if(pstart != nullptr) {
		delete[] pstart;
		pstart = nullptr;
	}
	if(edges != nullptr) {
		delete[] edges;
		edges = nullptr;
	}
	if(degree != nullptr) {
		delete[] degree;
		degree = nullptr;
	}
	if(rid != nullptr) {
		delete[] rid;
		rid = nullptr;
	}
	if(matrix != nullptr) {
		delete[] matrix;
		matrix = nullptr;
	}
	if(vis != nullptr) {
		delete[] vis;
		vis = nullptr;
	}
	if(dp != nullptr) {
		delete[] dp;
		dp = nullptr;
	}
	if(prob != nullptr) {
		delete[] prob;
		prob = nullptr;
	}
	if(alias != nullptr) {
		delete[] alias;
		alias = nullptr;
	}
	if(neighbors != nullptr) {
		delete[] neighbors;
		neighbors = nullptr;
	}
	if(order_by_color != nullptr) {
		delete[] order_by_color;
		order_by_color = nullptr;
	}
	if(heaps != nullptr) {
		delete[] heaps;
		heaps = nullptr;
	}
	if(ubs != nullptr) {
		delete[] ubs;
		ubs = nullptr;
	}
}

void Graph::read_graph_binary() {
	printf("# Start reading graph, Require files \"b_degree.bin\" and \"b_adj.bin\"\n");
	FILE *f = Utility::open_file((dir + string("/b_degree.bin")).c_str(), "rb");

	ui tt;
	fread(&tt, sizeof(int), 1, f);
	if(tt != sizeof(int)) {
		printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, f);
	fread(&m, sizeof(int), 1, f);

	printf("\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	ui *degree = new ui[n];
	fread(degree, sizeof(int), n, f);

#ifndef NDEBUG
	long long sum = 0;
	for(ui i = 0;i < n;i ++) sum += degree[i];
	if(sum != m) printf("m not equal sum of degrees\n");
#endif

	fclose(f);

	f = Utility::open_file((dir + string("/b_adj.bin")).c_str(), "rb");

	if(pstart == nullptr) pstart = new ept[n+1];
	if(edges == nullptr) edges = new ui[m];

	pstart[0] = 0;
	for(ui i = 0;i < n;i ++) {
		if(degree[i] > 0) {
			fread(edges+pstart[i], sizeof(int), degree[i], f);

			// remove self loops and parallel edges
			ui *buff = edges+pstart[i];
			sort(buff, buff+degree[i]);
			ui idx = 0;
			for(ui j = 0;j < degree[i];j ++) {
				if(buff[j] >= n) printf("vertex id %u wrong\n", buff[j]);
				if(buff[j] == i||(j > 0&&buff[j] == buff[j-1])) continue;
				buff[idx ++] = buff[j];
			}
			degree[i] = idx;
		}

		pstart[i+1] = pstart[i] + degree[i];
	}

	fclose(f);

	delete[] degree;
}

void Graph::read_graph() {
	printf("# Start reading graph from an edgelist file, Require files \"edges.txt\"\n");
	printf("# Note that this function is not optimized. Reading from a binary file will be faster\n");
	FILE *f = Utility::open_file((dir + string("/edges.txt")).c_str(), "r");

	fscanf(f, "%u%u", &n, &m);
	m *= 2;
	printf("*\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	vector<pair<ui,ui> > vp;
	for(ui i = 0;i < m/2;i ++) {
		ui a, b;
		fscanf(f, "%u%u", &a, &b);
		if(a >= n || b >= n) {
			printf("!!! Vertex IDs must be between 0 and n-1. Exit !!!\n");
			return ;
		}
		vp.pb(mp(a,b));
		vp.pb(mp(b,a));
	}
	sort(vp.begin(), vp.end());

	if(pstart != nullptr) delete[] pstart;
	pstart = new ept[n+1];
	if(edges != nullptr) delete[] edges;
	edges = new ui[m];

	pstart[0] = 0;
	ui idx = 0;
	for(ui i = 0;i < n;i ++) {
		pstart[i+1] = pstart[i];
		while(idx < vp.size()&&vp[idx].first == i) edges[pstart[i+1] ++] = vp[idx ++].second;
	}

	fclose(f);

#ifndef NDEBUG
	printf("Finished reading graph\n");
#endif
}

void Graph::kCliqueCounting(ui K, string alg, double epsilon) {
	if(K <= 3 || K >= n) {
		printf("K(=%u) is either too small or too large!\n", K);
		return ;
	}

	Timer t;

	if(degree == nullptr) degree = new ui[n];
	if(vis == nullptr) vis = new char[n];
	if(rid == nullptr) rid = new ui[n];

	ui *peel_sequence = new ui[n];
	ui *core = new ui[n];
	ListLinearHeap *heap = new ListLinearHeap(n, n-1);
	ui *out_mapping = new ui[n];

	ui max_core = degen(K-1, peel_sequence, core, heap);
	core_shrink_graph(K-1, n, m, peel_sequence, core, pstart, edges);
	if(m&&false) {
		ept tm = m;
		CTPrune::truss_pruning(n, m, K-2, peel_sequence, out_mapping, rid, pstart, edges, degree, true);
		if(m < tm/2) {
			ui *edges_new = new ui[m];
			memcpy(edges_new, edges, sizeof(ui)*m);
			delete[] edges; edges = edges_new;
		}
		degen(K-1, peel_sequence, core, heap);
		assert(n == 0||core[peel_sequence[0]] >= K-1);
	}

	printf("*** Preprocessing time: %s\n", Utility::integer_to_string(t.elapsed()).c_str());
	fflush(stdout);

	if(m) {
		if(matrix == nullptr) matrix = new char[max_core*max_core];

		ept *pend = new ept[n+1];
		reorganize_adjacency_lists(n, peel_sequence, pend);

		// get_subgraph_statistics(pstart, pend, edges, vis);

		Pivotor_Matrix *pm = new Pivotor_Matrix();
		pm->init_C(max_core, K);

		ui *ids = new ui[max_core];

		memset(vis, 0, sizeof(char)*n);
		if(alg.rfind("exact", 0) == 0) {
			total_cnt = 0;
			for(ui i = 0;i < n;i ++) if(pend[peel_sequence[i]]-pstart[peel_sequence[i]] >= K-1) {
				ui u = peel_sequence[i];
				ui subgraph_n = pend[u] - pstart[u];
				if(subgraph_n == n - i - 1) {
					total_cnt += pm->count(subgraph_n+1, K);
					break;
				}

				ept number_of_edges = extract_subgraph(u, out_mapping, pend);
				if(number_of_edges*2/subgraph_n == subgraph_n-1) {
					total_cnt += pm->count(subgraph_n, K-1);
					continue;
				}

				double local_cnt = 0;
				if(alg == "exact") local_cnt = pm->pivotor(subgraph_n, matrix, degree, rid, ids, K-1);
				else if(alg == "exact1") local_cnt = pm->pivotor2(subgraph_n, matrix, degree, rid, ids, K-1);
				else if(alg == "exact2") local_cnt = pm->k_clique_counting(subgraph_n, matrix, degree, rid, ids, vis, K-1);
				total_cnt += local_cnt;
			}
			printf("*** Total number of %u-cliques: %.2lf\n", K, total_cnt);
		}
		else if(alg == "density") {
			total_cnt = total_ub = 0;
			if(dp == nullptr) dp = new double[max_core*K];
			ui *local_color = core;

			FILE *fout = Utility::open_file("statistics.csv", "w");
			fprintf(fout, "n, m, avg_degree, avg_density, exact_cnt, ub_cnt, exact/ub, time\n");
			for(ui i = 0;i < n;i ++) if(pend[peel_sequence[i]]-pstart[peel_sequence[i]] >= K-1) {
				Timer tt;
				ui u = peel_sequence[i];
				ui subgraph_n = pend[u] - pstart[u];
				double local_cnt = 0, local_ub = 0;
				if(subgraph_n == n-i-1) {
					local_cnt = pm->count(subgraph_n+1, K);
					total_cnt += local_cnt;
					total_ub += local_cnt;
					break;
				}
				ept number_of_edges = extract_subgraph(u, out_mapping, pend);
				if(number_of_edges*2/subgraph_n == subgraph_n-1) {
					local_cnt = pm->count(subgraph_n, K-1);
					total_cnt += local_cnt;
					total_ub += local_cnt;
					continue;
				}

				ui ids_n = subgraph_n;
				for(ui j = 0;j < subgraph_n;j ++) peel_sequence[j] = j;
				if(matrix_degen(subgraph_n, ids, ids_n, K-2, number_of_edges)) {
					local_cnt = pm->count(ids_n, K-1);
					total_cnt += local_cnt;
					total_ub += local_cnt;
					continue;
				}

				if(matrix_coloring(subgraph_n, ids, ids_n, local_color) < K-1) continue;

				reorder_peel_sequence_based_on_color(ids_n, ids, local_color, ids);
				local_ub = construct_DP_table_color_path(subgraph_n, ids_n, ids, K-1);
				total_ub += local_ub;

				local_cnt = pm->pivotor(subgraph_n, matrix, degree, rid, ids, K-1);
				total_cnt += local_cnt;

				assert(local_cnt <= local_ub);
				fprintf(fout, "%u, %u, %.5lf, %.5lf, %.2lf, %.2lf, %.5lf, %llu\n", subgraph_n, number_of_edges, number_of_edges*2.0/subgraph_n, number_of_edges*2.0/subgraph_n/(subgraph_n-1), local_cnt, local_ub, local_cnt/local_ub, tt.elapsed());
			}
			printf("Exact: %.2lf, UB: %.2lf, Density: %.8lf\n", total_cnt, total_ub, total_cnt/total_ub);
			fclose(fout);
		}
		else if(alg.rfind("TKDE23", 0) == 0) {
			if(dp == nullptr) dp = new double[max_core*K];
			ui *local_color = core;
			vector<double> ub(n);
			vector<ui> sample_nodes;
			total_cnt = total_ub = 0;
			for(ui i = 0;i < n;i ++) if(pend[peel_sequence[i]]-pstart[peel_sequence[i]] >= K-1) {
				ui u = peel_sequence[i];
				ui subgraph_n = pend[u] - pstart[u];
				if(subgraph_n == n-i-1) {
					total_cnt += pm->count(subgraph_n+1, K);
					break;
				}
				ept number_of_edges = extract_subgraph(u, out_mapping, pend);
				if(number_of_edges*2/subgraph_n == subgraph_n-1) {
					total_cnt += pm->count(subgraph_n, K-1);
					continue;
				}

				if(number_of_edges > (K-1)*subgraph_n) {
					ui ids_n = subgraph_n;
					for(ui j = 0;j < ids_n;j ++) ids[j] = j;
					if(matrix_degen(subgraph_n, ids, ids_n, K-2, number_of_edges)) {
						total_cnt += pm->count(ids_n, K-1);
						continue;
					}

					if(matrix_coloring(subgraph_n, ids, ids_n, local_color) < K-1) continue;

					sample_nodes.push_back(u);
					reorder_peel_sequence_based_on_color(ids_n, ids, local_color, ids);
					ub[u] = construct_DP_table_color_path(subgraph_n, ids_n, ids, K-1);
					total_ub += ub[u];
				}
				else total_cnt += pm->pivotor(subgraph_n, matrix, degree, rid, ids, K-1);
			}

#ifndef NDEBUG
			for(ui i = 0;i < n;i ++) assert(!vis[i]);
#endif

			printf("*** Preprocessing2 time: %s\n", Utility::integer_to_string(t.elapsed()).c_str());

			bool opt = false;
			if(alg == "TKDE23_opt") {
				opt = true;
				if(prob == nullptr) prob = new double[max_core*max_core*K];
				if(alias == nullptr) alias = new ui[max_core*max_core*K];
				if(neighbors == nullptr) neighbors = new ui[max_core*max_core];
			}

			double gamma = 15.22274/epsilon/epsilon; // gamma = 4(e-2)log(2/delta)/epsilon^2; delta=0.01
			double gamma1 = 1+(1+epsilon)*gamma;
			ul success_trails = 0;
			ul total_trails = 0;
			ul sample_size = (ul)gamma1;
			while(success_trails < gamma1&&!sample_nodes.empty()) {
				ul n_cliques, n_samples;
				printf("expected sample size: %llu\n", sample_size);
				fflush(stdout);
				if(opt) sample_a_batch_opt(sample_nodes, ub.data(), sample_size, n_cliques, n_samples, ids, out_mapping, local_color, pend, K);
				else sample_a_batch(sample_nodes, ub.data(), sample_size, n_cliques, n_samples, ids, out_mapping, local_color, pend, opt, K);
				printf("actual sample size: %llu, success trails: %llu, success rate: %.5lf\n", n_samples, n_cliques, double(n_cliques)/n_samples);
				success_trails += n_cliques;
				total_trails += n_samples;
				if(n_cliques > 0) sample_size = (gamma1/double(success_trails))*double(total_trails);
				else sample_size *= 10;
			}

			double estimate = 0;
			if(total_trails != 0) estimate = total_ub/total_trails*success_trails;
			printf("exact_cnt: %.2lf, estimate_cnt: %.2lf (UB: %.2lf), total_estimate_cnt: %.2lf\n", total_cnt, estimate, total_ub, total_cnt+estimate);
		}
		else if(alg.rfind("my", 0) == 0) {
			if(prob == nullptr) prob = new double[max_core*max_core*K];
			if(alias == nullptr) alias = new ui[max_core*max_core*K];
			if(neighbors == nullptr) neighbors = new ui[max_core*max_core];
			if(dp == nullptr) dp = new double[max_core*K];
			if(order_by_color == nullptr) order_by_color = new ui[max_core];
			ui *local_color = core;

			double a[] = {0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.1};
			vector<double> density_thresholds(a, a+7);
			// for(ui i = 0;i < density_thresholds.size();i ++) printf(" %.4lf", density_thresholds[i]);
			// printf("\n");
			if(heaps == nullptr) heaps = new vector<ui>[density_thresholds.size()]();
			if(ubs == nullptr) ubs = new vector<double>[density_thresholds.size()]();

			vector<ui> ids2(max_core);
			ul total_sample_time = 0, total_sample_n = 0;
			bool opt = false;
			if(alg == "my_opt") opt = true;
			total_cnt = total_ub = total_exp = 0;
			for(ui i = 0;i < n;i ++) if(pend[peel_sequence[i]]-pstart[peel_sequence[i]] >= K-1) {
				ui u = peel_sequence[i];
				ui subgraph_n = pend[u] - pstart[u];
				if(subgraph_n == n-i-1) {
					total_cnt += pm->count(subgraph_n+1, K);
					break;
				}
				ept number_of_edges = extract_subgraph(u, out_mapping, pend);
				if(number_of_edges*2/subgraph_n == subgraph_n-1) {
					total_cnt += pm->count(subgraph_n, K-1);
					continue;
				}

				ui ids_n = subgraph_n;
				for(ui j = 0;j < ids_n;j ++) ids[j] = j;
				process_a_subgraph(u, subgraph_n, ids_n, ids, out_mapping, K-1, pm, local_color, true,
						total_sample_n, total_sample_time, 0, density_thresholds, opt, ids2.data(), 0.0);
			}

			unordered_set<ui> *hash_tables = new unordered_set<ui>[n]();
			bool terminate = false;
			for(ui i = 0;i < density_thresholds.size()&&!terminate;i ++) {
				if(i) heaps[i-1].shrink_to_fit();
				if(heaps[i].empty()) continue;
				double epsilon_new = (1+total_cnt/total_ub)*epsilon;
				double gamma = 15.22274/epsilon_new/epsilon_new; // gamma = 4(e-2)log(2/delta)/epsilon^2; delta=0.01
				double gamma1 = 1+(1+epsilon_new)*gamma;
				double density_threshold = density_thresholds[i];
				if(i) density_threshold = density_thresholds[i-1];
				// printf("%.5lf %.5lf\n", total_exp, total_ub);
				//double expected_sampling_time = (1+epsilon_new)*gamma1*total_sample_time/density_threshold/total_sample_n;
				double expected_sampling_time = (1+epsilon_new)*gamma1/total_sample_n*total_sample_time/total_exp*total_ub;
				printf("dens_threshold: %.3lf (i=%u) %.3lf, eps_new: %.5lf, e_sample_time: %.0lf(s)\n", density_threshold, i, double(total_exp)/total_ub, epsilon_new, expected_sampling_time/1000000);

				vector<ui> &hh = heaps[i];
				ui cc = 0;
				while(!hh.empty()) {
					if(cc%1 == 0&&(i+1 == density_thresholds.size()||expected_sampling_time*2 < t.elapsed())) {
						printf("density is above %.5lf (i = %u)\n", density_threshold, i);
						terminate = true;
						break;
					}

					++ cc;

					ui ids_n = hh.back()-4; hh.pop_back();
					ui l = hh.back(), j, k;
					assert(hh.size() >= ids_n+3);
					for(j = 0, k = hh.size()-ids_n-3;j < ids_n;j ++, k++) ids[j] = hh[k];
					total_ub -= ubs[i].back();
					if(hh[k+1]) total_exp -= ubs[i].back()/hh[k+1]*hh[k];
					hh.resize(hh.size()-ids_n-3);
					ubs[i].pop_back();

					ui subgraph_n, peel_sequence_n;
					ept number_of_edges;
					if(ids_n == 1) {
						ui u = ids[0];
						subgraph_n = pend[u] - pstart[u];
						number_of_edges = extract_subgraph(u, out_mapping, pend);
						for(ui j = 0;j < subgraph_n;j ++) ids[j] = j;
						ids_n = subgraph_n;
						matrix_degen(subgraph_n, ids, ids_n, l-1, number_of_edges);
					}
					else {
						subgraph_n = ids_n;
						number_of_edges = get_subgraph(ids_n, ids, hash_tables, out_mapping, pend);
						for(ui j = 0;j < subgraph_n;j ++) ids[j] = j;
					}

					double density = double(number_of_edges)*2/subgraph_n/(subgraph_n-1);

					assert(l >= 3);
#ifndef NDEBUG
					assert_unique_and_bound("line 456", ids_n, ids, subgraph_n);
#endif
					for(ui j = 0;j < ids_n;j ++) {
						get_out_neighbors(matrix+ids[j]*subgraph_n, ids_n-j-1, ids+j+1, peel_sequence_n, peel_sequence);
						if(peel_sequence_n == ids_n-j-1) {
							total_cnt += pm->count(peel_sequence_n+1, l);
							break;
						}

#ifndef NDEBUG
						assert_unique_and_bound("line 459", ids_n, ids, subgraph_n);
						assert_unique_and_bound("line 460", peel_sequence_n, peel_sequence, subgraph_n);
#endif
						process_a_subgraph(n, subgraph_n, peel_sequence_n, peel_sequence, out_mapping, l-1, pm, local_color,
								false, total_sample_n, total_sample_time, i, density_thresholds, opt, ids2.data(), density);
					}
				}
			}

			printf("*** Preprocessing2 time: %s\n", Utility::integer_to_string(t.elapsed()).c_str());
			fflush(stdout);

			double estimate = 0;
			for(ui i = 0;i < density_thresholds.size();i ++) if(!heaps[i].empty()) {
				double epsilon_new = (1+total_cnt/total_ub)*epsilon;
				printf("pre total_ub: %.2lf, pre total_exp: %.2lf\n", total_ub, total_exp);

				ui nn = 0;
				for(ui j = i;j < density_thresholds.size();j ++) nn += ubs[j].size();
				//assert(nn == nodes_n);
				double *t_prob = new double[nn];
				ui *t_alias = new ui[nn];
				ui *lists = new ui[nn];
				ui cc = 0;
				double inferred_success = 0;
				for(ui j = i;j < density_thresholds.size();j ++) {
					vector<double> &ub = ubs[j];
					vector<ui> &hh = heaps[j];
					ul end = hh.size();
					for(ui k = ub.size();k > 0;k --) {
						t_prob[cc++] = ub[k-1];
						inferred_success += ub[k-1]*hh[end-4]/hh[end-3];
						assert(hh[end-1] <= end);
						end -= hh[end-1];
					}
					assert(end = 0);
				}
				total_ub = build_alias(nn, t_prob, t_alias, lists);
				printf("aft total_ub: %.2lf, aft total_exp: %.2lf\n", total_ub, total_exp);

				double gamma = 15.22274/epsilon_new/epsilon_new; // gamma = 4(e-2)log(2/delta)/epsilon^2; delta=0.01
				double gamma1 = 1+(1+epsilon_new)*gamma;
				double density_threshold = density_thresholds[i];
				if(i) density_threshold = density_thresholds[i-1];
				//if(opt) {
					double t_density = inferred_success/total_ub/1.1;
					if(t_density > density_threshold) density_threshold = t_density;
				//}
				ul sample_size = (1+epsilon_new)*gamma1/density_threshold;
				ul success_trails = 0;
				ul total_trails = 0;

				printf("nn: %u, density_threshold: %.5lf (i=%u), expected sample size: %llu, required success trails: %.2lf\n", nn, density_threshold, i, sample_size, gamma1);

				bool first_time = true;
				while(success_trails < gamma1) {
					for(ui j = 0;j < nn;j ++) lists[j] = 0;
					std::uniform_real_distribution<double> unif(0, nn);
					for(ui j = 0;j < sample_size;j ++) {
						double r = unif(gen);
						ui idx = ui(r);
						if(r - idx > t_prob[idx]+1e-10) idx = t_alias[idx];
						if(idx >= nn) idx = nn-1;
						assert(idx < nn);
						++ lists[idx];
					}
#ifndef NDEBUG
					for(ui j = nn;j > 0;j --) if(lists[j-1]) {
						printf(" %u:%u\n", j-1, lists[j-1]);
						break;
					}
					//for(ui j = 0;j < nn;j ++) if(lists[j]) printf(" %u:%u\n", j, lists[j]);
					//printf("\n");
#endif
					cc = 0;
					ul n_cliques = 0;
					for(ui j = i;j < density_thresholds.size();j ++) {
						vector<ui> &hh = heaps[j];
						ui end = hh.size();
#ifndef NDEBUG
						ui end_ub = ubs[j].size();
#endif
						while(end) {
							ui subgraph_n = hh[end-1]-4, l = hh[end-2];
							ui st = hh[end-3], ss = hh[end-4];
							if(first_time&&st <= lists[cc]) {
								lists[cc] -= st;
								n_cliques += ss;
							}
							if(lists[cc]) {
								for(ui k1 = 0, k = end-hh[end-1];k1 < subgraph_n;k1 ++, k++) ids[k1] = hh[k];

								ui peel_sequence_n;
								if(subgraph_n == 1) {
									ui u = ids[0];
									subgraph_n = pend[u] - pstart[u];
									ept number_of_edges = extract_subgraph(u, out_mapping, pend);
									for(ui k = 0;k < subgraph_n;k ++) peel_sequence[k] = k;
									peel_sequence_n = subgraph_n;
									matrix_degen(subgraph_n, peel_sequence, peel_sequence_n, l-1, number_of_edges);
								}
								else {
									get_subgraph(subgraph_n, ids, hash_tables, out_mapping, pend);
									for(ui k = 0;k < subgraph_n;k ++) peel_sequence[k] = k;
									peel_sequence_n = subgraph_n;
								}

								ui color_n = matrix_coloring(subgraph_n, peel_sequence, peel_sequence_n, local_color);
								assert(color_n >= l);
								reorder_peel_sequence_based_on_color(peel_sequence_n, peel_sequence, local_color, order_by_color);
								double t_ub = construct_DP_table_color_path_and_alias(subgraph_n, peel_sequence_n, order_by_color, l, local_color);
#ifndef NDEBUG
								if((ui)t_ub == 110&&(ui)ubs[j][end_ub-1] == 224) {
									print_array("ids", ids, subgraph_n);
									print_array("order_by_color", order_by_color, peel_sequence_n, out_mapping);
									print_dp(peel_sequence_n, l);
								}
								if(t_ub/ubs[j][end_ub-1] > 1+1e-5 || ubs[j][end_ub-1]/t_ub > 1+1e-5) {
									printf("t_ub: %.2lf, ubs[j][end_ub-1]: %.2lf\n", t_ub, ubs[j][end_ub-1]);
									//print_matrix(subgraph_n, matrix);
								}
								assert(t_ub/ubs[j][end_ub-1] < 1+1e-5&&ubs[j][end_ub-1]/t_ub < 1+1e-5);
#endif

								for(ui k = 0;k < lists[cc];k ++) n_cliques += sample_once_opt(subgraph_n, peel_sequence_n, order_by_color, local_color, l);
							}

							assert(hh[end-1] <= end);
							end -= hh[end-1];
							++ cc;
#ifndef NDEBUG
							-- end_ub;
							if(end == 0) assert(end_ub == 0);
#endif
						}
					}

					printf("a_sample_size: %llu, success: %llu, success rate: %.5lf\n", sample_size, n_cliques, double(n_cliques)/sample_size);
					success_trails += n_cliques;
					total_trails += sample_size;
					first_time = false;
				}
				delete[] lists;
				delete[] t_prob;
				delete[] t_alias;

				estimate = total_ub/total_trails*success_trails;
				break;
			}
			printf("exact_cnt: %.2lf, estimate_cnt: %.2lf (UB: %.2lf), total_estimate_cnt: %.2lf\n", total_cnt, estimate, total_ub, total_cnt+estimate);

			delete[] hash_tables;
		}
		else printf("The alg \"%s\" does not exists!\n", alg.c_str());

		delete[] ids;
		delete[] pend;
	}

	printf("*** Total time: %s\n", Utility::integer_to_string(t.elapsed()).c_str());

	delete heap;
	delete[] out_mapping;
	delete[] core;
	delete[] peel_sequence;
}

// compute core number and degeneracy ordering for the degree_threshold-core
ui Graph::degen(ui degree_threshold, ui *peel_sequence, ui *core, ListLinearHeap *heap) {
	for(ui i = 0;i < n;i ++) degree[i] = pstart[i+1] - pstart[i];

	ui queue_n = 0, new_size = 0;
	for(ui i = 0;i < n;i ++) if(degree[i] < degree_threshold) peel_sequence[queue_n ++] = i;
	for(ui i = 0;i < queue_n;i ++) {
		ui u = peel_sequence[i]; degree[u] = 0;
		for(ept j = pstart[u];j < pstart[u+1];j ++) if(degree[edges[j]] > 0) {
			if((degree[edges[j]] --) == degree_threshold) peel_sequence[queue_n ++] = edges[j];
		}
	}

	memset(vis, 0, sizeof(char)*n);
	for(ui i = 0;i < n;i ++) {
		if(degree[i] >= degree_threshold) peel_sequence[queue_n + (new_size ++)] = i;
		else {
			vis[i] = 1;
			core[i] = 0;
		}
	}
	assert(queue_n + new_size == n);

	ui max_core = 0;
	if(new_size != 0) {
		heap->init(new_size, new_size-1, peel_sequence+queue_n, degree);
		for(ui i = 0;i < new_size;i ++) {
			ui u, key;
			heap->pop_min(u, key);
			if(key > max_core) max_core = key;
			core[u] = max_core;
			peel_sequence[queue_n + i] = u;

			vis[u] = 1;

			for(ept j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j]] == 0) heap->decrement(edges[j], 1);
		}
	}

	printf("*** max_core: %u\n", max_core);
	return max_core;
}

// remove all vertices with core number < core_threshold
void Graph::core_shrink_graph(ui core_threshold, ui &n, ept &m, ui *peel_sequence, ui *core, ept *&pstart, ui *&edges) {
	ui cnt = 0;
	for(ui i = 0;i < n;i ++) if(core[i] >= core_threshold) rid[i] = cnt ++;

	if(cnt != n) {
		cnt = 0;
		ept pos = 0;
		for(ui i = 0;i < n;i ++) if(core[i] >= core_threshold) {
			ept t_start = pstart[i];
			pstart[cnt] = pos;
			for(ept j = t_start;j < pstart[i+1];j ++) if(core[edges[j]] >= core_threshold) edges[pos ++] = rid[edges[j]];
			++ cnt;
		}
		pstart[cnt] = pos;

		//assert(core[peel_sequence[n-cnt-1]] < core_threshold);
		//assert(cnt == 0||core[peel_sequence[n-cnt]] >= core_threshold);
		for(ui i = 0;i < cnt;i ++) peel_sequence[i] = rid[peel_sequence[n-cnt+i]];

		if(pos > 0&&pos < m/2) {
			ept *pstart_new = new ept[cnt+1];
			ui *edges_new = new ui[pos];
			memcpy(pstart_new, pstart, sizeof(ept)*(cnt+1));
			memcpy(edges_new, edges, sizeof(ui)*pos);
			delete[] pstart; pstart = pstart_new;
			delete[] edges; edges = edges_new;
		}

		n = cnt;
		m = pos;
	}

	printf("*** After core shrink: n = %s, m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());
}

void Graph::reorganize_adjacency_lists(ui n, ui *peel_sequence, ui *pend) {
	for(ui i = 0;i < n;i ++) rid[peel_sequence[i]] = i;
	for(ui i = 0;i < n;i ++) {
		ui &end = pend[i] = pstart[i];
		for(ui j = pstart[i];j < pstart[i+1];j ++) if(rid[edges[j]] > rid[i]) edges[end ++] = edges[j];
	}
	for(ui i = n;i > 0;i --) {
		ui u = peel_sequence[i-1];
		for(ui j = pstart[u];j < pend[u]&&rid[edges[j]] > rid[u];j ++) {
			ui v = edges[j];
			edges[pend[v] ++] = u;
			assert(pend[v] <= pstart[v+1]);
		}
	}
#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(pend[i] == pstart[i+1]);
#endif
	for(ui i = 0;i < n;i ++) {
		ui &end = pend[i] = pstart[i];
		while(end < pstart[i+1]&&rid[edges[end]] > rid[i]) ++ end;
	}
}

ept Graph::get_number_of_edges(ui u, char *exists, ept *pstart, ept *pend, ui *edges) {
	for(ept i = pstart[u];i < pend[u];i ++) exists[edges[i]] = 1;
	ept number_of_edges = 0;
	for(ept i = pstart[u];i < pend[u];i ++) {
		ui v = edges[i];
		for(ept j = pstart[v];j < pend[v];j ++) if(exists[edges[j]]) ++ number_of_edges;
	}
	for(ept i = pstart[u];i < pend[u];i ++) exists[edges[i]] = 0;
	return number_of_edges;
}

void Graph::get_subgraph_statistics(ept *pstart, ept *pend, ui *edges, ui K, char *vis) {
	FILE *fout = Utility::open_file("subgraphs.txt", "w");
	double min_density = 1, total_density = 0;
	ui cnt = 0, cntt[10];
	for(ui i = 0;i < 10;i ++) cntt[i] = 0;
	memset(vis, 0, sizeof(char)*n);
	for(ui i = 0;i < n;i ++) if(pend[i]-pstart[i] >= K-1) {
		ui number_of_vertices = pend[i]-pstart[i];
		ui number_of_edges = get_number_of_edges(i, vis, pstart, pend, edges);
		fprintf(fout, "%u %u %.5lf %.5lf\n", number_of_vertices, number_of_edges, number_of_edges*2.0/number_of_vertices, number_of_edges*2.0/number_of_vertices/(number_of_vertices-1));
		double density = number_of_edges*2.0/number_of_vertices/(number_of_vertices-1);
		++ cnt;
		total_density += density;
		if(density < min_density) min_density = density;
		for(ui j = 1;j < 10;j ++) if(density*10 <= j) ++ cntt[j];
	}
	printf("min_density: %.5lf, ave_density: %.5lf\n", min_density, total_density/cnt);
	for(ui i = 1;i < 10;i ++) printf("density <= 0.%u: %.3lf%%\n", i, cntt[i]*100.0/cnt);
	fclose(fout);
}

ept Graph::extract_subgraph(ui u, ui *out_mapping, ept *pend) {
	for(ept i = pstart[u];i < pend[u];i ++) {
		vis[edges[i]] = 1;
		rid[edges[i]] = i-pstart[u];
		out_mapping[i-pstart[u]] = edges[i];
	}
	ui tn = pend[u] - pstart[u];
	memset(matrix, 0, sizeof(char)*tn*tn);
	ept number_of_edges = 0;
	for(ept i = pstart[u];i < pend[u];i ++) {
		ui v = edges[i];
		for(ept j = pstart[v];j < pend[v];j ++) if(vis[edges[j]]) {
			++ number_of_edges;
			assert(v != edges[j]);
			assert(rid[v] != rid[edges[j]]);
			matrix[rid[v]*tn + rid[edges[j]]] = 1;
			matrix[rid[edges[j]]*tn + rid[v]] = 1;
		}
	}
	for(ept i = pstart[u];i < pend[u];i ++) vis[edges[i]] = 0;
	return number_of_edges;
}

bool Graph::matrix_degen(ui subgraph_n, ui *peel_sequence, ui &peel_sequence_n, ui threshold, ept &number_of_edges) {
#ifndef NDEBUG
	for(ui i = 0;i < peel_sequence_n;i ++) assert(peel_sequence[i] < subgraph_n);
#endif
	for(ui i = 0;i < peel_sequence_n;i ++) degree[i] = 0;
	for(ui i = 0;i < peel_sequence_n;i ++) {
		char *t_matrix = matrix + peel_sequence[i]*subgraph_n;
		for(ui j = i+1;j < peel_sequence_n;j ++) if(t_matrix[peel_sequence[j]]) {
			++ degree[i];
			++ degree[j];
		}
	}
	memset(vis, 0, sizeof(char)*peel_sequence_n);
	ui max_core = 0, ids_n = 0;
	// print_matrix(subgraph_n, matrix, peel_sequence_n, peel_sequence);
	ui *ids = rid;
	for(ui i = 0;i < peel_sequence_n;i ++) {
		ui min_degree = peel_sequence_n, idx;
		for(ui j = 0;j < peel_sequence_n;j ++) if(!vis[j]&&degree[j] < min_degree) {
			min_degree = degree[j];
			idx = j;
		}
		assert(min_degree < peel_sequence_n);

		if(ids_n == 0&&min_degree == peel_sequence_n-i-1) {
			for(ui j = 0;j < peel_sequence_n;j ++) if(!vis[j]) peel_sequence[ids_n++] = peel_sequence[j];
			memset(vis, 0, sizeof(char)*peel_sequence_n);
			peel_sequence_n = ids_n;
			return 1;
		}

		if(ids_n == 0) {
			number_of_edges = 0;
			for(ui j = 0;j < peel_sequence_n;j ++) if(!vis[j]) number_of_edges += degree[j];
			number_of_edges /= 2;
		}

		if(min_degree > max_core) max_core = min_degree;
		if(max_core >= threshold) ids[ids_n++] = peel_sequence[idx];
		vis[idx] = 1;
		char *t_matrix = matrix + peel_sequence[idx]*subgraph_n;
		for(ui j = 0;j < peel_sequence_n;j ++) if(!vis[j]&&t_matrix[peel_sequence[j]]) {
			assert(degree[j]);
			-- degree[j];
		}
	}
	memset(vis, 0, sizeof(char)*peel_sequence_n);
	for(ui i = 0;i < ids_n;i ++) peel_sequence[i] = ids[i];
	peel_sequence_n = ids_n;
	return 0;
}

ui Graph::matrix_coloring(ui subgraph_n, ui *peel_sequence, ui peel_sequence_n, ui *local_color) {
	ui max_color = 0;
	memset(local_color, 0, sizeof(ui)*subgraph_n);
	for(ui i = peel_sequence_n;i > 0;i --) {
		ui u = peel_sequence[i-1];
		char *t_matrix = matrix + u*subgraph_n;
		for(ui j = i;j < peel_sequence_n;j ++) if(t_matrix[peel_sequence[j]]) vis[local_color[peel_sequence[j]]] = 1;
		for(ui j = 0;;j ++) if(!vis[j]) {
			local_color[u] = j;
			if(j > max_color) max_color = j;
			break;
		}
		for(ui j = i;j < peel_sequence_n;j ++) if(t_matrix[peel_sequence[j]]) vis[local_color[peel_sequence[j]]] = 0;
	}
	return max_color+1;
}

void Graph::reorder_peel_sequence_based_on_color(ui peel_sequence_n, ui *peel_sequence, ui *local_color, ui *order_by_color) {
	ui *cnt = degree;
	memset(cnt, 0, sizeof(ui)*peel_sequence_n);
	for(ui i = 0;i < peel_sequence_n;i ++) ++ cnt[local_color[peel_sequence[i]]];
	for(ui i = 1;i < peel_sequence_n;i ++) cnt[i] += cnt[i-1];

	for(ui i = 0;i < peel_sequence_n;i ++) rid[i] = cnt[local_color[peel_sequence[i]]] --;
	for(ui i = 0;i < peel_sequence_n;i ++) cnt[peel_sequence_n-rid[i]] = peel_sequence[i];
	for(ui i = 0;i < peel_sequence_n;i ++) order_by_color[i] = cnt[i];
}

double Graph::construct_DP_table_color_path(ui subgraph_n, ui peel_sequence_n, ui *peel_sequence, ui K) {
	for(ui i = peel_sequence_n;i > 0;i --) {
		double *t_dp = dp + (i-1)*K;
		t_dp[0] = 1;
		for(ui k = 1;k < K;k ++) t_dp[k] = 0;
		char *t_matrix = matrix + peel_sequence[i-1]*subgraph_n;
		for(ui j = i;j < peel_sequence_n;j ++) if(t_matrix[peel_sequence[j]]) {
			for(ui k = 1;k < K;k ++) t_dp[k] += dp[j*K + k-1];
		}
	}

	double local_ub = 0;
	for(ui j = 0;j < peel_sequence_n;j ++) local_ub += dp[j*K + K-1];

	return local_ub;
}

double Graph::construct_DP_table_color_path_and_alias(ui subgraph_n, ui peel_sequence_n, ui *peel_sequence, ui K, ui *lists) {
	for(ui i = peel_sequence_n;i > 0;i --) {
		double *t_dp = dp + (i-1)*K;
		t_dp[0] = 1; // choose 1 vertex
		for(ui k = 1;k < K;k ++) t_dp[k] = 0;
		char *t_matrix = matrix + peel_sequence[i-1]*subgraph_n;
		for(ui j = i;j < peel_sequence_n;j ++) if(t_matrix[peel_sequence[j]]) {
			for(ui k = 1;k < K;k ++) t_dp[k] += dp[j*K + k-1];
		}
	}

	double local_ub = 0;
	for(ui i = 0;i < peel_sequence_n;i ++) {
		prob[i] = dp[i*K+ (K-1)];
		local_ub += prob[i];
	}
	build_alias(peel_sequence_n, prob, alias, lists);

	for(ui i = 0;i < peel_sequence_n;i ++) {
		char *t_matrix = matrix + peel_sequence[i]*subgraph_n;
		ui *t_neighbors = neighbors + i*peel_sequence_n;
		degree[i] = 0;
		for(ui j = i+1;j < peel_sequence_n;j ++) if(t_matrix[peel_sequence[j]]) t_neighbors[degree[i]++] = j;

		for(ui k = 1;k < K;k ++) {
			double *t_prob = prob+(i*K+k)*peel_sequence_n;
			ui *t_alias = alias+(i*K+k)*peel_sequence_n;
			for(ui j = 0;j < degree[i];j ++) t_prob[j] = dp[t_neighbors[j]*K+k-1];
			build_alias(degree[i], t_prob, t_alias, lists);
		}
	}

	return local_ub;
}

double Graph::build_alias(ui num, double *prob, ui *alias, ui *lists) {
	double total_prob = 0;
	for(ui i = 0;i < num;i ++) total_prob += prob[i];
	ui S_end = 0, L_start = num;
	for(ui i = 0;i < num;i ++) {
		prob[i] = prob[i]*num/total_prob;
		if(prob[i] < 1) lists[S_end++] = i;
		else lists[--L_start] = i;
	}
	while(S_end > 0&&L_start < num) {
		ui a = lists[--S_end], b = lists[L_start++];
		alias[a] = b;
		prob[b] = (prob[a]+prob[b])-1;
		if(prob[b] < 1) lists[S_end++] = b;
		else lists[--L_start] = b;
	}
	while(L_start < num) prob[lists[L_start++]] = 1;
	while(S_end > 0) prob[lists[--S_end]] = 1;

	return total_prob;
}

void Graph::print_matrix(ui n, char *matrix) {
	for(ui i = 0;i < n;i ++) {
		for(ui j = 0;j < n;j ++) {
			if(matrix[i*n+j]) printf(" 1");
			else printf(" 0");
		}
		printf("\n");
	}
}

void Graph::print_matrix(ui subgraph_n, char *matrix, ui ids_n, ui *ids) {
	for(ui i = 0;i < ids_n;i ++) {
		for(ui j = 0;j < ids_n;j ++) {
			if(matrix[ids[i]*subgraph_n+ids[j]]) printf(" 1");
			else printf(" 0");
		}
		printf("\n");
	}
}

void Graph::print_array(const char *str, ui *array, ui array_n) {
	printf("%s:", str);
	for(ui i = 0;i < array_n;i ++) printf(" %u", array[i]);
	printf("\n");
}

void Graph::print_array(const char *str, ui *array, ui array_n, ui *out_mapping) {
	printf("%s:", str);
	for(ui i = 0;i < array_n;i ++) printf(" %u", out_mapping[array[i]]);
	printf("\n");
}

void Graph::print_array(const char *str, char *array, ui array_n) {
	printf("%s:", str);
	for(ui i = 0;i < array_n;i ++) printf(" %u", array[i]);
	printf("\n");
}

void Graph::print_dp(ui peel_sequence_n, ui K) {
	printf("DP table\n");
	for(ui i = 0;i < peel_sequence_n;i ++) {
		for(ui j = 0;j < K;j ++) printf(" %.2lf", dp[i*K+j]);
		printf("\n");
	}
}

void Graph::assert_unique_and_bound(const char *str, ui ids_n, ui *ids, ui subgraph_n) {
	for(ui i = 0;i < ids_n;i ++) assert(ids[i] < subgraph_n&&!vis[ids[i]]);
	for(ui i = 0;i < ids_n;i ++) {
		if(vis[ids[i]]) print_array(str, ids, ids_n);
		assert(!vis[ids[i]]);
		vis[ids[i]] = 1;
	}
	for(ui i = 0;i < ids_n;i ++) vis[ids[i]] = 0;
}

void Graph::sample_a_batch(vector<ui> &nodes, double *ub, ul sample_size, ul &n_cliques, ul &n_samples,
		ui *peel_sequence, ui *out_mapping, ui *local_color, ept *pend, bool opt, ui K) {
	n_cliques = n_samples = 0;
#ifndef NDEBUG
	double tt = 0;
	for(ui i = 0;i < nodes.size();i ++) tt += ub[nodes[i]];
	//printf("tt: %.2lf, total_ub: %.2lf\n", tt, total_ub);
	assert(tt/total_ub < 1+1e-5&&total_ub/tt < 1+1e-5);
#endif
	//printf("sample_size: %.2lf, total_ub: %.2lf\n", sample_size, total_ub);
	for(ui i = 0;i < nodes.size();i ++) {
		ui u = nodes[i];
		ul samples = (ul)(sample_size/total_ub*ub[u] + 1e-6);
		if(samples == 0) continue;

		// printf("i:%u, samples: %llu\n", i, samples);

		ui subgraph_n = pend[u] - pstart[u];
		extract_subgraph(u, out_mapping, pend);
		ui peel_sequence_n = subgraph_n;
		for(ui j = 0;j < peel_sequence_n;j ++) peel_sequence[j] = j;
		ept number_of_edges;
		bool clique = matrix_degen(subgraph_n, peel_sequence, peel_sequence_n, K-2, number_of_edges);
		assert(!clique);
		ui color_n = matrix_coloring(subgraph_n, peel_sequence, peel_sequence_n, local_color);
		assert(color_n >= K-1);

		reorder_peel_sequence_based_on_color(peel_sequence_n, peel_sequence, local_color, peel_sequence);
		if(opt) construct_DP_table_color_path_and_alias(subgraph_n, peel_sequence_n, peel_sequence, K-1, local_color);
		else construct_DP_table_color_path(subgraph_n, peel_sequence_n, peel_sequence, K-1);

		for(ul j = 0;j < samples;j ++) {
			if(opt) n_cliques += sample_once_opt(subgraph_n, peel_sequence_n, peel_sequence, local_color, K-1);
			else n_cliques += sample_once(subgraph_n, peel_sequence_n, peel_sequence, local_color, K-1);
		}
		n_samples += samples;
	}
}

void Graph::sample_a_batch_opt(vector<ui> &nodes, double *ub, ul sample_size, ul &n_cliques, ul &n_samples,
		ui *peel_sequence, ui *out_mapping, ui *local_color, ept *pend, ui K) {
	n_cliques = n_samples = 0;
#ifndef NDEBUG
	double tt = 0;
	for(ui i = 0;i < nodes.size();i ++) tt += ub[nodes[i]];
	//printf("tt: %.2lf, total_ub: %.2lf\n", tt, total_ub);
	assert(tt/total_ub < 1+1e-5&&total_ub/tt < 1+1e-5);
#endif

	double *t_prob = new double[nodes.size()];
	ui *t_alias = new ui[nodes.size()];
	ui *lists = new ui[nodes.size()];
	double inferred_success = 0;
	for(ui i = 0;i < nodes.size();i ++) t_prob[i] = ub[nodes[i]];
	total_ub = build_alias(nodes.size(), t_prob, t_alias, lists);

	for(ui i = 0;i < nodes.size();i ++) lists[i] = 0;
	std::uniform_real_distribution<double> unif(0, nodes.size());
	for(ui j = 0;j < sample_size;j ++) {
		double r = unif(gen);
		ui idx = ui(r);
		if(r - idx > t_prob[idx]+1e-10) idx = t_alias[idx];
		if(idx >= nodes.size()) idx = nodes.size()-1;
		assert(idx < nn);
		++ lists[idx];
	}

	for(ui i = 0;i < nodes.size();i ++) {
		ui u = nodes[i];
		ul samples = lists[i];
		if(samples == 0) continue;

		// printf("i:%u, samples: %llu\n", i, samples);

		ui subgraph_n = pend[u] - pstart[u];
		extract_subgraph(u, out_mapping, pend);
		ui peel_sequence_n = subgraph_n;
		for(ui j = 0;j < peel_sequence_n;j ++) peel_sequence[j] = j;
		ept number_of_edges;
		bool clique = matrix_degen(subgraph_n, peel_sequence, peel_sequence_n, K-2, number_of_edges);
		assert(!clique);
		ui color_n = matrix_coloring(subgraph_n, peel_sequence, peel_sequence_n, local_color);
		assert(color_n >= K-1);

		reorder_peel_sequence_based_on_color(peel_sequence_n, peel_sequence, local_color, peel_sequence);
		construct_DP_table_color_path_and_alias(subgraph_n, peel_sequence_n, peel_sequence, K-1, local_color);

		for(ul j = 0;j < samples;j ++) {
			n_cliques += sample_once_opt(subgraph_n, peel_sequence_n, peel_sequence, local_color, K-1);
		}
		n_samples += samples;
	}

	delete[] t_prob;
	delete[] t_alias;
	delete[] lists;
}

ui Graph::sample_once(ui subgraph_n, ui peel_sequence_n, ui *peel_sequence, ui *cliques, ui K) {
	std::uniform_real_distribution<double> unif(0, 1);

	double sum = 0;
	for(ui i = 0;i < peel_sequence_n;i ++) sum += dp[i*K+K-1];
	ui pre_idx = peel_sequence_n-1;
	double r = unif(gen), partial_sum = 0;
	for(ui i = 0;i < peel_sequence_n;i ++) {
		partial_sum += dp[i*K+K-1];
		if(partial_sum/sum + 1e-6 > r) {
			pre_idx = i;
			break;
		}
	}
	cliques[0] = peel_sequence[pre_idx];

	for(ui i = 1;i < K;i ++) {
		char *t_matrix = matrix + peel_sequence[pre_idx]*subgraph_n;
		sum = 0;
		ui now_idx = 0;
		for(ui j = pre_idx+1;j < peel_sequence_n;j ++) if(t_matrix[peel_sequence[j]]) {
			sum += dp[j*K+(K-1-i)];
			now_idx = j;
		}
		r = unif(gen);
		partial_sum = 0;
		for(ui j = pre_idx+1;j < peel_sequence_n;j ++) if(t_matrix[peel_sequence[j]]) {
			partial_sum += dp[j*K+(K-1-i)];
			if(partial_sum/sum + 1e-6 > r) {
				now_idx = j;
				break;
			}
		}

		t_matrix = matrix + peel_sequence[now_idx]*subgraph_n;
		for(ui j = 0;j < i;j ++) if(!t_matrix[cliques[j]]) return 0;
		cliques[i] = peel_sequence[now_idx];
		pre_idx = now_idx;
	}

	return 1;
}

ui Graph::sample_once_opt(ui subgraph_n, ui peel_sequence_n, ui *peel_sequence, ui *cliques, ui K) {
	std::uniform_real_distribution<double> unif(0, 1);

	double r = unif(gen)*peel_sequence_n;
	ui pre_idx = ui(r);
	if(r - pre_idx > prob[pre_idx]+1e-10) pre_idx = alias[pre_idx];
	if(pre_idx >= peel_sequence_n) pre_idx = peel_sequence_n-1;
	cliques[0] = peel_sequence[pre_idx];

	for(ui i = 1;i < K;i ++) {
		double *t_prob = prob + (pre_idx*K+K-i)*peel_sequence_n;
		ui *t_alias = alias + (pre_idx*K+K-i)*peel_sequence_n;

		r = unif(gen)*degree[pre_idx];
		ui now_idx = ui(r);
		if(r - now_idx > t_prob[now_idx]+1e-10) now_idx = t_alias[now_idx];
		if(now_idx >= degree[pre_idx]) now_idx = degree[pre_idx]-1;
		now_idx = neighbors[pre_idx*peel_sequence_n+now_idx];

		char *t_matrix = matrix + peel_sequence[now_idx]*subgraph_n;
		for(ui j = 0;j < i;j ++) if(!t_matrix[cliques[j]]) return 0;
		cliques[i] = peel_sequence[now_idx];
		pre_idx = now_idx;
	}

	return 1;
}

void Graph::get_out_neighbors(char *adj, ui ids_n, ui *ids, ui &peel_sequence_n, ui *peel_sequence) {
	peel_sequence_n = 0;
	for(ui i = 0;i < ids_n;i ++) if(adj[ids[i]]) peel_sequence[peel_sequence_n++] = ids[i];
}

ept Graph::count_number_of_edges(ui subgraph_n, ui ids_n, ui *ids) {
	ept number_of_edges = 0;
	for(ui i = 0;i < ids_n;i ++) {
		char *t_matrix = matrix + ids[i]*subgraph_n;
		for(ui j = i+1;j < ids_n;j ++) if(t_matrix[ids[j]]) ++ number_of_edges;
	}
	return number_of_edges;
}

void Graph::push_to_heap(vector<ui> &hh, ui peel_sequence_n, ui *peel_sequence, ui *out_mapping, ui success, ui trails, ui l) {
	for(ui i = 0;i < peel_sequence_n;i ++) {
		if(out_mapping == nullptr) hh.push_back(peel_sequence[i]);
		else hh.push_back(out_mapping[peel_sequence[i]]);
	}
	hh.push_back(success);
	hh.push_back(trails);
	hh.push_back(l);
	hh.push_back(peel_sequence_n+4);
}

void Graph::sample_a_batch_for_one_subgraph(ui subgraph_n, ui peel_sequence_n, ui *order_by_color, ui *clique, ui l, ul &success, ul &total) {
	ui beta = 10;
	if(peel_sequence_n/l < beta) beta = peel_sequence_n/l;
	total = peel_sequence_n*beta;
	if(total == 0) printf("WA %u %u\n", l, peel_sequence_n);
	success = 0;
	for(ul j = 0;j < total;j ++) success += sample_once_opt(subgraph_n, peel_sequence_n, order_by_color, clique, l);
}

ept Graph::get_subgraph(ui ids_n, ui *ids, unordered_set<ui> *hash_tables, ui *out_mapping, ept *pend) {
	ept number_of_edges;
	memset(matrix, 0, sizeof(char)*ids_n*ids_n);
	for(ui i = 0;i < ids_n;i ++) {
		ui u = ids[i];
		unordered_set<ui> &hash_table = hash_tables[u];
		if(pend[u] > pstart[u]&&hash_table.empty()) {
			for(ept j = pstart[u];j < pend[u];j ++) hash_table.insert(edges[j]);
		}

		out_mapping[i] = ids[i];
		for(ui j = 0;j < ids_n;j ++) if(hash_table.find(ids[j]) != hash_table.end()) {
			matrix[i*ids_n + j] = matrix[j*ids_n + i] = 1;
			++ number_of_edges;
		}
	}
#ifndef NDEBUG
	for(ui i = 0;i < ids_n;i ++) {
		assert(!matrix[i*ids_n+i]);
		for(ui j = i+1;j < ids_n;j ++) assert(matrix[i*ids_n+j] == matrix[j*ids_n+i]);
	}
#endif
	return number_of_edges;
}

void Graph::process_a_subgraph(ui leading_vertex, ui subgraph_n, ui peel_sequence_n, ui *peel_sequence, ui *out_mapping, ui l, Pivotor_Matrix *pm,
		ui *local_color, bool record_time, ul &total_sample_n, ul &total_sample_time, ui start_idx, vector<double> &density_thresholds, bool opt, ui *ids, double dd) {
#ifndef NDEBUG
	assert_unique_and_bound("line 1113", peel_sequence_n, peel_sequence, subgraph_n);
#endif

	if(peel_sequence_n < l) return ;
	if(l == 2) {
		total_cnt += count_number_of_edges(subgraph_n, peel_sequence_n, peel_sequence);
		return ;
	}
	assert(l > 2);
	if(peel_sequence_n <= l+10&&peel_sequence_n <= 2*l) {
		total_cnt += pm->pivotor(subgraph_n, matrix, degree, rid, peel_sequence, peel_sequence_n, l);
		return ;
	}
	ept number_of_edges;
	if(matrix_degen(subgraph_n, peel_sequence, peel_sequence_n, l-1, number_of_edges)) {
		total_cnt += pm->count(peel_sequence_n, l);
		return ;
	}
	if(peel_sequence_n <= l+10&&peel_sequence_n <= 2*l) {
		total_cnt += pm->pivotor(subgraph_n, matrix, degree, rid, peel_sequence, peel_sequence_n, l);
		return ;
	}

	double t_ub = 0;
	ul success=0, total=0;

	if(!opt||double(number_of_edges)*2/peel_sequence_n/(peel_sequence_n-1) > dd*density_factor) {
		if(matrix_coloring(subgraph_n, peel_sequence, peel_sequence_n, local_color) < l) return ;

		reorder_peel_sequence_based_on_color(peel_sequence_n, peel_sequence, local_color, order_by_color);
		t_ub = construct_DP_table_color_path_and_alias(subgraph_n, peel_sequence_n, order_by_color, l, local_color);
		total_ub += t_ub;

		if(record_time) {
			Timer tt;
			sample_a_batch_for_one_subgraph(subgraph_n, peel_sequence_n, order_by_color, local_color, l, success, total);
			total_sample_n += total;
			total_sample_time += tt.elapsed();
		}
		else sample_a_batch_for_one_subgraph(subgraph_n, peel_sequence_n, order_by_color, local_color, l, success, total);
		total_exp += t_ub/total*success;
		// printf("1 %.5lf\n", total_exp);
	}
	if(dd < 0.0001) dd = double(number_of_edges)*2/peel_sequence_n/(peel_sequence_n-1);

	double density = 0;
	if(total) density = double(success)/total;
	ul old_size = heaps[start_idx].size();
#ifndef NDEBUG
	assert_unique_and_bound("line 1153", peel_sequence_n, peel_sequence, subgraph_n);
	assert_unique_and_bound("line 1154", subgraph_n, out_mapping, n);
#endif
	for(ui i = start_idx;i < density_thresholds.size();i ++) if(density < density_thresholds[i]) {
		if(!opt&&leading_vertex < n) {
			peel_sequence[0] = leading_vertex;
			push_to_heap(heaps[i], 1, peel_sequence, nullptr, success, total, l);
		}
		else if(i == start_idx&&opt) push_to_heap(heaps[i], peel_sequence_n, peel_sequence, nullptr, success, total, l);
		else push_to_heap(heaps[i], peel_sequence_n, peel_sequence, out_mapping, success, total, l);
		ubs[i].push_back(t_ub);
#ifndef NDEBUG
		if((ui)t_ub == 224&&out_mapping[peel_sequence[0]] == 9009) {
			print_array("subgraph1:", peel_sequence, peel_sequence_n, out_mapping);
			print_array("order_by_color", order_by_color, peel_sequence_n, out_mapping);
			print_dp(peel_sequence_n, l);
			//print_matrix(subgraph_n, matrix, peel_sequence_n, peel_sequence);
		}
#endif
		break;
	}

	while(opt&&heaps[start_idx].size() > old_size) {
		vector<ui> &hh = heaps[start_idx];
		ui ids_n = hh.back()-4; hh.pop_back();
		ui l = hh.back(), j, k;
		for(j = 0, k = hh.size()-ids_n-3;j < ids_n;j ++, k++) ids[j] = hh[k];
		total_ub -= ubs[start_idx].back();
		if(hh[k+1]) total_exp -= ubs[start_idx].back()/hh[k+1]*hh[k];
		//printf("2 %.5lf %u %u\n", total_exp, hh[k], hh[k+1]);
		hh.resize(hh.size()-ids_n-3);
		ubs[start_idx].pop_back();

#ifndef NDEBUG
		assert_unique_and_bound("line 1180", ids_n, ids, subgraph_n);
#endif

		assert(l >= 3);
		for(ui j = 0;j < ids_n;j ++) {
			get_out_neighbors(matrix+ids[j]*subgraph_n, ids_n-j-1, ids+j+1, peel_sequence_n, peel_sequence);
			if(peel_sequence_n < l-1) continue;
			if(peel_sequence_n == ids_n-j-1) {
#ifndef NDEBUG
				for(ui k = 0;k < peel_sequence_n;k ++) for(ui k1 = k+1;k1 < peel_sequence_n;k1 ++) {
					assert(matrix[peel_sequence[k]*subgraph_n+peel_sequence[k1]]);
				}
#endif
				total_cnt += pm->count(peel_sequence_n+1, l);
				break;
			}
			if(l == 3) {
				total_cnt += count_number_of_edges(subgraph_n, peel_sequence_n, peel_sequence);
				continue;
			}
			assert(l > 3);
			if(peel_sequence_n <= l-1+10&&peel_sequence_n <= 2*l-2) {
				total_cnt += pm->pivotor(subgraph_n, matrix, degree, rid, peel_sequence, peel_sequence_n, l-1);
				continue;
			}
			if(matrix_degen(subgraph_n, peel_sequence, peel_sequence_n, l-2, number_of_edges)) {
				total_cnt += pm->count(peel_sequence_n, l-1);
				continue;
			}
			if(peel_sequence_n <= l-1+10&&peel_sequence_n <= 2*l-2) {
				total_cnt += pm->pivotor(subgraph_n, matrix, degree, rid, peel_sequence, peel_sequence_n, l-1);
				continue;
			}

			t_ub = 0;
			success = total = 0;
			if(!opt||double(number_of_edges)*2/peel_sequence_n/(peel_sequence_n-1) > dd*density_factor) {
				if(matrix_coloring(subgraph_n, peel_sequence, peel_sequence_n, local_color) < l-1) continue;

				reorder_peel_sequence_based_on_color(peel_sequence_n, peel_sequence, local_color, order_by_color);
				t_ub = construct_DP_table_color_path_and_alias(subgraph_n, peel_sequence_n, order_by_color, l-1, local_color);
				total_ub += t_ub;

				sample_a_batch_for_one_subgraph(subgraph_n, peel_sequence_n, order_by_color, local_color, l-1, success, total);
				total_exp += t_ub/total*success;
				// printf("1 %.5lf %u %u\n", total_exp, success, total);
			}

#ifndef NDEBUG
			assert_unique_and_bound("line 1223", peel_sequence_n, peel_sequence, subgraph_n);
			assert_unique_and_bound("line 1224", subgraph_n, out_mapping, n);
#endif
			double density = 0;
			if(total) density = double(success)/total;
			for(ui i = start_idx;i < density_thresholds.size();i ++) if(density < density_thresholds[i]) {
				if(i == start_idx) push_to_heap(heaps[i], peel_sequence_n, peel_sequence, nullptr, success, total, l-1);
				else push_to_heap(heaps[i], peel_sequence_n, peel_sequence, out_mapping, success, total, l-1);
#ifndef NDEBUG
				if((ui)t_ub == 224&&out_mapping[peel_sequence[0]] == 9009) {
					print_array("subgraph2:", peel_sequence, peel_sequence_n, out_mapping);
					print_array("order_by_color", order_by_color, peel_sequence_n, out_mapping);
					print_dp(peel_sequence_n, l-1);
					//print_matrix(subgraph_n, matrix, peel_sequence_n, peel_sequence);
				}
#endif
				ubs[i].push_back(t_ub);
				break;
			}
		}
	}
}
