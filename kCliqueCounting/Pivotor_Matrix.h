/*
 * Pivotor_Matrix.h
 *
 *  Created on: 26 Sep 2023
 *      Author: ljchang
 */

#ifndef PIVOTOR_MATRIX_H_
#define PIVOTOR_MATRIX_H_

#include "Utility.h"

class Pivotor_Matrix {
private:
	ui n;
	char *matrix;
	ui *degree;
	ui *rid;
	double cnt;
	std::vector<ui> branching_vertices_buf;

	double *C;
	ui K;

public:
	Pivotor_Matrix() {
		matrix = nullptr;
		degree = nullptr;
		rid = nullptr;
		cnt = 0;
		n = 0;

		C = nullptr;
		K = 0;
	}

	~Pivotor_Matrix() {
		if(C != nullptr) {
			delete[] C;
			C = nullptr;
		}
	}

	void init_C(ui degen, ui K_) {
		K = K_;
		C = new double[(degen+1)*(K+1)];
		C[0] = 1;
		for(ui i = 1;i <= degen;i ++) {
			C[i*(K+1)] = 1;
			if(i <= K) C[i*(K+1)+i] = 1;
			for(ui j = 1;j < i&&j <= K;j ++) C[i*(K+1)+j] = C[(i-1)*(K+1)+j] + C[(i-1)*(K+1)+j-1];
		}
	}

	double k_clique_counting(ui n_, char *matrix_, ui *degree_, ui *rid_, ui *ids, char *vis, ui K) {
		n = n_;
		matrix = matrix_;
		degree = degree_;
		rid = rid_;

		for(ui i = 0;i < n;i ++) ids[i] = rid[i] = i;

		cnt = 0;
		memset(vis, 0, sizeof(char)*n);
		recursive_kc(ids, n, n, K, 0, vis);

		return cnt;
	}

	double pivotor(ui n_, char *matrix_, ui *degree_, ui *rid_, ui *ids, ui K) {
		n = n_;
		matrix = matrix_;
		degree = degree_;
		rid = rid_;

		for(ui i = 0;i < n;i ++) ids[i] = rid[i] = i;
		for(ui i = 0;i < n;i ++) {
			degree[i] = 0;
			char *t_matrix = matrix + i*n;
			for(ui j = 0;j < n;j ++) if(t_matrix[j]) ++ degree[i];
		}

		cnt = 0;
		recursive(ids, n, n, K, 0, 0);

		return cnt;
	}

	double pivotor(ui n_, char *matrix_, ui *degree_, ui *rid_, ui *ids, ui ids_n, ui K) {
		n = n_;
		matrix = matrix_;
		degree = degree_;
		rid = rid_;

		for(ui i = 0;i < ids_n;i ++) rid[ids[i]] = i;
		for(ui i = 0;i < ids_n;i ++) {
			degree[ids[i]] = 0;
			char *t_matrix = matrix + ids[i]*n;
			for(ui j = 0;j < ids_n;j ++) if(t_matrix[ids[j]]) ++ degree[ids[i]];
		}

		cnt = 0;
		recursive(ids, ids_n, n, K, 0, 0);

		return cnt;
	}

	double pivotor2(ui n_, char *matrix_, ui *degree_, ui *rid_, ui *ids, ui K) {
		n = n_;
		matrix = matrix_;
		degree = degree_;
		rid = rid_;

		for(ui i = 0;i < n;i ++) ids[i] = rid[i] = i;
		for(ui i = 0;i < n;i ++) {
			degree[i] = 0;
			char *t_matrix = matrix + i*n;
			for(ui j = 0;j < n;j ++) if(t_matrix[j]) ++ degree[i];
		}

		cnt = 0;
		recursive2(ids, n, n, K, 0, 0);

		return cnt;
	}

	double count(ui a, ui b) {
		if(a < b) return 0;
		return C[a*(K+1)+b];
	}

private:
	void recursive_kc(ui *ids, ui ids_n, ui u, ui K, ui branching_vertices_start, char *vis) {
		if(u < n) {
			char *t_matrix = matrix + u*n;
			ui new_ids_n = ids_n;
			for(ui i = 0;i < new_ids_n;) {
				if(!t_matrix[ids[i]]) {
					-- new_ids_n;
					std::swap(rid[ids[i]], rid[ids[new_ids_n]]);
					std::swap(ids[i], ids[new_ids_n]);
				}
				else ++ i;
			}
			ids_n = new_ids_n;
		}

		if(K == 1) {
			cnt += ids_n;
			return ;
		}

		for(ui i = 0;i < ids_n;i ++) {
			degree[ids[i]] = 0;
			char *t_matrix = matrix + ids[i]*n;
			for(ui j = 0;j < ids_n;j ++) if(t_matrix[ids[j]]) ++ degree[ids[i]];
		}

		ui branching_vertices_size = 0;
		for(ui i = 0;i < ids_n&&ids_n-i >= K;i ++) {
			ui min_degree = n, v;
			for(ui j = 0;j < ids_n;j ++) if(!vis[ids[j]]&&degree[ids[j]] < min_degree) {
				min_degree = degree[ids[j]];
				v = ids[j];
			}
			if(min_degree == ids_n-i-1) {
				cnt += count(min_degree+1, K);
				break;
			}
			vis[v] = 1;
			add_to_branching_vertices(v, branching_vertices_start, branching_vertices_size);
			char *t_matrix = matrix + v*n;
			for(ui j = 0;j < ids_n;j ++) if(!vis[ids[j]]&&t_matrix[ids[j]]) -- degree[ids[j]];
		}
		for(ui i = 0;i < ids_n;i ++) vis[ids[i]] = 0;

		for(ui i = 0;i < branching_vertices_size;i ++) {
			ui v = branching_vertices_buf[branching_vertices_start+i];
			ui idx = rid[v];
			std::swap(rid[v], rid[ids[ids_n-1-i]]);
			std::swap(ids[idx], ids[ids_n-1-i]);

			recursive_kc(ids, ids_n-1-i, v, K-1, branching_vertices_start+branching_vertices_size, vis);

			assert(ids[ids_n-1-i] == v);
		}
	}

	void recursive(ui *ids, ui ids_n, ui u, ui K, ui p, ui branching_vertices_start) {
		ui new_ids_n = ids_n;
		if(u < n) {
			char *t_matrix = matrix + u*n;
			for(ui i = 0;i < new_ids_n;) {
				if(!t_matrix[ids[i]]) {
					-- new_ids_n;
					std::swap(rid[ids[i]], rid[ids[new_ids_n]]);
					std::swap(ids[i], ids[new_ids_n]);
				}
				else ++ i;
			}
			for(ui i = 0;i < new_ids_n;i ++) {
				char *t_matrix = matrix + ids[i]*n;
				for(ui j = new_ids_n;j < ids_n;j ++) if(t_matrix[ids[j]]) -- degree[ids[i]];
			}
		}

		if(K == 1) cnt += p + new_ids_n;
		else if(new_ids_n == 0) cnt += count(p, K);
		else if(new_ids_n == 1) cnt += count(p+1, K);
		else if(p+new_ids_n >= K) {
			ui max_degree = 0, pv = n;
			for(ui i = 0;i < new_ids_n;i ++) if(degree[ids[i]] >= max_degree) {
				max_degree = degree[ids[i]];
				pv = ids[i];
			}

			ui branching_vertices_size = 0;
			add_to_branching_vertices(pv, branching_vertices_start, branching_vertices_size);

			char *t_matrix = matrix + pv*n;
			for(ui i = 0;i < new_ids_n;i ++) if(ids[i] != pv&&!t_matrix[ids[i]]) {
				add_to_branching_vertices(ids[i], branching_vertices_start, branching_vertices_size);
			}

			for(ui i = 0;i < branching_vertices_size;i ++) {
				ui v = branching_vertices_buf[branching_vertices_start+i];
				ui idx = rid[v];
				std::swap(rid[v], rid[ids[new_ids_n-1-i]]);
				std::swap(ids[idx], ids[new_ids_n-1-i]);

				t_matrix = matrix + v*n;
				for(ui j = 0;j < new_ids_n-1-i;j ++) if(t_matrix[ids[j]]) -- degree[ids[j]];

				if(i == 0) recursive(ids, new_ids_n-1-i, v, K, p+1, branching_vertices_start+branching_vertices_size);
				else recursive(ids, new_ids_n-1-i, v, K-1, p, branching_vertices_start+branching_vertices_size);

				assert(ids[new_ids_n-1-i] == v);
			}

			for(ui i = 0;i < branching_vertices_size;i ++) {
				t_matrix = matrix + branching_vertices_buf[branching_vertices_start+i]*n;
				assert(branching_vertices_buf[branching_vertices_start+i] == ids[new_ids_n-1-i]);
				for(ui j = 0;j < new_ids_n-1-i;j ++) if(t_matrix[ids[j]]) ++ degree[ids[j]];
			}
		}

		for(ui i = 0;i < new_ids_n;i ++) {
			char *t_matrix = matrix + ids[i]*n;
			for(ui j = new_ids_n;j < ids_n;j ++) if(t_matrix[ids[j]]) ++ degree[ids[i]];
		}
	}

	void recursive2(ui *ids, ui ids_n, ui u, ui K, ui p, ui branching_vertices_start) {
		ui new_ids_n = ids_n;
		if(u < n) {
			char *t_matrix = matrix + u*n;
			for(ui i = 0;i < new_ids_n;) {
				if(!t_matrix[ids[i]]) {
					-- new_ids_n;
					std::swap(rid[ids[i]], rid[ids[new_ids_n]]);
					std::swap(ids[i], ids[new_ids_n]);
				}
				else ++ i;
			}
			for(ui i = 0;i < new_ids_n;i ++) {
				degree[ids[i]] = 0;
				char *t_matrix = matrix + ids[i]*n;
				for(ui j = 0;j < new_ids_n;j ++) if(t_matrix[ids[j]]) ++ degree[ids[i]];
			}
		}

		if(K == 1) cnt += p + new_ids_n;
		else if(new_ids_n == 0) cnt += count(p, K);
		else if(new_ids_n == 1) cnt += count(p+1, K);
		else if(p+new_ids_n >= K) {
			ui max_degree = 0, pv = n;
			for(ui i = 0;i < new_ids_n;i ++) if(degree[ids[i]] >= max_degree) {
				max_degree = degree[ids[i]];
				pv = ids[i];
			}

			ui branching_vertices_size = 0;
			add_to_branching_vertices(pv, branching_vertices_start, branching_vertices_size);

			char *t_matrix = matrix + pv*n;
			for(ui i = 0;i < new_ids_n;i ++) if(ids[i] != pv&&!t_matrix[ids[i]]) {
				add_to_branching_vertices(ids[i], branching_vertices_start, branching_vertices_size);
			}

			for(ui i = 0;i < branching_vertices_size;i ++) {
				ui v = branching_vertices_buf[branching_vertices_start+i];
				ui idx = rid[v];
				std::swap(rid[v], rid[ids[new_ids_n-1-i]]);
				std::swap(ids[idx], ids[new_ids_n-1-i]);

				if(i == 0) recursive2(ids, new_ids_n-1-i, v, K, p+1, branching_vertices_start+branching_vertices_size);
				else recursive2(ids, new_ids_n-1-i, v, K-1, p, branching_vertices_start+branching_vertices_size);

				assert(ids[new_ids_n-1-i] == v);
			}
		}
	}

	void add_to_branching_vertices(ui u, ui branching_vertices_start, ui &branching_vertices_size) {
		if(branching_vertices_start + branching_vertices_size < branching_vertices_buf.size()) {
			branching_vertices_buf[branching_vertices_start+branching_vertices_size] = u;
		}
		else branching_vertices_buf.push_back(u);
		++ branching_vertices_size;
	}
};


#endif /* PIVOTOR_MATRIX_H_ */
