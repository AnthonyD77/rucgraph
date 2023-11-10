#pragma once
#include <build_in_progress/HL/dynamic/two_hop_labels_base.h>


/*
in the following codes, it is faster to make each vertex a task in the thread pool, not a fixed vector of vertices as a task.

a fixed vector of vertices as a task:
there are imbalances between workloads of vectors of vertices, and most threads would wait for a few threads to finish in the end
*/


/*
clean L
下面的代码貌似因为内存碎片化导致L占用的内存越来越多
*/

vector<vector<double>> T_clean_L_dynamic;
queue<int> Qid_clean_PPR;

void clean_L_dynamic(vector<vector<two_hop_label_v1>>& L, PPR_type& PPR, int thread_num) {

	int N = L.size();

	vector<vector<double>>().swap(T_clean_L_dynamic);
	T_clean_L_dynamic.resize(thread_num);
	queue<int>().swap(Qid_clean_PPR);
	for (int i = 0; i < thread_num; i++) {
		T_clean_L_dynamic[i].resize(N, std::numeric_limits<weightTYPE>::max());
		Qid_clean_PPR.push(i);
	}

	int div = 1e3;
	if (N < div) {
		div = N;
	}
	for (int jj = 0; jj < N / div; jj++) { // to save RAM of ThreadPool

		cout << "initialize clean_L_dynamic ThreadPool " << jj << endl;
		ThreadPool pool_dynamic(thread_num);
		std::vector<std::future<int>> results_dynamic;

		for (int q = 0; q < div; q++) {
			int v = jj * div + q;

			results_dynamic.emplace_back(
				pool_dynamic.enqueue([v, &L, &PPR] { // pass const type value j to thread; [] can be empty

					mtx_595[max_N_ID_for_mtx_595 - 1].lock();
					int used_id = Qid_clean_PPR.front();
					Qid_clean_PPR.pop();
					mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

					vector<two_hop_label_v1> Lv_final;

					mtx_595[v].lock_shared();
					vector<two_hop_label_v1> Lv = L[v];
					mtx_595[v].unlock_shared();

					auto& T = T_clean_L_dynamic[used_id];

					int size = Lv.size();
					for (int i = 0; i < size; i++) {
						int u = Lv[i].vertex;
						if (Lv[i].distance == std::numeric_limits<weightTYPE>::max()) {
							continue;
						}
						if (v == u) {
							Lv_final.push_back(Lv[i]);
							T[Lv[i].vertex] = Lv[i].distance;
							continue;
						}
						mtx_595[u].lock_shared();
						auto Lu = L[u];
						mtx_595[u].unlock_shared();

						double min_dis = std::numeric_limits<weightTYPE>::max();
						int min_dis_v;
						for (auto label : Lu) {
							double query_dis = label.distance + T[label.vertex];
							if (query_dis < min_dis) {
								min_dis = query_dis;
								min_dis_v = label.vertex;
							}
						}

						if (min_dis > Lv[i].distance + 1e-5) {
							Lv_final.push_back(Lv[i]);
							T[Lv[i].vertex] = Lv[i].distance;
						}
						//else { // with the following clean_PPR, this is not required (otherwise there will be errors)
						//	if (min_dis_v != u) {
						//		mtx_5952[v].lock();
						//		PPR_insert(PPR, v, min_dis_v, u);
						//		mtx_5952[v].unlock();
						//	}
						//	if (min_dis_v != v) {
						//		mtx_5952[u].lock();
						//		PPR_insert(PPR, u, min_dis_v, v);
						//		mtx_5952[u].unlock();
						//	}
						//}
					}

					for (auto label : Lv_final) {
						T[label.vertex] = std::numeric_limits<weightTYPE>::max();
					}

					mtx_595[v].lock();
					L[v] = Lv_final;
					L[v].shrink_to_fit();
					//if (L[v].size() != L[v].capacity()) {
					//	cout << "ssssssssssssss" << endl;
					//	getchar();
					//}
					mtx_595[v].unlock();

					mtx_595[max_N_ID_for_mtx_595 - 1].lock();
					Qid_clean_PPR.push(used_id);
					mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

					return 1; // return to results; the return type must be the same with results
					}));
		}

		for (auto&& result : results_dynamic)
			result.get(); // all threads finish here
		results_dynamic.clear();
	}	
}



/*
clean PPR

it is slower to use thread_vertices_clean_L_dynamic in the following codes
*/

vector<vector<weightTYPE>> d_clean_PPR;
vector<vector<weightTYPE>> T_clean_PPR;
vector<vector<graph_hash_of_mixed_weighted_HL_PLL_v1_handle_t_for_sp>> Q_handles_clean_PPR;

void clean_PPR(graph_v_of_v_idealID& ideal_g, vector<vector<two_hop_label_v1>>& L, std::vector<std::vector<std::pair<int, std::vector<int>>>>& PPR,
	ThreadPool& pool, std::vector<std::future<int>>& results, int thread_num) {

	int N = L.size();
	std::vector<std::vector<std::pair<int, std::vector<int>>>>().swap(PPR);
	PPR.resize(N);

	vector<vector<weightTYPE>>().swap(d_clean_PPR);
	vector<vector<weightTYPE>>().swap(T_clean_PPR);
	d_clean_PPR.resize(thread_num);
	T_clean_PPR.resize(thread_num);
	Q_handles_clean_PPR.resize(thread_num);
	queue<int>().swap(Qid_clean_PPR);
	for (int i = 0; i < thread_num; i++) {
		d_clean_PPR[i].resize(N, std::numeric_limits<weightTYPE>::max());
		T_clean_PPR[i].resize(N, std::numeric_limits<weightTYPE>::max());
		Q_handles_clean_PPR[i].resize(N);
		Qid_clean_PPR.push(i);
	}

	for (int u = 0; u < N; u++) {
		results.emplace_back(
			pool.enqueue([u, &ideal_g, &L, &PPR] { // pass const type value j to thread; [] can be empty

				mtx_595[max_N_ID_for_mtx_595 - 1].lock();
				int used_id = Qid_clean_PPR.front();
				Qid_clean_PPR.pop();
				mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

				auto& d = d_clean_PPR[used_id];
				auto& T = T_clean_PPR[used_id];
				queue<int> d_changed_vertices;
				auto& Q_handles = Q_handles_clean_PPR[used_id];
				PLL_dynamic_node_for_sp node;
				boost::heap::fibonacci_heap<PLL_dynamic_node_for_sp> Q;

				for (auto label : L[u]) {
					T[label.vertex] = label.distance;
				}

				node.vertex = u;
				node.priority_value = 0;
				Q_handles[u] = Q.push(node);
				d[u] = 0;
				d_changed_vertices.push(u);

				while (Q.size() > 0) {
					node = Q.top();
					Q.pop();
					int v = node.vertex;
					if (u <= v) {
						double Lvu = search_sorted_two_hop_label(L[v], u);
						double dv = d[v];
						if (abs(Lvu - dv) < 1e-4) {
							for (auto adj : ideal_g[v]) {
								if (d[adj.first] == std::numeric_limits<weightTYPE>::max()) {
									node.vertex = adj.first;
									node.priority_value = dv + adj.second;
									Q_handles[adj.first] = Q.push(node);
									d[adj.first] = node.priority_value;
									d_changed_vertices.push(adj.first);
								}
								else {
									if (d[adj.first] > dv + adj.second) {
										node.vertex = adj.first;
										node.priority_value = dv + adj.second;
										Q.update(Q_handles[adj.first], node);
										d[adj.first] = node.priority_value;
									}
								}
							}
						}
						else {
							double min_dis = std::numeric_limits<weightTYPE>::max();
							int min_dis_v;
							for (auto label : L[v]) {
								double query_dis = label.distance + T[label.vertex];
								if (query_dis < min_dis) {
									min_dis = query_dis;
									min_dis_v = label.vertex;
								}
							}
							if (min_dis_v != v) {
								mtx_595[u].lock();
								PPR_insert(PPR, u, min_dis_v, v);
								mtx_595[u].unlock();
							}
							if (min_dis_v != u) {
								mtx_595[v].lock();
								PPR_insert(PPR, v, min_dis_v, u);
								mtx_595[v].unlock();
							}
						}
					}
				}

				while (d_changed_vertices.size() > 0) {
					d[d_changed_vertices.front()] = std::numeric_limits<weightTYPE>::max(); // reverse-allocate T values
					d_changed_vertices.pop();
				}

				for (auto label : L[u]) {
					T[label.vertex] = std::numeric_limits<weightTYPE>::max();
				}

				mtx_595[max_N_ID_for_mtx_595 - 1].lock();
				Qid_clean_PPR.push(used_id);
				mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

				return 1; // return to results; the return type must be the same with results
				}));
	}
	for (auto&& result : results)
		result.get(); // all threads finish here
	results.clear();
}

void clean_PPR(graph_v_of_v_idealID& ideal_g, vector<vector<two_hop_label_v1>>& L, std::vector<std::vector<std::pair<int, std::vector<int>>>>& PPR, int thread_num) {

	int N = L.size();
	std::vector<std::vector<std::pair<int, std::vector<int>>>>().swap(PPR);
	PPR.resize(N);

	vector<vector<weightTYPE>>().swap(d_clean_PPR);
	vector<vector<weightTYPE>>().swap(T_clean_PPR);
	d_clean_PPR.resize(thread_num);
	T_clean_PPR.resize(thread_num);
	Q_handles_clean_PPR.resize(thread_num);
	queue<int>().swap(Qid_clean_PPR);
	for (int i = 0; i < thread_num; i++) {
		d_clean_PPR[i].resize(N, std::numeric_limits<weightTYPE>::max());
		T_clean_PPR[i].resize(N, std::numeric_limits<weightTYPE>::max());
		Q_handles_clean_PPR[i].resize(N);
		Qid_clean_PPR.push(i);
	}

	int div = 1e3;
	if (N < div) {
		div = N;
	}
	for (int jj = 0; jj < N / div; jj++) { // to save RAM of ThreadPool

		cout << "initialize clean_PPR ThreadPool " << jj << endl;
		ThreadPool pool_dynamic(thread_num);
		std::vector<std::future<int>> results_dynamic;

		for (int q = 0; q < div; q++) {
			int u = jj * div + q;

			results_dynamic.emplace_back(
				pool_dynamic.enqueue([u, &ideal_g, &L, &PPR] { // pass const type value j to thread; [] can be empty

					mtx_595[max_N_ID_for_mtx_595 - 1].lock();
					int used_id = Qid_clean_PPR.front();
					Qid_clean_PPR.pop();
					mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

					auto& d = d_clean_PPR[used_id];
					auto& T = T_clean_PPR[used_id];
					queue<int> d_changed_vertices;
					auto& Q_handles = Q_handles_clean_PPR[used_id];
					PLL_dynamic_node_for_sp node;
					boost::heap::fibonacci_heap<PLL_dynamic_node_for_sp> Q;

					for (auto label : L[u]) {
						T[label.vertex] = label.distance;
					}

					node.vertex = u;
					node.priority_value = 0;
					Q_handles[u] = Q.push(node);
					d[u] = 0;
					d_changed_vertices.push(u);

					while (Q.size() > 0) {
						node = Q.top();
						Q.pop();
						int v = node.vertex;
						if (u <= v) {
							double Lvu = search_sorted_two_hop_label(L[v], u);
							double dv = d[v];
							if (abs(Lvu - dv) < 1e-4) {
								for (auto adj : ideal_g[v]) {
									if (d[adj.first] == std::numeric_limits<weightTYPE>::max()) {
										node.vertex = adj.first;
										node.priority_value = dv + adj.second;
										Q_handles[adj.first] = Q.push(node);
										d[adj.first] = node.priority_value;
										d_changed_vertices.push(adj.first);
									}
									else {
										if (d[adj.first] > dv + adj.second) {
											node.vertex = adj.first;
											node.priority_value = dv + adj.second;
											Q.update(Q_handles[adj.first], node);
											d[adj.first] = node.priority_value;
										}
									}
								}
							}
							else {
								double min_dis = std::numeric_limits<weightTYPE>::max();
								int min_dis_v;
								for (auto label : L[v]) {
									double query_dis = label.distance + T[label.vertex];
									if (query_dis < min_dis) {
										min_dis = query_dis;
										min_dis_v = label.vertex;
									}
								}
								if (min_dis_v != v) {
									mtx_595[u].lock();
									PPR_insert(PPR, u, min_dis_v, v);
									mtx_595[u].unlock();
								}
								if (min_dis_v != u) {
									mtx_595[v].lock();
									PPR_insert(PPR, v, min_dis_v, u);
									mtx_595[v].unlock();
								}
							}
						}
					}

					while (d_changed_vertices.size() > 0) {
						d[d_changed_vertices.front()] = std::numeric_limits<weightTYPE>::max(); // reverse-allocate T values
						d_changed_vertices.pop();
					}

					for (auto label : L[u]) {
						T[label.vertex] = std::numeric_limits<weightTYPE>::max();
					}

					mtx_595[max_N_ID_for_mtx_595 - 1].lock();
					Qid_clean_PPR.push(used_id);
					mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

					return 1; // return to results; the return type must be the same with results
					}));
		}

		for (auto&& result : results_dynamic)
			result.get(); // all threads finish here
		results_dynamic.clear();
	}	
}