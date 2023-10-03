#pragma once
#include <build_in_progress/HL/dynamic/two_hop_labels_base.h>




/*clean L*/

pair<weightTYPE, int> clean_L_extract_distance(vector<two_hop_label_v1>& L_s, vector<two_hop_label_v1>& L_t) {

	/*return std::numeric_limits<double>::max() is not connected*/

	weightTYPE distance = std::numeric_limits<weightTYPE>::max(); // if disconnected, return this large value
	int common_hub;

	auto vector1_check_pointer = L_s.begin(), vector2_check_pointer = L_t.begin();
	auto pointer_L_s_end = L_s.end(), pointer_L_t_end = L_t.end();
	while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end) {
		if (vector1_check_pointer->vertex == vector2_check_pointer->vertex) {
			weightTYPE dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
			if (distance > dis) {
				distance = dis;
				common_hub = vector1_check_pointer->vertex;
			}
			vector1_check_pointer++;
		}
		else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex) {
			vector2_check_pointer++;
		}
		else {
			vector1_check_pointer++;
		}
	}

	return { distance , common_hub };
}

void clean_L_element(int v, vector<vector<two_hop_label_v1>>& L, PPR_type& PPR) {

	vector<two_hop_label_v1> Lv_final;
	mtx_595[v].lock();
	vector<two_hop_label_v1> Lv = L[v];
	mtx_595[v].unlock();

	int size = Lv.size();
	for (int i = 0; i < size; i++) {
		int u = Lv[i].vertex;
		if (v == u) {
			Lv_final.push_back(Lv[i]);
			continue;
		}
		mtx_595[u].lock();
		pair<weightTYPE, int> query = clean_L_extract_distance(Lv_final, L[u]);
		mtx_595[u].unlock();
		if (query.first > Lv[i].distance + 1e-5) {
			Lv_final.push_back(Lv[i]);
		}
		else {
			if (query.second != u) {
				mtx_5952[v].lock();
				PPR_insert(PPR, v, query.second, u);
				mtx_5952[v].unlock();
			}
			if (query.second != v) {
				mtx_5952[u].lock();
				PPR_insert(PPR, u, query.second, v);
				mtx_5952[u].unlock();
			}
		}
	}
	mtx_595[v].lock();
	vector<two_hop_label_v1>(Lv_final).swap(L[v]);
	mtx_595[v].unlock();
}

void clean_L_dynamic(vector<vector<two_hop_label_v1>>& L, PPR_type& PPR, ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	int N = L.size();

	for (int target_v = 0; target_v < N; target_v++) {
		results_dynamic.emplace_back(
			pool_dynamic.enqueue([target_v, &L, &PPR] { // pass const type value j to thread; [] can be empty
				clean_L_element(target_v, L, PPR);
				return 1; // return to results; the return type must be the same with results
				}));
	}
	for (auto&& result : results_dynamic)
		result.get(); // all threads finish here
	results_dynamic.clear();
}







vector<vector<weightTYPE>> d_clean_PPR;
vector<vector<weightTYPE>> T_clean_PPR;
vector<vector<graph_hash_of_mixed_weighted_HL_PLL_v1_handle_t_for_sp>> Q_handles_clean_PPR;
queue<int> Qid_clean_PPR;

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
			pool.enqueue([u, &ideal_g , &L, &PPR] { // pass const type value j to thread; [] can be empty

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
								if (d[adj.first] == std::numeric_limits<weightTYPE>::max()) { //尚未到达的点
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
							for (auto label: L[v]) {
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