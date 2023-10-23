#pragma once

#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

void WeightDecreaseMaintenance_improv_step1(int v1, int v2, weightTYPE w_new, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR, std::vector<affected_label>* CL,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	for (int sl = 0; sl < 2; sl++) {
		if (sl == 1) {
			swap(v1, v2);
		}
		for (auto it : (*L)[v1]) {
			if (it.vertex <= v2) {
				results_dynamic.emplace_back(pool_dynamic.enqueue([it, v2, L, PPR, w_new, CL] {

					auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, it.vertex, v2); // query_result is {distance, common hub}
					if (query_result.first > it.distance + w_new) {
						mtx_595_1.lock();
						CL->push_back(affected_label{ v2 , it.vertex, it.distance + w_new });
						mtx_595_1.unlock();
					}
					else {
						auto search_result = search_sorted_two_hop_label((*L)[v2], it.vertex);
						if (search_result > it.distance + w_new && search_result != MAX_VALUE) {
							mtx_595_1.lock();
							CL->push_back(affected_label{ v2, it.vertex, it.distance + w_new });
							mtx_595_1.unlock();
						}
						if (query_result.second != it.vertex) {
							mtx_5952[v2].lock();
							PPR_insert(*PPR, v2, query_result.second, it.vertex);
							mtx_5952[v2].unlock();
						}
						if (query_result.second != v2) {
							mtx_5952[it.vertex].lock();
							PPR_insert(*PPR, it.vertex, query_result.second, v2);
							mtx_5952[it.vertex].unlock();
						}
					}

					return 1; }));
			}
		}
	}

	for (auto&& result : results_dynamic) {
		result.get();
	}
	results_dynamic.clear();
}

void DIFFUSE(graph_v_of_v_idealID& instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR, std::vector<affected_label>& CL,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	//cout << "CL.size(): " << CL.size() << endl;

	for (auto it : CL) {
		results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, &instance_graph, PPR] {

			mtx_595_1.lock();
			int current_tid = Qid_595.front();
			Qid_595.pop();
			mtx_595_1.unlock();

			int u = it.first, v = it.second;
			weightTYPE du = it.dis;

			mtx_595[v].lock();
			auto Lv = (*L)[v]; // to avoid interlocking
			mtx_595[v].unlock();

			vector<int> Dis_changed;
			auto& DIS = Dis[current_tid];
			auto& Q_HANDLES = Q_handles[current_tid];
			auto& Q_VALUE = Q_value[current_tid];
			DIS[u] = { du, v }; // <distance, hub responsible for this distance>
			Dis_changed.push_back(u);

			boost::heap::fibonacci_heap<node_for_DIFFUSE> Q;
			Q_HANDLES[u] = Q.push(node_for_DIFFUSE(u, du));
			Q_VALUE[u] = du;

			while (!Q.empty()) {

				node_for_DIFFUSE temp2 = Q.top();
				int x = temp2.index;
				weightTYPE dx = temp2.disx;
				Q.pop();
				Q_VALUE[x] = MAX_VALUE;

				mtx_595[x].lock();
				insert_sorted_two_hop_label((*L)[x], v, dx);
				mtx_595[x].unlock();

				for (auto& nei : instance_graph[x]) {
					int xnei = nei.first;
					weightTYPE d_new = dx + nei.second;

					if (v < xnei) {
						if (DIS[xnei].first == -1) {
							mtx_595[xnei].lock();
							DIS[xnei] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4((*L)[xnei], Lv);
							mtx_595[xnei].unlock();
							Dis_changed.push_back(xnei);
						}
						if (DIS[xnei].first > d_new + 1e-5) {
							DIS[xnei] = { d_new, v };
							if (Q_VALUE[xnei] == MAX_VALUE) {
								Q_HANDLES[xnei] = Q.push(node_for_DIFFUSE(xnei, d_new));
							}
							else {
								Q.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
							}
							Q_VALUE[xnei] = d_new;
						}
						else {
							mtx_595[xnei].lock();
							auto search_result = search_sorted_two_hop_label2((*L)[xnei], v);
							mtx_595[xnei].unlock();
							if (search_result.second != -1 && std::min(search_result.first, Q_VALUE[xnei]) > d_new + 1e-5) {
								//cout << "here " << Q_VALUE[xnei] << endl;
								if (Q_VALUE[xnei] == MAX_VALUE) {
									Q_HANDLES[xnei] = Q.push(node_for_DIFFUSE(xnei, d_new));
								}
								else {
									Q.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
								}
								Q_VALUE[xnei] = d_new;
							}
							if (DIS[xnei].second != v) {
								mtx_5952[xnei].lock();
								PPR_insert(*PPR, xnei, DIS[xnei].second, v);
								mtx_5952[xnei].unlock();
							}
							if (DIS[xnei].second != xnei) {
								mtx_5952[v].lock();
								PPR_insert(*PPR, v, DIS[xnei].second, xnei);
								mtx_5952[v].unlock();
							}
						}
					}
				}
			}

			for (int i : Dis_changed) {
				DIS[i] = { -1, -1 };
			}

			mtx_595_1.lock();
			Qid_595.push(current_tid);
			mtx_595_1.unlock();

			return 1; }));
	}

	for (auto&& result : results_dynamic) {
		result.get();
	}
	results_dynamic.clear();
}

void WeightDecreaseMaintenance_improv(graph_v_of_v_idealID& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	std::vector<affected_label> CL;
	WeightDecreaseMaintenance_improv_step1(v1, v2, w_new, &mm.L, &mm.PPR, &CL, pool_dynamic, results_dynamic);

	DIFFUSE(instance_graph, &mm.L, &mm.PPR, CL, pool_dynamic, results_dynamic);
}
