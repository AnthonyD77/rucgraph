#pragma once

#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

void PI11(graph_v_of_v_idealID& instance_graph, vector<vector<two_hop_label_v1>>* L,
	std::vector<affected_label>& al1_curr, std::vector<affected_label>* al1_next, ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	for (auto it : al1_curr) {
		results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, al1_next, &instance_graph] {
			for (auto nei : instance_graph[it.first]) {
				weightTYPE search_weight = search_sorted_two_hop_label((*L)[nei.first], it.second);
				if (abs(it.dis + nei.second - search_weight) < 1e-5) {
					mtx_595_1.lock();
					al1_next->push_back(affected_label(nei.first, it.second, search_weight));
					mtx_595_1.unlock();
				}
			}
			insert_sorted_two_hop_label((*L)[it.first], it.second, MAX_VALUE); // this does not change the size of L[it->first] here, so does not need to lock here

			return 1; }));
	}

	for (auto&& result : results_dynamic) {
		result.get();
	}
	results_dynamic.clear();
}

void PI12(graph_v_of_v_idealID& instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR,
	std::vector<affected_label>& al1_curr, std::vector<pair_label>* al2_next, ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	for (auto it : al1_curr) {
		results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, PPR, al2_next, &instance_graph] {

			int v = it.first, u = it.second;
			mtx_5952[v].lock();
			std::vector<int> temp = PPR_retrieve(*PPR, v, u);
			mtx_5952[v].unlock();
			PPR_binary_operations_insert(temp, u);

			mtx_595[v].lock();
			auto Lv = (*L)[v]; // to avoid interlocking
			mtx_595[v].unlock();

			for (auto t : temp) {
				if (v < t) {
					weightTYPE d1 = MAX_VALUE;
					for (auto nei : instance_graph[t]) {
						mtx_595[nei.first].lock();
						d1 = min(d1, search_sorted_two_hop_label((*L)[nei.first], v) + (weightTYPE)nei.second);
						mtx_595[nei.first].unlock();
					}
					mtx_595[t].lock();
					auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4((*L)[t], Lv);
					mtx_595[t].unlock();
					if (query_result.first > d1) {
						mtx_595[t].lock();
						insert_sorted_two_hop_label((*L)[t], v, d1);
						mtx_595[t].unlock();
						mtx_595_1.lock();
						al2_next->push_back(pair_label(t, v));
						mtx_595_1.unlock();
					}
					else {
						if (query_result.second != v) {
							mtx_5952[t].lock();
							PPR_insert(*PPR, t, query_result.second, v);
							mtx_5952[t].unlock();
						}
						if (query_result.second != t) {
							mtx_5952[v].lock();
							PPR_insert(*PPR, v, query_result.second, t);
							mtx_5952[v].unlock();
						}
					}
				}
				if (t < v) {
					weightTYPE d1 = MAX_VALUE;
					for (auto nei : instance_graph[v]) {
						mtx_595[nei.first].lock();
						d1 = min(d1, search_sorted_two_hop_label((*L)[nei.first], t) + (weightTYPE)nei.second);
						mtx_595[nei.first].unlock();
					}
					mtx_595[t].lock();
					auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4((*L)[t], Lv);
					mtx_595[t].unlock();
					if (query_result.first > d1) {
						mtx_595[v].lock();
						insert_sorted_two_hop_label((*L)[v], t, d1);
						mtx_595[v].unlock();
						mtx_595_1.lock();
						al2_next->push_back(pair_label(v, t));
						mtx_595_1.unlock();
					}
					else {
						if (query_result.second != v) {
							mtx_5952[t].lock();
							PPR_insert(*PPR, t, query_result.second, v);
							mtx_5952[t].unlock();
						}
						if (query_result.second != t) {
							mtx_5952[v].lock();
							PPR_insert(*PPR, v, query_result.second, t);
							mtx_5952[v].unlock();
						}
					}
				}
			}

			return 1; }));
	}

	for (auto&& result : results_dynamic) {
		result.get();
	}
	results_dynamic.clear();
}

void PI22(graph_v_of_v_idealID& instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR,
	std::vector<pair_label>& al2_curr, std::vector<pair_label>* al2_next, ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	for (auto it = al2_curr.begin(); it != al2_curr.end(); it++) {
		results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, PPR, al2_next, &instance_graph] {

			mtx_595[it->second].lock();
			auto Lxx = (*L)[it->second]; // to avoid interlocking
			mtx_595[it->second].unlock();

			for (auto nei : instance_graph[it->first]) {
				if (nei.first > it->second) {
					mtx_595[it->first].lock();
					weightTYPE search_result = search_sorted_two_hop_label((*L)[it->first], it->second) + nei.second;
					mtx_595[it->first].unlock();
					mtx_595[nei.first].lock();
					auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4((*L)[nei.first], Lxx); 
					mtx_595[nei.first].unlock();
					if (query_result.first + 1e-3 >= search_result) {
						mtx_595[nei.first].lock();
						insert_sorted_two_hop_label((*L)[nei.first], it->second, search_result);
						mtx_595[nei.first].unlock();
						mtx_595_1.lock();
						al2_next->push_back(pair_label(nei.first, it->second));
						mtx_595_1.unlock();
					}
					else {
						if (query_result.second != it->second) {
							mtx_5952[nei.first].lock();
							PPR_insert(*PPR, nei.first, query_result.second, it->second);
							mtx_5952[nei.first].unlock();
						}
						if (query_result.second != nei.first) {
							mtx_5952[it->second].lock();
							PPR_insert(*PPR, it->second, query_result.second, nei.first);
							mtx_5952[it->second].unlock();
						}
					}
				}
			}

			return 1; }));
	}

	for (auto&& result : results_dynamic) {
		result.get();
	}
	results_dynamic.clear();
}

void WeightIncrease2021(graph_v_of_v_idealID& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	std::vector<affected_label> al1_curr, al1_next;
	std::vector<pair_label> al2_curr, al2_next;

	for (auto it : mm.L[v1]) {
		if (it.vertex <= v2 && abs(search_sorted_two_hop_label(mm.L[v2], it.vertex) - it.distance - w_old) < 1e-5) {
			al1_curr.push_back(affected_label(v2, it.vertex, it.distance + w_old));
		}
	}
	for (auto it : mm.L[v2]) {
		if (it.vertex <= v1 && abs(search_sorted_two_hop_label(mm.L[v1], it.vertex) - it.distance - w_old) < 1e-5) {
			al1_curr.push_back(affected_label(v1, it.vertex, it.distance + w_old));
		}
	}

	while (al1_curr.size() || al2_curr.size()) {
		PI11(instance_graph, &mm.L, al1_curr, &al1_next, pool_dynamic, results_dynamic);
		PI12(instance_graph, &mm.L, &mm.PPR, al1_curr, &al2_next, pool_dynamic, results_dynamic);
		PI22(instance_graph, &mm.L, &mm.PPR, al2_curr, &al2_next, pool_dynamic, results_dynamic);
		al1_curr = al1_next;
		al2_curr = al2_next;
		std::vector<affected_label>().swap(al1_next);
		std::vector<pair_label>().swap(al2_next);
	}
}

