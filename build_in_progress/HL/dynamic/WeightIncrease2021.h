#pragma once

#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

void PI11(graph_hash_of_mixed_weighted* instance_graph, vector<vector<two_hop_label_v1>>* L,
	std::vector<affected_label>& al1_curr, std::vector<affected_label>* al1_next) {

	for (auto it : al1_curr) {
		auto v_neis = instance_graph->adj_v_and_ec(it.first);
		for (auto nei = v_neis.begin(); nei != v_neis.end(); nei++) {
			weightTYPE search_weight = search_sorted_two_hop_label((*L)[nei->first], it.second);
			if (abs(it.dis + nei->second - search_weight) < 1e-5) {
				al1_next->push_back(affected_label(nei->first, it.second, search_weight));
			}
		}
		insert_sorted_two_hop_label((*L)[it.first], it.second, MAX_VALUE); // this does not change the size of L[it->first] here, so does not need to lock here
	}
}

void PI12(graph_hash_of_mixed_weighted* instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR,
	std::vector<affected_label>& al1_curr, std::vector<pair_label>* al2_next) {

	for (auto it : al1_curr) {
		int v = it.first, u = it.second;
		std::vector<int> temp = PPR_retrieve(*PPR, v, u);
		PPR_binary_operations_insert(temp, u);

		for (auto t : temp) {
			if (v < t) {
				weightTYPE d1 = MAX_VALUE;
				auto neis = instance_graph->adj_v_and_ec(t);
				for (auto nei : neis) {
					d1 = min(d1, search_sorted_two_hop_label((*L)[nei.first], v) + (weightTYPE)nei.second);
				}
				auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, t, v);
				if (query_result.first > d1) {
					insert_sorted_two_hop_label((*L)[t], v, d1);
					al2_next->push_back(pair_label(t, v));
				}
				else {
					if (query_result.second != v) {
						PPR_insert(*PPR, t, query_result.second, v);
					}
					if (query_result.second != t) {
						PPR_insert(*PPR, v, query_result.second, t);
					}
				}
			}
			if (t < v) {
				weightTYPE d1 = MAX_VALUE;
				auto neis = instance_graph->adj_v_and_ec(v);
				for (auto nei : neis) {
					d1 = min(d1, search_sorted_two_hop_label((*L)[nei.first], t) + (weightTYPE)nei.second);
				}
				auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, v, t);
				if (query_result.first > d1) {
					insert_sorted_two_hop_label((*L)[v], t, d1);
					al2_next->push_back(pair_label(v, t));
				}
				else {
					if (query_result.second != v) {
						PPR_insert(*PPR, t, query_result.second, v);
					}
					if (query_result.second != t) {
						PPR_insert(*PPR, v, query_result.second, t);
					}
				}
			}
		}

	}
}

void PI22(graph_hash_of_mixed_weighted* instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR,
	std::vector<pair_label>& al2_curr, std::vector<pair_label>* al2_next) {

	for (auto it = al2_curr.begin(); it != al2_curr.end(); it++) {
		auto v_neis = instance_graph->adj_v_and_ec(it->first);
		for (auto nei = v_neis.begin(); nei != v_neis.end(); nei++) {
			if (nei->first > it->second) {
				weightTYPE search_result = search_sorted_two_hop_label((*L)[it->first], it->second) + nei->second;
				auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, nei->first, it->second);
				if (query_result.first + 1e-3 >= search_result) {
					insert_sorted_two_hop_label((*L)[nei->first], it->second, search_result);
					al2_next->push_back(pair_label(nei->first, it->second));
				}
				else {
					if (query_result.second != it->second) {
						PPR_insert(*PPR, nei->first, query_result.second, it->second);
					}
					if (query_result.second != nei->first) {
						PPR_insert(*PPR, it->second, query_result.second, nei->first);
					}
				}
			}
		}

	}
}

void WeightIncrease2021(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_old) {

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
		PI11(&instance_graph, &mm.L, al1_curr, &al1_next);
		PI12(&instance_graph, &mm.L, &mm.PPR, al1_curr, &al2_next);
		PI22(&instance_graph, &mm.L, &mm.PPR, al2_curr, &al2_next);
		al1_curr = al1_next;
		al2_curr = al2_next;
		std::vector<affected_label>().swap(al1_next);
		std::vector<pair_label>().swap(al2_next);
	}
}

