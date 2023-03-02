#pragma once

#include <build_in_progress/HL/dynamic/graph_hash_of_mixed_weighted_PLL_dynamic.h>

//#define Q(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(mm.L, mm.reduction_measures_2019R2, mm.reduction_measures_2019R1,mm.f_2019R1, instance_graph, x, y)
#define Q(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(mm.L, x, y) // reduction is not used here
#define Q2(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(mm.L, x, y) // reduction is not used here

class changed_label {
public:
	int left, right;
	weightTYPE dis;
};

void ProDecrease(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	std::vector<changed_label>& CL_curr, std::vector<changed_label>& CL_next, weightTYPE weight_change) {

	changed_label xx;

	for (auto it = CL_curr.begin(); it != CL_curr.end(); it++) {
		auto neis = instance_graph.adj_v_and_ec(it->left);
		for (auto nei = neis.begin(); nei != neis.end(); nei++) {
			if (it->right < nei->first) {
				auto query_result = Q2(nei->first, it->right); // query_result is {distance, common hub}
				if (query_result.first > it->dis + nei->second || abs(search_sorted_two_hop_label(mm.L[nei->first], it->right) - it->dis - nei->second - weight_change) < 1e-5) {
					insert_sorted_two_hop_label(mm.L[nei->first], it->right, it->dis + nei->second);
					xx.left = nei->first, xx.right = it->right, xx.dis = it->dis + nei->second;
					CL_next.push_back(xx);
				}
				else { // update PPR
					PPR_insert(mm.PPR, nei->first, query_result.second, it->right);
					PPR_insert(mm.PPR, it->right, query_result.second, nei->first);
				}
			}
		}
	}
}

void WeightDecreaseMaintenance(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new) {

	std::vector<changed_label> CL_curr, CL_next;
	changed_label xx;

	for (auto it = mm.L[v1].begin(); it != mm.L[v1].end(); it++) {
		if (it->vertex <= v2 && (Q(it->vertex, v2) > it->distance + w_new || abs(search_sorted_two_hop_label(mm.L[v2], it->vertex) - it->distance - w_old) < 1e-5)) {
			insert_sorted_two_hop_label(mm.L[v2], it->vertex, it->distance + w_new);
			xx.left = v2, xx.right = it->vertex, xx.dis = it->distance + w_new;
			CL_curr.push_back(xx);
		}
	}

	for (auto it = mm.L[v2].begin(); it != mm.L[v2].end(); it++) {
		if (it->vertex <= v1 && (Q(it->vertex, v1) > it->distance + w_new || abs(search_sorted_two_hop_label(mm.L[v1], it->vertex) - it->distance - w_old) < 1e-5)) { // comparing rank is required here, since the ranks of v1 and v2 are not sorted in this function
			insert_sorted_two_hop_label(mm.L[v1], it->vertex, it->distance + w_new);
			xx.left = v1, xx.right = it->vertex, xx.dis = it->distance + w_new;
			CL_curr.push_back(xx);
		}
	}

	while (CL_curr.size()) {
		ProDecrease(instance_graph, mm, CL_curr, CL_next, w_old - w_new);
		CL_curr = CL_next;
		std::vector<changed_label>().swap(CL_next);
	}
}


