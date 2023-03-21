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

void CLEAR(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	std::vector<changed_label>& CL_curr, std::vector<changed_label>& CL_next, weightTYPE weight_change) {

	changed_label xx;

	for (auto it = CL_curr.begin(); it != CL_curr.end(); it++) {
		auto neis = instance_graph.adj_v_and_ec(it->left);
		for (auto nei = neis.begin(); nei != neis.end(); nei++) {
			if (it->right < nei->first) {
				weightTYPE search_weight = search_sorted_two_hop_label(mm.L[nei->first], it->right);
				if (abs(it->dis + nei->second + weight_change - search_weight) < 1e-5) {
					insert_sorted_two_hop_label(mm.L[nei->first], it->right, std::numeric_limits<weightTYPE>::max());
					xx.left = nei->first, xx.right = it->right, xx.dis = it->dis + nei->second;
					CL_next.push_back(xx);
				}
			}
		}
	}
}

void ProDecrease(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	std::vector<changed_label>& CL_curr, std::vector<changed_label>& CL_next, weightTYPE weight_change) {

	changed_label xx;

	for (auto it = CL_curr.begin(); it != CL_curr.end(); it++) {
		auto neis = instance_graph.adj_v_and_ec(it->left);
		for (auto nei = neis.begin(); nei != neis.end(); nei++) {
			if (it->right < nei->first) {
				auto query_result = Q2(nei->first, it->right); // query_result is {distance, common hub}
				if (query_result.first > it->dis + nei->second) {
					insert_sorted_two_hop_label(mm.L[nei->first], it->right, it->dis + nei->second);
					xx.left = nei->first, xx.right = it->right, xx.dis = it->dis + nei->second;
					CL_next.push_back(xx);
				}
				else {
					auto search_result = search_sorted_two_hop_label2(mm.L[nei->first], it->right);
					if (search_result.first != std::numeric_limits<weightTYPE>::max() && search_result.first > it->dis + nei->second) {
						mm.L[nei->first][search_result.second].distance = it->dis + nei->second;
						xx.left = nei->first, xx.right = it->right, xx.dis = it->dis + nei->second;
						CL_next.push_back(xx);
					}
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
		if (it->vertex <= v2) {
			if (abs(search_sorted_two_hop_label(mm.L[v2], it->vertex) - it->distance - w_old) < 1e-5) {
				insert_sorted_two_hop_label(mm.L[v2], it->vertex, std::numeric_limits<weightTYPE>::max());
				xx.left = v2, xx.right = it->vertex, xx.dis = it->distance + w_new;
				CL_curr.push_back(xx);
			}
		}
	}
	for (auto it = mm.L[v2].begin(); it != mm.L[v2].end(); it++) {
		if (it->vertex <= v1) {
			if (abs(search_sorted_two_hop_label(mm.L[v1], it->vertex) - it->distance - w_old) < 1e-5) {
				insert_sorted_two_hop_label(mm.L[v1], it->vertex, std::numeric_limits<weightTYPE>::max());
				xx.left = v1, xx.right = it->vertex, xx.dis = it->distance + w_new;
				CL_curr.push_back(xx);
			}
		}
	}
	while (CL_curr.size()) {
		CLEAR(instance_graph, mm, CL_curr, CL_next, w_old - w_new);
		CL_curr = CL_next;
		std::vector<changed_label>().swap(CL_next);
	}


	for (auto it = mm.L[v1].begin(); it != mm.L[v1].end(); it++) {	
		if (it->vertex <= v2) {
			auto query_result = Q2(it->vertex, v2); // query_result is {distance, common hub}
			if (query_result.first > it->distance + w_new) {
				insert_sorted_two_hop_label(mm.L[v2], it->vertex, it->distance + w_new);
				xx.left = v2, xx.right = it->vertex, xx.dis = it->distance + w_new;
				CL_curr.push_back(xx);
			}
			else {
				PPR_insert(mm.PPR, v2, query_result.second, it->vertex);
				PPR_insert(mm.PPR, it->vertex, query_result.second, v2);
			}
		}
	}
	for (auto it = mm.L[v2].begin(); it != mm.L[v2].end(); it++) {
		if (it->vertex <= v1) {
			auto query_result = Q2(it->vertex, v1); // query_result is {distance, common hub}
			if (query_result.first > it->distance + w_new) {
				insert_sorted_two_hop_label(mm.L[v1], it->vertex, it->distance + w_new);
				xx.left = v1, xx.right = it->vertex, xx.dis = it->distance + w_new;
				CL_curr.push_back(xx);
			}
			else {
				PPR_insert(mm.PPR, v1, query_result.second, it->vertex);
				PPR_insert(mm.PPR, it->vertex, query_result.second, v1);
			}
		}
	}
	while (CL_curr.size()) {
		ProDecrease(instance_graph, mm, CL_curr, CL_next, w_old - w_new);
		CL_curr = CL_next;
		std::vector<changed_label>().swap(CL_next);
	}
}


