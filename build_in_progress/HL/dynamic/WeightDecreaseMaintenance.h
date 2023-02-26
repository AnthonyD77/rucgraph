#pragma once

#include <build_in_progress/HL/dynamic/graph_hash_of_mixed_weighted_PLL_dynamic.h>
#include <build_in_progress/HL/sort_v/graph_hash_of_mixed_weighted_update_vertexIDs_by_degrees.h>
#include <graph_hash_of_mixed_weighted/two_graphs_operations/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID_2.h>
#include <graph_hash_of_mixed_weighted/random_graph/graph_hash_of_mixed_weighted_generate_random_graph.h>
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_read_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_save_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/common_algorithms/graph_hash_of_mixed_weighted_shortest_paths.h>


//#define Q(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(mm.L, mm.reduction_measures_2019R2, mm.reduction_measures_2019R1,mm.f_2019R1, instance_graph, x, y)
#define Q(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(mm.L, x, y) // reduction is not used here
#define Q2(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(mm.L, x, y) // reduction is not used here

class changed_label {
public:
	int left, right;
	weightTYPE dis;
};

void ProDecrease(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	std::vector<changed_label>& CL_curr, std::vector<changed_label>& CL_next) {

	changed_label xx;

	for (auto it = CL_curr.begin(); it != CL_curr.end(); it++) {
		auto neis = instance_graph.adj_v_and_ec(it->left);
		for (auto nei = neis.begin(); nei != neis.end(); nei++) {
			if (it->right < nei->first) {
				auto query_result = Q2(nei->first, it->right);
				if (query_result.first > it->dis + nei->second) {
					sorted_two_hop_label_v1_vector_binary_insert_or_update(mm.L[nei->first], it->right, it->dis + nei->second);
					xx.left = nei->first, xx.right = it->right, xx.dis = it->dis + nei->second;
					CL_next.push_back(xx);
				}
				else {
					PPR_insert(mm.PPR, nei->first, query_result.second, it->right);
					PPR_insert(mm.PPR, it->right, query_result.second, nei->first);
				}
			}
		}
	}
}

void WeightDecreaseMaintenance(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new) {

	if (Q(v1, v2) < w_new)
		return;

	std::vector<changed_label> CL_curr, CL_next;
	changed_label xx;

	for (auto it = mm.L[v1].begin(); it != mm.L[v1].end(); it++) {
		if (it->vertex < v2 && Q(it->vertex, v2) > it->distance + w_new) {
			sorted_two_hop_label_v1_vector_binary_insert_or_update(mm.L[v2], it->vertex, it->distance + w_new);
			xx.left = v2, xx.right = it->vertex, xx.dis = it->distance + w_new;
			CL_curr.push_back(xx);
		}
	}

	for (auto it = mm.L[v2].begin(); it != mm.L[v2].end(); it++) {

		//if (it->vertex == 1) {
		//	cout << it->vertex << " " << it->distance << " " << v1 << endl;
		//}

		if (it->vertex < v1 && Q(it->vertex, v1) > it->distance + w_new) {
			sorted_two_hop_label_v1_vector_binary_insert_or_update(mm.L[v1], it->vertex, it->distance + w_new);
			xx.left = v1, xx.right = it->vertex, xx.dis = it->distance + w_new;
			CL_curr.push_back(xx);
		}
	}

	while (CL_curr.size()) {
		ProDecrease(instance_graph, mm, CL_curr, CL_next);
		CL_curr = CL_next;
		std::vector<changed_label>().swap(CL_next);
	}
}


