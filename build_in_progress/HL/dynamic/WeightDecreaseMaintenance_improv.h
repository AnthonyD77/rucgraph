#pragma once

#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

void WeightDecreaseMaintenance_improv_step1(int v1, int v2, weightTYPE w_new, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR, std::vector<affected_label>* CL) {

	for (int sl = 0; sl < 2; sl++) {
		if (sl == 1) {
			swap(v1, v2);
		}
		for (auto it : (*L)[v1]) {
			if (it.vertex <= v2) {
				auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, it.vertex, v2); // query_result is {distance, common hub}
				if (query_result.first > it.distance + w_new) {
					CL->push_back(affected_label{ v2 , it.vertex, it.distance + w_new });
				}
				else {
					auto search_result = search_sorted_two_hop_label((*L)[v2], it.vertex);
					if (search_result > it.distance + w_new && search_result != MAX_VALUE) {
						CL->push_back(affected_label{ v2, it.vertex, it.distance + w_new });
					}
					if (query_result.second != it.vertex) {
						PPR_insert(*PPR, v2, query_result.second, it.vertex);
					}
					if (query_result.second != v2) {
						PPR_insert(*PPR, it.vertex, query_result.second, v2);
					}
				}
			}
		}
	}
}

void DIFFUSE(graph_hash_of_mixed_weighted* instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR, std::vector<affected_label>& CL) {

	//cout << "CL.size(): " << CL.size() << endl;

	for (auto it : CL) {
		int u = it.first, v = it.second;
		weightTYPE du = it.dis;

		vector<int> Dis_changed;
		int N=instance_graph->hash_of_vectors.size();
		vector<pair<weightTYPE,int> > DIS(N,{-1,-1});
		vector<handle_t_for_DIFFUSE> Q_HANDLES(N);
		vector<weightTYPE> Q_VALUE(N,MAX_VALUE);
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

			insert_sorted_two_hop_label((*L)[x], v, dx);

			auto neis = instance_graph->adj_v_and_ec(x);
			for (auto nei : neis) {
				int xnei = nei.first;
				weightTYPE d_new = dx + nei.second;

				if (v < xnei) {
					if (DIS[xnei].first == -1) {
						DIS[xnei] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, xnei, v);
						Dis_changed.push_back(xnei);
					}
					if (DIS[xnei].first > d_new - 1e-5) {
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
						auto search_result = search_sorted_two_hop_label2((*L)[xnei], v);
						if (search_result.second != -1 && std::min(search_result.first, Q_VALUE[xnei]) > d_new) {
							if (Q_VALUE[xnei] == MAX_VALUE) {
								Q_HANDLES[xnei] = Q.push(node_for_DIFFUSE(xnei, d_new));
							}
							else {
								Q.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
							}
							Q_VALUE[xnei] = d_new;
						}
						if (DIS[xnei].second != v) {
							PPR_insert(*PPR, xnei, DIS[xnei].second, v);
						}
						if (DIS[xnei].second != xnei) {
							PPR_insert(*PPR, v, DIS[xnei].second, xnei);
						}
					}
				}
			}
		}
	}
}

void WeightDecreaseMaintenance_improv(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_new) {

	std::vector<affected_label> CL;
	WeightDecreaseMaintenance_improv_step1(v1, v2, w_new, &mm.L, &mm.PPR, &CL);

	DIFFUSE(&instance_graph, &mm.L, &mm.PPR, CL);
}
