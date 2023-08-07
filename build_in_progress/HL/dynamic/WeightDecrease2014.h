#pragma once

#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

weightTYPE PrefixalQuery(int s, int t, vector<vector<two_hop_label_v1>>* L, int k) {
	weightTYPE min_ans = MAX_VALUE;
	for (int i = 0; i <= k;i++) {
		auto search_result_s = search_sorted_two_hop_label2((*L)[s], i);
		auto search_result_t = search_sorted_two_hop_label2((*L)[t], i);
		/*cout << "search_result_s: " << search_result_s.first << " " <<"search_result_t: " << search_result_t.first << endl;*/
		if (search_result_s.first + search_result_t.first < min_ans) {
			min_ans = search_result_s.first + search_result_t.first;
		}
	}
	return min_ans;
}

void ResumePBFS(graph_hash_of_mixed_weighted* instance_graph, vector<vector<two_hop_label_v1>>* L, int vk, int u, weightTYPE delta) {
	std::queue<affected_label> Q;
	affected_label x;
	x.first = u;
	x.dis = delta;
	Q.push(x);

	while (Q.size()) {
		affected_label now;
		now = Q.front();
		Q.pop();

		auto v = now.first;
		weightTYPE now_delta = now.dis;

		weightTYPE PrefixalQueryAns = PrefixalQuery(vk, v, L, vk);
		/*cout << "PreQ: " << PrefixalQueryAns << " " << "delta: " << now_delta <<" "<<"vk: "<<vk<< endl;*/

		if (PrefixalQueryAns <= now_delta)
		{
			continue;
		}
			
		if (v >= vk)
		{
			/*std::cout << "insert: " << vk <<" "<< now_delta << endl;*/
			insert_sorted_two_hop_label((*L)[v], vk, now_delta);

			auto neis = instance_graph->adj_v_and_ec(v);

			for (auto nei : neis) {
				auto vnei = nei.first;
				affected_label now_nei;
				now_nei.first = vnei;
				now_nei.dis = now_delta + nei.second;
				Q.push(now_nei);
			}
		}
		
	}
}

void WeightDecrease2014(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_new) {

	auto& L = mm.L;

	for (int sl = 0; sl < 2; sl++) {
		if (sl == 1) {
			swap(v1, v2);
		}
		for (auto it : L[v1]) {
			int v = it.vertex;
			weightTYPE dis = it.distance + w_new;

			ResumePBFS(&instance_graph, &mm.L, v, v2, dis);
		}
	}
}
