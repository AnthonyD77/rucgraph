#pragma once

#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

//weightTYPE PrefixalQuery(int s, int t, vector<vector<two_hop_label_v1>>& L, int k) {
//	weightTYPE min_ans = MAX_VALUE;
//	for (int i = 0; i <= k;i++) {
//		auto search_result_s = search_sorted_two_hop_label2(L[s], i);
//		auto search_result_t = search_sorted_two_hop_label2(L[t], i);
//		/*cout << "search_result_s: " << search_result_s.first << " " <<"search_result_t: " << search_result_t.first << endl;*/
//		if (search_result_s.first + search_result_t.first < min_ans) {
//			min_ans = search_result_s.first + search_result_t.first;
//		}
//	}
//	return min_ans;
//}

weightTYPE PrefixalQuery2(vector<vector<two_hop_label_v1>>& L, int source, int terminal, int k) {

	/*return std::numeric_limits<double>::max() is not connected*/

	if (source == terminal) {
		return 0;
	}

	weightTYPE distance = std::numeric_limits<weightTYPE>::max(); // if disconnected, return this large value

	auto vector1_check_pointer = L[source].begin();
	auto vector2_check_pointer = L[terminal].begin();
	auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
	while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end) {
		if (vector1_check_pointer->vertex == vector2_check_pointer->vertex) {
			weightTYPE dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
			if (distance > dis) {
				distance = dis;
			}
			vector1_check_pointer++;
		}
		else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex) {
			vector2_check_pointer++;
			if (vector2_check_pointer->vertex > k) {
				vector2_check_pointer = pointer_L_t_end;
			}
		}
		else {
			vector1_check_pointer++;
			if (vector1_check_pointer->vertex > k) {
				vector1_check_pointer = pointer_L_s_end;
			}
		}
	}

	return distance;

}

void ResumePBFS(graph_hash_of_mixed_weighted& instance_graph, vector<vector<two_hop_label_v1>>& L, int vk, int u, weightTYPE delta) {
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
		//if (vk == v) {
		//	mtx_595[v].lock();
		//}
		//else {
		//	mtx_595[vk].lock(), mtx_595[v].lock();  // 这里有锁的互斥的难题（2个线程加锁是交叉的）
		//}
		weightTYPE PrefixalQueryAns = PrefixalQuery2(L, vk, v, vk);
		//if (vk == v) {
		//	mtx_595[v].unlock();
		//}
		//else {
		//	mtx_595[vk].unlock(), mtx_595[v].unlock();
		//}
		/*cout << "PreQ: " << PrefixalQueryAns << " " << "delta: " << now_delta <<" "<<"vk: "<<vk<< endl;*/

		if (PrefixalQueryAns <= now_delta)
		{
			continue;
		}
			
		/*std::cout << "insert: " << vk <<" "<< now_delta << endl;*/
		//mtx_595[v].lock();
		insert_sorted_two_hop_label(L[v], vk, now_delta);
		//mtx_595[v].unlock();

		auto neis = instance_graph.adj_v_and_ec(v);

		for (auto nei : neis) {
			auto vnei = nei.first;
			affected_label now_nei;
			now_nei.first = vnei;
			now_nei.dis = now_delta + nei.second;
			Q.push(now_nei);
		}
		
	}
}

void WeightDecrease2014(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_new,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	auto& L = mm.L;

	for (int sl = 0; sl < 2; sl++) {
		if (sl == 1) {
			swap(v1, v2);
		}
		for (auto it : L[v1]) {
			int v = it.vertex;
			weightTYPE dis = it.distance + w_new;

			/*the following code cause errors even when there is only 1 thread, why?*/
			//results_dynamic.emplace_back(pool_dynamic.enqueue([&instance_graph, &L, v, v2, dis] {
			//	ResumePBFS(instance_graph, L, v, v2, dis);
			//	return 1; }));

			ResumePBFS(instance_graph, L, v, v2, dis);
		}
	}

	for (auto&& result : results_dynamic) {
		result.get();
	}
	results_dynamic.clear();
}
