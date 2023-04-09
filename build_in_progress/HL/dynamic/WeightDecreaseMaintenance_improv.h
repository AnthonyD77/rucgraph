#pragma once

#include <build_in_progress/HL/dynamic/graph_hash_of_mixed_weighted_PLL_dynamic.h>

void WeightDecreaseMaintenance_improv_step1(int v1, int v2, weightTYPE w_new, vector<vector<two_hop_label_v1>>& L, PPR_type& PPR, std::vector<affected_label>& CL, int thread_num) {

	for (int sl = 0; sl < 2; sl++) {
		if (sl == 1) {
			swap(v1, v2);
		}
		for (auto it : L[v1]) {
			if (it.vertex <= v2) {
				auto query_result = Query2(it.vertex, v2); // query_result is {distance, common hub}
				if (query_result.first > it.distance + w_new) {
					CL.push_back(affected_label{ v2 , it.vertex, it.distance + w_new });
				}
				else {
					auto search_result = search_sorted_two_hop_label2(L[v2], it.vertex);
					if (search_result.first > it.distance + w_new && search_result.first != std::numeric_limits<weightTYPE>::max()) {
						CL.push_back(affected_label{ v2, it.vertex, it.distance + w_new });
					}
					PPR_insert(PPR, v2, query_result.second, it.vertex);
					PPR_insert(PPR, it.vertex, query_result.second, v2);
				}
			}
		}
	}

}

void WeightDecreaseMaintenance_improv_step1_parallelVersion(int v1, int v2, weightTYPE w_new, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR, std::vector<affected_label>* CL, int thread_num) {

	ThreadPool pool(thread_num);
	std::vector<std::future<int>> results;

	for (int sl = 0; sl < 2; sl++) {
		if (sl == 1) {
			int x = v1;
			v1 = v2;
			v2 = x;
		}
		for (auto it : (*L)[v1]) {
			if (it.vertex <= v2) {
				results.emplace_back(pool.enqueue([it, v2, L, PPR, w_new, CL] {

					auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, it.vertex, v2); // query_result is {distance, common hub}
					if (query_result.first > it.distance + w_new) {
						mtx_595_1.lock();
						CL->push_back(affected_label{ v2 , it.vertex, it.distance + w_new });
						mtx_595_1.unlock();
					}
					else {
						auto search_result = search_sorted_two_hop_label2((*L)[v2], it.vertex);
						if (search_result.first > it.distance + w_new && search_result.first != std::numeric_limits<weightTYPE>::max()) {
							mtx_595_1.lock();
							CL->push_back(affected_label{ v2, it.vertex, it.distance + w_new });
							mtx_595_1.unlock();
						}
						mtx_5952[v2].lock();
						PPR_insert(*PPR, v2, query_result.second, it.vertex);
						mtx_5952[v2].unlock();
						mtx_5952[it.vertex].lock();
						PPR_insert(*PPR, it.vertex, query_result.second, v2);
						mtx_5952[it.vertex].unlock();
					}


					return 1; }));
			}
		}
	}

	for (auto&& result : results) {
		result.get();
	}
	results.clear();
}


void DIFFUSE(graph_hash_of_mixed_weighted* instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR, std::vector<affected_label>& CL, int thread_num) {

	ThreadPool pool(thread_num);
	std::vector<std::future<int>> results;

	//cout << "CL.size(): " << CL.size() << endl;

	for (auto it : CL) {
		results.emplace_back(pool.enqueue([it, L, instance_graph, PPR] {

			mtx_595_1.lock();
			int current_tid = Qid_595.front();
			Qid_595.pop();
			mtx_595_1.unlock();

			int u = it.first, v = it.second;
			weightTYPE du = it.dis;
			vector<int> Dis_changed;
			Dis[current_tid][u] = pair(du, v);
			Dis_changed.push_back(u);

			boost::heap::fibonacci_heap<node_for_DIFFUSE> Q;
			Q_handles[current_tid][u] = Q.push(node_for_DIFFUSE(u, du));
			Q_value[current_tid][u] = du;

			while (!Q.empty()) {
				node_for_DIFFUSE temp2 = Q.top();
				int x = temp2.index;
				weightTYPE dx = temp2.disx;
				Q.pop();
				Q_value[current_tid][x] = MAX_VALUE;

				mtx_595[x].lock();
				insert_sorted_two_hop_label((*L)[x], v, dx);
				mtx_595[x].unlock();

				auto neis = instance_graph->adj_v_and_ec(x);
				for (auto nei : neis) {
					int xnei = nei.first;
					double dxxnei = nei.second;

					if (v < xnei) {
						if (Dis[current_tid][xnei].first == -1) {
							mtx_595[xnei].lock(), mtx_595[v].lock();
							Dis[current_tid][xnei] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, xnei, v);
							mtx_595[xnei].unlock(), mtx_595[v].unlock();
							Dis_changed.push_back(xnei);
						}
						if (Dis[current_tid][xnei].first > dx + dxxnei - 1e-5) {
							Dis[current_tid][xnei] = pair(dx + dxxnei, v);
							if (Q_value[current_tid][xnei] == MAX_VALUE) {
								Q_handles[current_tid][xnei] = Q.push(node_for_DIFFUSE(xnei, dx + dxxnei));
							}
							else {
								Q.update(Q_handles[current_tid][xnei], node_for_DIFFUSE(xnei, dx + dxxnei));
							}
							Q_value[current_tid][xnei] = dx + dxxnei;
						}
						else {
							mtx_595[xnei].lock();
							auto search_result = search_sorted_two_hop_label2((*L)[xnei], v);
							mtx_595[xnei].unlock();
							if (search_result.second != -1 && std::min(search_result.first, Q_value[current_tid][xnei]) > dx + dxxnei) {
								if (Q_value[current_tid][xnei] == MAX_VALUE) {
									Q_handles[current_tid][xnei] = Q.push(node_for_DIFFUSE(xnei, dx + dxxnei));
								}
								else {
									Q.update(Q_handles[current_tid][xnei], node_for_DIFFUSE(xnei, dx + dxxnei));
								}
								Q_value[current_tid][xnei] = dx + dxxnei;
							}
							mtx_5952[xnei].lock();
							PPR_insert(*PPR, xnei, Dis[current_tid][xnei].second, v);
							mtx_5952[xnei].unlock();
							mtx_5952[v].lock();
							PPR_insert(*PPR, v, Dis[current_tid][xnei].second, xnei);
							mtx_5952[v].unlock();
						}
					}
				}
			}


			for (int i : Dis_changed) {
				Dis[current_tid][i] = { -1, -1 };
			}

			mtx_595_1.lock();
			Qid_595.push(current_tid);
			mtx_595_1.unlock();

			return 1; }));
	}

	for (auto&& result : results) {
		result.get();
	}
	results.clear();
}

void WeightDecreaseMaintenance_improv(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_new, int thread_num) {

	//auto begin = std::chrono::high_resolution_clock::now();

	std::vector<affected_label> CL;
	//WeightDecreaseMaintenance_improv_step1(v1, v2, w_new, mm.L, mm.PPR, CL, thread_num);
	WeightDecreaseMaintenance_improv_step1_parallelVersion(v1, v2, w_new, &mm.L, &mm.PPR, &CL, thread_num);

	//auto end = std::chrono::high_resolution_clock::now();
	//double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	//cout << "runningtime1 = " << runningtime << "s" << endl;
	//begin = std::chrono::high_resolution_clock::now();

	DIFFUSE(&instance_graph, &mm.L, &mm.PPR, CL, thread_num);

	//end = std::chrono::high_resolution_clock::now();
	//runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	//cout << "runningtime2 = " << runningtime << "s" << endl;
}
