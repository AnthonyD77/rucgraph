#pragma once

#include <build_in_progress/HL/dynamic/graph_hash_of_mixed_weighted_PLL_dynamic.h>

/*this function has bugs*/
void WeightDecreaseMaintenance_step1(int v1, int v2, weightTYPE w_new, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR, std::vector<affected_label>* CL,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	for (int sl = 0; sl < 2; sl++) {
		if (sl == 1) {
			swap(v1, v2);
		}
		for (auto it : (*L)[v1]) {
			if (it.vertex <= v2) {
				results_dynamic.emplace_back(pool_dynamic.enqueue([it, v2, L, PPR, w_new, CL] {

					pair<weightTYPE, int> query_result;
					if (it.vertex == v2) {
						query_result = { 0,v2 };
					}
					else {
						mtx_595[it.vertex].lock(), mtx_595[v2].lock(); // you cannot lock the same lock twice!
						query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, it.vertex, v2); // query_result is {distance, common hub}
						mtx_595[it.vertex].unlock(), mtx_595[v2].unlock();
					}

					if (query_result.first > it.distance + w_new) {
						mtx_595[v2].lock();
						insert_sorted_two_hop_label((*L)[v2], it.vertex, it.distance + w_new);
						mtx_595[v2].unlock();
						mtx_595_1.lock();
						CL->push_back(affected_label{ v2 , it.vertex, it.distance + w_new });
						mtx_595_1.unlock();
					}
					else {
						mtx_595[v2].lock();
						auto search_result = search_sorted_two_hop_label2((*L)[v2], it.vertex);
						mtx_595[v2].unlock();
						if (search_result.first > it.distance + w_new && search_result.first != MAX_VALUE) {
							mtx_595[v2].lock();
							(*L)[v2][search_result.second].distance = it.distance + w_new;
							mtx_595[v2].unlock();
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

/*this function has bugs*/
void ProDecrease_parallel(graph_hash_of_mixed_weighted* instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR,
	std::vector<affected_label>& CL_curr, std::vector<affected_label>* CL_next, ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	for (auto it = CL_curr.begin(); it != CL_curr.end(); it++) {
		results_dynamic.emplace_back(pool_dynamic.enqueue([it, instance_graph, CL_next, L, PPR] {

			auto neis = instance_graph->adj_v_and_ec(it->first);
			for (auto nei = neis.begin(); nei != neis.end(); nei++) {
				if (it->second < nei->first) {

					mtx_595[nei->first].lock(), mtx_595[it->second].lock();
					auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, nei->first, it->second); // query_result is {distance, common hub}
					mtx_595[nei->first].unlock(), mtx_595[it->second].unlock();

					if (query_result.first > it->dis + nei->second) {

						mtx_595[nei->first].lock();
						insert_sorted_two_hop_label((*L)[nei->first], it->second, it->dis + nei->second);
						mtx_595[nei->first].unlock();
						mtx_595_1.lock();
						CL_next->push_back(affected_label(nei->first, it->second, it->dis + nei->second));
						mtx_595_1.unlock();
					}
					else {
						mtx_595[nei->first].lock();
						auto search_result = search_sorted_two_hop_label2((*L)[nei->first], it->second);
						mtx_595[nei->first].unlock();
						if (search_result.first != std::numeric_limits<weightTYPE>::max() && search_result.first > it->dis + nei->second) {
							mtx_595[nei->first].lock();
							(*L)[nei->first][search_result.second].distance = it->dis + nei->second;
							mtx_595[nei->first].unlock();
							mtx_595_1.lock();
							CL_next->push_back(affected_label(nei->first, it->second, it->dis + nei->second));
							mtx_595_1.unlock();
						}
						if (query_result.second != it->second) {
							mtx_5952[nei->first].lock();
							PPR_insert(*PPR, nei->first, query_result.second, it->second);
							mtx_5952[nei->first].unlock();
						}
						if (query_result.second != nei->first) {
							mtx_5952[it->second].lock();
							PPR_insert(*PPR, it->second, query_result.second, nei->first);
							mtx_5952[it->second].unlock();
						}
					}
				}
			}

			return 1; }));
	}
}


void ProDecrease(graph_hash_of_mixed_weighted& instance_graph, vector<vector<two_hop_label_v1>>& L, PPR_type& PPR,
	std::vector<affected_label>& CL_curr, std::vector<affected_label>& CL_next) {

	affected_label xx;

	for (auto it = CL_curr.begin(); it != CL_curr.end(); it++) {
		auto neis = instance_graph.adj_v_and_ec(it->first);
		for (auto nei = neis.begin(); nei != neis.end(); nei++) {
			if (it->second < nei->first) {
				auto query_result = Query2(nei->first, it->second); // query_result is {distance, common hub}
				if (query_result.first > it->dis + nei->second) {
					insert_sorted_two_hop_label(L[nei->first], it->second, it->dis + nei->second);
					xx.first = nei->first, xx.second = it->second, xx.dis = it->dis + nei->second;
					CL_next.push_back(xx);
				}
				else {
					auto search_result = search_sorted_two_hop_label2(L[nei->first], it->second);
					if (search_result.first != std::numeric_limits<weightTYPE>::max() && search_result.first > it->dis + nei->second) {
						L[nei->first][search_result.second].distance = it->dis + nei->second;
						xx.first = nei->first, xx.second = it->second, xx.dis = it->dis + nei->second;
						CL_next.push_back(xx);
					}
					PPR_insert(PPR, nei->first, query_result.second, it->second);
					PPR_insert(PPR, it->second, query_result.second, nei->first);
				}
			}
		}
	}
}



void WeightDecreaseMaintenance(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_new,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	std::vector<affected_label> CL_curr, CL_next;
	affected_label xx;

	//WeightDecreaseMaintenance_step1(v1, v2, w_new, &mm.L, &mm.PPR, &CL_curr, pool_dynamic, results_dynamic);


	/*
	the following part does not suit parallel computation:
	the reason is that L is changed below, and as a result, in each following loop, L[v2] or L[v1] is locked at each step, 
	which means that following loops cannot be actually parallized
	*/
	for (auto it = mm.L[v1].begin(); it != mm.L[v1].end(); it++) {	
		if (it->vertex <= v2) {
			auto query_result = Q2(it->vertex, v2); // query_result is {distance, common hub}
			if (query_result.first > it->distance + w_new) {
				insert_sorted_two_hop_label(mm.L[v2], it->vertex, it->distance + w_new);
				xx.first = v2, xx.second = it->vertex, xx.dis = it->distance + w_new;
				CL_curr.push_back(xx);
			}
			else {
				auto search_result = search_sorted_two_hop_label2(mm.L[v2], it->vertex);
				if (search_result.first != std::numeric_limits<weightTYPE>::max() && search_result.first > it->distance + w_new) {
					mm.L[v2][search_result.second].distance = it->distance + w_new;
					xx.first = v2, xx.second = it->vertex, xx.dis = it->distance + w_new;
					CL_next.push_back(xx);
				}
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
				xx.first = v1, xx.second = it->vertex, xx.dis = it->distance + w_new;
				CL_curr.push_back(xx);
			}
			else {
				auto search_result = search_sorted_two_hop_label2(mm.L[v1], it->vertex);
				if (search_result.first != std::numeric_limits<weightTYPE>::max() && search_result.first > it->distance + w_new) {
					mm.L[v1][search_result.second].distance = it->distance + w_new;
					xx.first = v1, xx.second = it->vertex, xx.dis = it->distance + w_new;
					CL_next.push_back(xx);
				}
				PPR_insert(mm.PPR, v1, query_result.second, it->vertex);
				PPR_insert(mm.PPR, it->vertex, query_result.second, v1);
			}
		}
	}




	while (CL_curr.size()) {
		//ProDecrease(instance_graph, mm.L, mm.PPR, CL_curr, CL_next);
		ProDecrease_parallel(&instance_graph, &mm.L, &mm.PPR, CL_curr, &CL_next, pool_dynamic, results_dynamic);
		CL_curr = CL_next;
		std::vector<affected_label>().swap(CL_next);
	}
}


