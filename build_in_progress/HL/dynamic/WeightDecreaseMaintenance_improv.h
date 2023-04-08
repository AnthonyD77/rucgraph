#pragma once

#include <build_in_progress/HL/dynamic/graph_hash_of_mixed_weighted_PLL_dynamic.h>
#include <queue>
#include <boost/heap/fibonacci_heap.hpp> 
#include <unordered_map>
#include <vector>

//#define Q31(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(mm.L, mm.reduction_measures_2019R2, mm.reduction_measures_2019R1,mm.f_2019R1, instance_graph, x, y)
#define Q31(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, x, y) // reduction is not used here
#define Q32(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(L, x, y) // reduction is not used here

class changed_label3 {
public:
	int left, right;
	weightTYPE dis;
};

class dij {
public:
	int x;
	weightTYPE disx;
};
struct cmp {
	bool operator()(dij& a, dij& b) {
		return a.disx < b.disx;
	}
};
priority_queue <dij, std::vector<dij>, cmp> q; //Min Heap PQ


// dij initialization
struct node_for_DIFFUSE {
	int index;
	weightTYPE disx;
	node_for_DIFFUSE() {
	}
	node_for_DIFFUSE(int _u, weightTYPE _dis) {
		index = _u;
		disx = _dis;
	}
}; // define the node in the queue

bool operator<(node_for_DIFFUSE const& x, node_for_DIFFUSE const& y) {
	return x.disx > y.disx; // < is the max-heap; > is the min heap
}
typedef typename boost::heap::fibonacci_heap<node_for_DIFFUSE>::handle_type handle_t_for_DIFFUSE;

weightTYPE searchPQvalue(boost::heap::fibonacci_heap<node_for_DIFFUSE>Q, int index) {
	for (auto it = Q.begin(); it != Q.end(); it++) {
		if (it->index == index) {
			return it->disx;
		}
	}
}

void DIFFUSE(graph_hash_of_mixed_weighted& instance_graph, vector<vector<two_hop_label_v1>>& L, PPR_type& PPR, std::vector<changed_label3>& CL) {

	changed_label3 xx;
	node_for_DIFFUSE temp1;

	for (auto it : CL) {
		int u = it.left, v = it.right;
		weightTYPE du = it.dis;

		unordered_map<int, pair<weightTYPE, int>> Dis(L.size());
		Dis[u] = pair(du, v);

		boost::heap::fibonacci_heap<node_for_DIFFUSE> Q;
		unordered_map<int, handle_t_for_DIFFUSE> Q_handles;
		unordered_map<int, weightTYPE> Q_value;

		temp1.index = u;
		temp1.disx = du;
		Q_handles[u] = Q.push(temp1);
		Q_value[u] = du;

		while (!Q.empty()) {

			node_for_DIFFUSE temp2;
			temp2 = Q.top();
			int x = temp2.index;
			weightTYPE dx = temp2.disx;
			Q.pop();
			Q_handles.erase(x);
			insert_sorted_two_hop_label(L[x], v, dx);

			auto neis = instance_graph.adj_v_and_ec(x);
			for (auto nei = neis.begin(); nei != neis.end(); nei++) {
				if (v < nei->first) {
					auto query_result = Q32(nei->first, v);

					if (Dis.count(nei->first) == 0) {
						auto query_result = Q32(nei->first, v);
						Dis[nei->first] = pair(query_result.first, query_result.second);
					}
					if (Dis[nei->first].first > dx + nei->second - 1e-5) {
						Dis[nei->first] = pair(dx + nei->second, v);
						if (Q_handles.count(nei->first) == 0) {
							Q_handles[nei->first] = Q.push(node_for_DIFFUSE(nei->first, dx + nei->second));
						}
						else {
							Q.update(Q_handles[nei->first], node_for_DIFFUSE(nei->first, dx + nei->second));
						}
						Q_value[nei->first] = dx + nei->second;
					}
					else {
						auto search_result = search_sorted_two_hop_label2(L[nei->first], v);
						if (search_result.second != -1 && std::min(search_result.first, Q_value[nei->first]) > dx + nei->second) {
							if (Q_handles.count(nei->first) == 0) {
								Q_handles[nei->first] = Q.push(node_for_DIFFUSE(nei->first, dx + nei->second));
							}
							else {
								Q.update(Q_handles[nei->first], node_for_DIFFUSE(nei->first, dx + nei->second));
							}
							Q_value[nei->first] = dx + nei->second;
						}
						PPR_insert(PPR, nei->first, Dis[nei->first].second, v);
						PPR_insert(PPR, v, Dis[nei->first].second, nei->first);
					}
				}
			}
		}
	}
}

void WeightDecreaseMaintenance_improv(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new) {

	std::vector<changed_label3> CL;
	changed_label3 xx;
	auto& L = mm.L;


	for (auto it : L[v1]) {
		if (it.vertex <= v2) {
			auto query_result = Q32(it.vertex, v2); // query_result is {distance, common hub}
			if (query_result.first > it.distance + w_new) {
				xx.left = v2, xx.right = it.vertex, xx.dis = it.distance + w_new;
				CL.push_back(xx);
			}
			else {
				auto search_result = search_sorted_two_hop_label2(L[v2], it.vertex);
				if (search_result.first > it.distance + w_new && search_result.first != std::numeric_limits<weightTYPE>::max()) {
					xx.left = v2, xx.right = it.vertex, xx.dis = it.distance + w_new;
					CL.push_back(xx);
				}
				PPR_insert(mm.PPR, v2, query_result.second, it.vertex);
				PPR_insert(mm.PPR, it.vertex, query_result.second, v2);
			}
		}
	}
	for (auto it = mm.L[v2].begin(); it != mm.L[v2].end(); it++) {
		if (it->vertex <= v1) {
			auto query_result = Q32(it->vertex, v1); // query_result is {distance, common hub}
			if (query_result.first > it->distance + w_new) {
				xx.left = v1, xx.right = it->vertex, xx.dis = it->distance + w_new;
				CL.push_back(xx);
			}
			else {
				auto search_result = search_sorted_two_hop_label2(mm.L[v1], it->vertex);
				if (search_result.first > it->distance + w_new && search_result.first != std::numeric_limits<weightTYPE>::max()) {
					xx.left = v1, xx.right = it->vertex, xx.dis = it->distance + w_new;
					CL.push_back(xx);
				}
				PPR_insert(mm.PPR, v1, query_result.second, it->vertex);
				PPR_insert(mm.PPR, it->vertex, query_result.second, v1);
			}
		}
	}
	DIFFUSE(instance_graph, mm.L, mm.PPR, CL);
}
