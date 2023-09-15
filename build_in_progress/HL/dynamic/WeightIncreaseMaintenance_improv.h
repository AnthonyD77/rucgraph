#pragma once

#include <queue>
#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

void SPREAD1(graph_hash_of_mixed_weighted* instance_graph, vector<vector<two_hop_label_v1>>* L,
	std::vector<affected_label>& al1, std::vector<pair_label>* al2) {

	for (auto it : al1) {
		queue<pair<int, weightTYPE> > q; //(u,d)
		int v = it.second;
		q.push(pair<int, weightTYPE>(it.first, it.dis));
		while (!q.empty()) {
			int x = q.front().first;
			weightTYPE dx = q.front().second;
			q.pop();
			insert_sorted_two_hop_label((*L)[x], v, MAX_VALUE); // this does not change the size of L[x] here, so does not need to lock here
			al2->push_back(pair_label(x, v));
			auto v_neis = instance_graph->adj_v_and_ec(x);
			for (auto nei : v_neis) {
				weightTYPE search_weight = search_sorted_two_hop_label((*L)[nei.first], v);
				if (abs(dx + nei.second - search_weight) < 1e-5) { 
					q.push(pair<int, weightTYPE>(nei.first, dx + nei.second));
				}
			}
		}
	}
}

void SPREAD2(graph_hash_of_mixed_weighted* instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR,
	std::vector<pair_label>& al2, std::vector<affected_label>* al3) {

	for (auto it : al2) {
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
				if (query_result.first > d1 + 1e-5) { // only add new label when it's absolutely necessary
					al3->push_back(affected_label(t, v, d1));
					//cout<<"spread21: "<<t<<' '<<v << ' ' << d1 << ' ' << query_result.first << endl;
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
			else if (t < v) {
				weightTYPE d1 = MAX_VALUE;
				auto neis = instance_graph->adj_v_and_ec(v);
				for (auto nei : neis) {
					d1 = min(d1, search_sorted_two_hop_label((*L)[nei.first], t) + (weightTYPE)nei.second);
				}
				auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, v, t);
				if (query_result.first > d1 + 1e-5) {
					al3->push_back(affected_label(v, t, d1));
					//cout << "spread22: " << v << ' ' << t << ' ' << d1 << ' ' << query_result.first << endl;
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

void SPREAD3(graph_hash_of_mixed_weighted* instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR, std::vector<affected_label>& al3) {

	for (auto it : al3) {

		int u = it.first, v = it.second;
		weightTYPE du = it.dis;

		auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, u, v);

		if (query_result.first < du + 1e-5) {
			if (query_result.second != v) {
				PPR_insert(*PPR, u, query_result.second, v);
			}
			if (query_result.second != u) {
				PPR_insert(*PPR, v, query_result.second, u);
			}
			continue;
		}

		vector<int> Dis_changed;
		int N=instance_graph->hash_of_vectors.size();
		vector<pair<weightTYPE,int> > DIS(N,{-1,-1});
		vector<handle_t_for_DIFFUSE> Q_HANDLES(N);
		vector<weightTYPE> Q_VALUE(N,MAX_VALUE);
		DIS[u]={du,u}; // <distance, hub responsible for this distance>
		Dis_changed.push_back(u);

		boost::heap::fibonacci_heap<node_for_DIFFUSE> pq;
		Q_HANDLES[u] = pq.push(node_for_DIFFUSE(u, du));
		Q_VALUE[u] = du;

		//cout << "spread3 0: " << u << ' ' << v << " " << du << endl;
		while (!pq.empty()) {
			int x = pq.top().index;
			weightTYPE dx = pq.top().disx;
			pq.pop();
			Q_VALUE[x] = MAX_VALUE;

			insert_sorted_two_hop_label((*L)[x], v, std::min(dx, search_sorted_two_hop_label((*L)[x], v)));
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
							Q_HANDLES[xnei] = pq.push(node_for_DIFFUSE(xnei, d_new));
						}
						else {
							pq.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
						}
						Q_VALUE[xnei] = d_new;
					}
					else {
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

void WeightIncreaseMaintenance_improv(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_old) {

	std::vector<affected_label> al1, al3;
	std::vector<pair_label> al2;

	/*it's slow to paralize the following part*/
	for (auto it : mm.L[v1]) {
		if (it.vertex <= v2 && abs(search_sorted_two_hop_label(mm.L[v2], it.vertex) - it.distance - w_old) < 1e-5) {
			al1.push_back(affected_label(v2, it.vertex, it.distance + w_old));
		}
	}
	for (auto it : mm.L[v2]) {
		if (it.vertex <= v1 && abs(search_sorted_two_hop_label(mm.L[v1], it.vertex) - it.distance - w_old) < 1e-5) {
			al1.push_back(affected_label(v1, it.vertex, it.distance + w_old));
		}
	}

	//cout << "al1.size() " << al1.size() << endl;

	SPREAD1(&instance_graph, &mm.L, al1, &al2);
	SPREAD2(&instance_graph, &mm.L, &mm.PPR, al2, &al3);
	SPREAD3(&instance_graph, &mm.L, &mm.PPR, al3);

	//for (auto it : al2) {
	//	cout << "al2 " << it.first << " " << it.second << endl;
	//}
	//for (auto it : al3) {
	//	cout << "al3 " << it.first << " " << it.second << " " << it.dis << endl;
	//}
}

