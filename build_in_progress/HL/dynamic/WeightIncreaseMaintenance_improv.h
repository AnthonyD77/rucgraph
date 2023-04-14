#pragma once

#include <queue>
#include <unordered_map>
#include <build_in_progress/HL/dynamic/graph_hash_of_mixed_weighted_PLL_dynamic.h>


class pq_label {
public:
	int u;
	weightTYPE dis;
	pq_label(int _u, weightTYPE _dis) {
		u = _u;
		dis = _dis;
	}
	bool operator < (const pq_label other)const {
		return dis < other.dis;
	}
};
typedef typename boost::heap::fibonacci_heap<pq_label>::handle_type handle_for_SPREAD3; // pairing heap has a similar speed with fibonacci_heap here


void SPREAD1(graph_hash_of_mixed_weighted& instance_graph, vector<vector<two_hop_label_v1>>& L,
	std::vector<affected_label>& al1, std::vector<pair_label>& al2, weightTYPE w_change) {

	for (auto it : al1) {
		queue<pair<int, weightTYPE> > q; //(u,d)
		int v = it.second;
		q.push(pair<int, weightTYPE>(it.first, it.dis));
		while (!q.empty()) {
			pair<int, weightTYPE> fr = q.front();
			int x = fr.first;
			weightTYPE dx = fr.second;
			q.pop();
			insert_sorted_two_hop_label(L[x], v, MAX_VALUE);
			//cout<<"spread1: "<<v<<' '<<x<<' '<<dx+w_change<<endl;
			al2.push_back(pair_label(x, v));
			//cout<<"al2: "<<x<<' '<<v<<endl;
			auto v_neis = instance_graph.adj_v_and_ec(x);
			for (auto nei : v_neis) {
				pair_label xnv(nei.first, v);
				weightTYPE search_weight = search_sorted_two_hop_label(L[nei.first], v);
				if (find(al2.begin(), al2.end(), xnv) == al2.end() && abs(dx + nei.second - search_weight) < 1e-5) { // al2 find is extremely slow here! (maybe just remove it?)
					q.push(pair<int, weightTYPE>(nei.first, dx + nei.second));
				}
			}
		}
	}

}

void SPREAD2(graph_hash_of_mixed_weighted& instance_graph, vector<vector<two_hop_label_v1>>& L, PPR_type& PPR,
	std::vector<pair_label>& al2, std::vector<affected_label>& al3) {

	for (auto it : al2) {
		int v = it.first, u = it.second;
		std::vector<int> temp = PPR_retrieve(PPR, v, u);
		PPR_binary_operations_insert(temp, u);
		//cout<<"sp2: "<<it->first<<' '<<it->second<<endl;

		for (auto t : temp) {
			if (v < t) {
				weightTYPE d1 = std::numeric_limits<weightTYPE>::max();
				auto neis = instance_graph.adj_v_and_ec(t);
				for (auto nei : neis) {
					d1 = min(d1, search_sorted_two_hop_label(L[nei.first], v) + (weightTYPE)nei.second);
				}
				auto query_result = Query2(t, v);
				if (query_result.first > d1 + 1e-5) { // only add new label when it's absolutely necessary
					al3.push_back(affected_label(t, v, d1));
					//cout<<"spread21: "<<t<<' '<<v << ' ' << d1 << ' ' << query_result.first << endl;
				}
				else {
					PPR_insert(PPR, t, query_result.second, v);
					PPR_insert(PPR, v, query_result.second, t);
				}
			}
			if (t < v) {
				weightTYPE d1 = std::numeric_limits<weightTYPE>::max();
				auto neis = instance_graph.adj_v_and_ec(v);
				for (auto nei : neis) {
					d1 = min(d1, search_sorted_two_hop_label(L[nei.first], t) + (weightTYPE)nei.second);
				}
				auto query_result = Query2(v, t);
				if (query_result.first > d1 + 1e-5) {
					al3.push_back(affected_label(v, t, d1));
					//cout << "spread22: " << v << ' ' << t << ' ' << d1 << ' ' << query_result.first << endl;
				}
				else {
					PPR_insert(PPR, t, query_result.second, v);
					PPR_insert(PPR, v, query_result.second, t);
				}
			}
		}
	}
}

void SPREAD3(graph_hash_of_mixed_weighted& instance_graph, vector<vector<two_hop_label_v1>>& L, PPR_type& PPR, std::vector<affected_label>& al3) {

	for (auto it : al3) {
		int u = it.first, v = it.second;
		weightTYPE du = it.dis;
		auto query_result = Query2(u, v);
		if (query_result.first < du + 1e-5) {
			PPR_insert(PPR, u, query_result.second, v);
			PPR_insert(PPR, v, query_result.second, u);
			continue;
		}
		unordered_map<int, pair<weightTYPE, int> > dis; //pair<distance,hub>
		dis[u] = pair(du, v);
		boost::heap::fibonacci_heap<pq_label> pq;
		unordered_map <int, handle_for_SPREAD3> handles;
		handles[u] = pq.push(pq_label(u, du));
		//cout << "spread3 0: " << u << ' ' << v << " " << du << endl;

		while (!pq.empty()) {
			int x = pq.top().u;
			weightTYPE dx = pq.top().dis;
			pq.pop();
			handles.erase(x);
			insert_sorted_two_hop_label(L[x], v, dx);
			//cout << "spread3 1: " << x << ' ' << v << ' ' << dx << endl;
			auto neis = instance_graph.adj_v_and_ec(x);
			for (auto nei : neis) {
				if (v < nei.first) {
					if (dis.count(nei.first) == 0) {
						auto query_result = Query2(nei.first, v);
						dis[nei.first] = pair(query_result.first, query_result.second);

					}
					if (dis[nei.first].first > dx + nei.second - 1e-5) {
						dis[nei.first] = pair(dx + nei.second, v);
						//pq.push(pq_label(nei.first, dx + nei.second));
						if (handles.count(nei.first) == 0) {
							handles[nei.first] = pq.push(pq_label(nei.first, dx + nei.second));
						}
						else {
							pq.update(handles[nei.first], pq_label(nei.first, dx + nei.second));
						}
					}
					else {
						PPR_insert(PPR, nei.first, dis[nei.first].second, v);
						PPR_insert(PPR, v, dis[nei.first].second, nei.first);
					}
				}
			}
		}
	}

}

void WeightIncreaseMaintenance_improv(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new) {

	std::vector<affected_label> al1, al3;
	std::vector<pair_label> al2;

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

	SPREAD1(instance_graph, mm.L, al1, al2, w_new - w_old);
	SPREAD2(instance_graph, mm.L, mm.PPR, al2, al3);
	SPREAD3(instance_graph, mm.L, mm.PPR, al3);
}

