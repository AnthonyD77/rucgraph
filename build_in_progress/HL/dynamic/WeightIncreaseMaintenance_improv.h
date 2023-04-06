#pragma once

#include <set>
#include <queue>
#include <map>
#include <boost/heap/fibonacci_heap.hpp> 
#include <build_in_progress/HL/dynamic/graph_hash_of_mixed_weighted_PLL_dynamic.h>

#define Q(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(mm.L, x, y) // reduction is not used here
#define Q2(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(mm.L, x, y) // reduction is not used here

class affected_label2 {
public:
	int first, second;
	weightTYPE dis;
	affected_label2(int _first, int _second, weightTYPE _dis) {
		first = _first;
		second = _second;
		dis = _dis;
	}
};

class pair_label2 { // pair_label2 is stored in NoP
public:
	int first, second;
	pair_label2(int _first, int _second) {
		first = _first;
		second = _second;
	}
	bool operator == (const pair_label2 other) const {
		return (first == other.first && second == other.second);
	}
	bool operator < (const pair_label2 other) const { // used to sort/search pair_label2 in set
		if (first != other.first)
			return first < other.first;
		return second < other.second;
	}
};

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



void SPREAD1(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	std::vector<affected_label2>& al1, std::vector<pair_label2>& al2, weightTYPE w_change, std::set<pair_label2>& nop) {

	for (auto it = al1.begin(); it != al1.end(); it++) {
		queue<pair<int, weightTYPE> > q; //(u,d)
		int v = it->second;
		q.push(pair<int, weightTYPE>(it->first, it->dis));
		while (!q.empty()) {
			pair<int, weightTYPE> fr = q.front();
			int x = fr.first;
			weightTYPE dx = fr.second;
			q.pop();
			insert_sorted_two_hop_label(mm.L[x], v, dx + w_change);
			//cout<<"spread1: "<<v<<' '<<x<<' '<<dx+w_change<<endl;
			al2.push_back(pair_label2(x, v));
			//cout<<"al2: "<<x<<' '<<v<<endl;
			auto v_neis = instance_graph.adj_v_and_ec(x);
			for (auto nei = v_neis.begin(); nei != v_neis.end(); nei++) {
				pair_label2 xnv(nei->first, v);
				weightTYPE search_weight = search_sorted_two_hop_label(mm.L[nei->first], v);
				if (find(al2.begin(), al2.end(), xnv) == al2.end() && abs(dx + nei->second - search_weight) < 1e-5) {
					q.push(pair<int, weightTYPE>(nei->first, dx + nei->second));
				}
			}
		}
	}

}

void SPREAD2(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	std::vector<pair_label2>& al2, std::vector<affected_label2>& al3, std::set<pair_label2>& nop) {
	for (auto it = al2.begin(); it != al2.end(); it++) {
		std::vector<int> temp = PPR_retrieve(mm.PPR, it->first, it->second);
		PPR_binary_operations_insert(temp, it->second);
		//cout<<"sp2: "<<it->first<<' '<<it->second<<endl;

		for (auto t = temp.begin(); t != temp.end(); t++) {
			if (it->first < *t && (!nop.count(pair_label2(*t, it->first)))) {
				weightTYPE d1 = std::numeric_limits<weightTYPE>::max();
				auto neis = instance_graph.adj_v_and_ec(*t);
				for (auto nei = neis.begin(); nei != neis.end(); nei++) {
					d1 = min(d1, search_sorted_two_hop_label(mm.L[nei->first], it->first) + (weightTYPE)nei->second);
				}
				auto query_result = Q2(*t, it->first);
				if (query_result.first - d1 > -1e-5) {
					al3.push_back(affected_label2(*t, it->first, d1));
					//cout<<"spread21: "<<*t<<' '<<it->first<<endl;
					nop.insert(pair_label2(*t, it->first));
				}
			}
			else if (*t < it->first && (!nop.count(pair_label2(it->first, *t)))) {
				weightTYPE d1 = std::numeric_limits<weightTYPE>::max();
				auto neis = instance_graph.adj_v_and_ec(it->first);
				for (auto nei = neis.begin(); nei != neis.end(); nei++) {
					d1 = min(d1, search_sorted_two_hop_label(mm.L[nei->first], *t) + (weightTYPE)nei->second);
				}
				auto query_result = Q2(it->first, *t);
				if (query_result.first - d1 > -1e-5) {
					al3.push_back(affected_label2(it->first, *t, d1));
					//cout<<"spread22: "<<it->first<<' '<<*t<<endl;
					nop.insert(pair_label2(it->first, *t));
				}
			}
		}
	}
}

void SPREAD3(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, std::vector<affected_label2>& al3, std::set<pair_label2>& nop) {

	for (auto it = al3.begin(); it != al3.end(); it++) {
		map<int, pair<weightTYPE, int> > dis; //pair<distance,hub>
		dis[it->first] = pair(it->dis, it->second);
		boost::heap::fibonacci_heap<pq_label> pq;
		pq.push(pq_label(it->first, it->dis));
		int v = it->second;
		while (!pq.empty()) {
			int x = pq.top().u;
			weightTYPE dx = pq.top().dis;
			insert_sorted_two_hop_label(mm.L[x], v, dx);
			//cout<<"spread3: "<<x<<' '<<v<<' '<<dx<<endl;
			pq.pop();
			auto neis = instance_graph.adj_v_and_ec(x);
			for (auto nei = neis.begin(); nei != neis.end(); nei++) {
				if (v < nei->first) {
					//cout<<"spread nei: "<<nei->first<<' '<<v<<endl;
					if (!dis.count(nei->first)) {
						auto query_result = Q2(nei->first, v);
						dis[nei->first] = pair(query_result.first, query_result.second);
						//cout<<"1: "<<query_result.first<<' '<<query_result.second<<endl;
					}
					if (dis[nei->first].first - dx - nei->second > -1e-5) {
						dis[nei->first] = pair(dx + nei->second, v);
						pq.push(pq_label(nei->first, dx + nei->second));
						//cout<<"2: "<<dx+nei->second<<' '<<v<<endl;
					}
					else {
						PPR_insert(mm.PPR, nei->first, dis[nei->first].second, v);
						PPR_insert(mm.PPR, v, dis[nei->first].second, nei->first);
					}
				}
			}
		}
	}

}

void WeightIncreaseMaintenance_improv(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new) {

	std::set<pair_label2> nop;
	std::vector<affected_label2> al1, al3;
	std::vector<pair_label2> al2;

	for (auto it : mm.L[v1]) {
		if (it.vertex <= v2 && abs(search_sorted_two_hop_label(mm.L[v2], it.vertex) - it.distance - w_old) < 1e-5) {
			al1.push_back(affected_label2(v2, it.vertex, it.distance + w_old));
			//cout<<"main1: "<<v2<<' '<<it->vertex<<endl;
		}
	}
	for (auto it : mm.L[v2]) {
		if (it.vertex <= v1 && abs(search_sorted_two_hop_label(mm.L[v1], it.vertex) - it.distance - w_old) < 1e-5) {
			al1.push_back(affected_label2(v1, it.vertex, it.distance + w_old));
			//cout<<"main2: "<<v1<<' '<<it->vertex<<endl;
		}
	}

	SPREAD1(instance_graph, mm, al1, al2, w_new - w_old, nop);
	SPREAD2(instance_graph, mm, al2, al3, nop);
	SPREAD3(instance_graph, mm, al3, nop);
}

