#pragma once

#include <set>
#include <build_in_progress/HL/dynamic/graph_hash_of_mixed_weighted_PLL_dynamic.h>

#define Q(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(mm.L, x, y) // reduction is not used here
#define Q2(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(mm.L, x, y) // reduction is not used here

class affected_label {
public:
	int first, second;
	weightTYPE dis;
	affected_label(int _first, int _second, weightTYPE _dis) {
		first = _first;
		second = _second;
		dis = _dis;
	}
};

class pair_label { // pair_label is stored in NoP
public:
	int first, second;
	pair_label(int _first, int _second) {
		first = _first;
		second = _second;
	}
	bool operator < (const pair_label other) const { // used to sort pair_label in set
		if (first != other.first)
			return first < other.first;
		return second < other.second;
	}
	bool operator == (const pair_label other) const { // used to search pair_label in NoP
		if (first != other.first)
			return false;
		return second == other.second;
	}
};


void PI11(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	std::vector<affected_label>& al1_curr, std::vector<affected_label>& al1_next, weightTYPE w_change) {

	//weightTYPE d_nei_max = std::numeric_limits<weightTYPE>::max();
	for (auto it = al1_curr.begin(); it != al1_curr.end(); it++) {
		auto v_neis = instance_graph.adj_v_and_ec(it->first);
		for (auto nei = v_neis.begin(); nei != v_neis.end(); nei++) {
			//if (search_sorted_two_hop_label(mm.L[nei->first], it->second) + nei->second < d_nei_max) {
			//	d_nei_max = search_sorted_two_hop_label(mm.L[nei->first], it->second) + nei->second;
			//}
			if (abs(it->dis + nei->second - search_sorted_two_hop_label(mm.L[nei->first], it->second)) < 1e-5) {
				al1_next.push_back(affected_label(nei->first, it->second, it->dis + nei->second));
			}
		}
		insert_sorted_two_hop_label(mm.L[it->first], it->second, it->dis + w_change);
		//insert_sorted_two_hop_label(mm.L[it->first], it->second, min(it->dis + w_change, d_nei_max));
	}

}

void PI12(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	std::vector<affected_label>& al1_curr, std::vector<pair_label>& al2_next, std::set<pair_label>& nop) {

	for (auto it = al1_curr.begin(); it != al1_curr.end(); it++) {
		std::vector<int> temp = PPR_retrieve(mm.PPR, it->first, it->second);

		temp.push_back(it->second);


		for (auto t = temp.begin(); t != temp.end(); t++) {

			


			if (it->first < *t) { //  && !nop.count(pair_label(*t, it->first))
				weightTYPE d1 = std::numeric_limits<weightTYPE>::max();
				auto neis = instance_graph.adj_v_and_ec(*t);
				for (auto nei = neis.begin(); nei != neis.end(); nei++) {
					d1 = min(d1, search_sorted_two_hop_label(mm.L[nei->first], it->first) + (weightTYPE)nei->second);
				}
				auto query_result = Q2(*t, it->first);
				if (query_result.first > d1) {
					insert_sorted_two_hop_label(mm.L[*t], it->first, d1);
					al2_next.push_back(pair_label(*t, it->first));
					//nop.insert(pair_label(*t, it->first));
					//PPR_erase(mm.PPR, it->first, it->second, *t); // update PPR
				}
				else {
					PPR_insert(mm.PPR, *t, query_result.second, it->first);
					PPR_insert(mm.PPR, it->first, query_result.second, *t);
				}
			}
			if (*t < it->first) { //  && !nop.count(pair_label(it->first, *t))
				weightTYPE d1 = std::numeric_limits<weightTYPE>::max();
				auto neis = instance_graph.adj_v_and_ec(it->first);
				for (auto nei = neis.begin(); nei != neis.end(); nei++) {
					d1 = min(d1, search_sorted_two_hop_label(mm.L[nei->first], *t) + (weightTYPE)nei->second);
				}
				auto query_result = Q2(it->first, *t);
				if (query_result.first > d1) {
					insert_sorted_two_hop_label(mm.L[it->first], *t, d1);
					al2_next.push_back(pair_label(it->first, *t));
					//nop.insert(pair_label(it->first, *t));
					//PPR_erase(mm.PPR, *t, it->second, it->first); // update PPR
				}
				else {
					PPR_insert(mm.PPR, *t, query_result.second, it->first);
					PPR_insert(mm.PPR, it->first, query_result.second, *t);
				}
			}

			//if (nop.count(pair_label(*t, it->first)) == 1) {
			//	cout << "here" << endl;
			//}

			//if (nop.count(pair_label(it->first, *t)) == 1) {
			//	cout << "here" << endl;
			//}
		}
	}

}

void PI22(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	std::vector<pair_label>& al2_curr, std::vector<pair_label>& al2_next, std::set<pair_label>& nop) {

	for (auto it = al2_curr.begin(); it != al2_curr.end(); it++) {
		auto v_neis = instance_graph.adj_v_and_ec(it->first);
		for (auto nei = v_neis.begin(); nei != v_neis.end(); nei++) {

			if (nei->first > it->second) {
				weightTYPE search_result = search_sorted_two_hop_label(mm.L[it->first], it->second) + nei->second;
				auto query_result = Q2(nei->first, it->second);
				if (query_result.first + 1e-5 >= search_result) {
					insert_sorted_two_hop_label(mm.L[nei->first], it->second, search_result);
					al2_next.push_back(pair_label(nei->first, it->second));
					//nop.insert(pair_label(nei->first, it->second));
				}
				else {
					PPR_insert(mm.PPR, nei->first, query_result.second, it->second);
					PPR_insert(mm.PPR, it->second, query_result.second, nei->first);
				}
			}
		}
	}

}


void WeightIncreaseMaintenance(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new) {


	//cout << "xx " << Q(v1, v2) << " " << w_old << endl;

	//if (Q(v1, v2) + 1e-5 < w_old)
	//	return;

	std::set<pair_label> nop;
	std::vector<affected_label> al1_curr, al1_next;
	std::vector<pair_label> al2_curr, al2_next;

	for (auto it = mm.L[v1].begin(); it != mm.L[v1].end(); it++) {
		if (it->vertex <= v2 && abs(search_sorted_two_hop_label(mm.L[v2], it->vertex) - it->distance - w_old) < 1e-5) {
			al1_curr.push_back(affected_label(v2, it->vertex, it->distance + w_old));
		}
	}

	for (auto it = mm.L[v2].begin(); it != mm.L[v2].end(); it++) {
		if (it->vertex <= v1 && abs(search_sorted_two_hop_label(mm.L[v1], it->vertex) - it->distance - w_old) < 1e-5) {
			al1_curr.push_back(affected_label(v1, it->vertex, it->distance + w_old));
		}
	}

	while (al1_curr.size() || al2_curr.size()) {
		PI11(instance_graph, mm, al1_curr, al1_next, w_new - w_old);
		PI12(instance_graph, mm, al1_curr, al2_next, nop);
		PI22(instance_graph, mm, al2_curr, al2_next, nop);
		al1_curr = al1_next;
		al2_curr = al2_next;
		std::vector<affected_label>().swap(al1_next);
		std::vector<pair_label>().swap(al2_next);
	}
}

