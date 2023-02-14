#pragma once

#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_binary_operations.h>

/*the pairs should be sorted from small to large based on pair.first, for fast querying*/
typedef std::vector<std::vector<std::pair<int, std::vector<int>>>> PPR_type;




void PPR_insert(PPR_type& PPR, int v1, int v2, int v3) {

	/*add v3 into PPR(v1, v2)*/

	int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(PPR[v1], v2);
	if (pos == -1) {
		std::vector<int> x = { v3 };
		graph_hash_of_mixed_weighted_binary_operations_insert(PPR[v1], v2, x);
	}
	else {
		PPR[v1][pos].second.push_back(v3);
	}

}


void PPR_replace(PPR_type& PPR, int v1, int v2, std::vector<int>& loads) {

	/*replace PPR(v1, v2) = loads*/

	int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(PPR[v1], v2);
	if (pos == -1) {
		graph_hash_of_mixed_weighted_binary_operations_insert(PPR[v1], v2, loads);
	}
	else {
		PPR[v1][pos].second = loads;
	}

}













