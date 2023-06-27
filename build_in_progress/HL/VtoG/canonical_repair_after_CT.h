#pragma once
#include <build_in_progress/HL/VtoG/graph_hash_of_mixed_weighted_CT.h>



class canonical_repair_after_CT_info {
public:

	long long int label_size_before_canonical_repair = 0;
	long long int label_size_after_canonical_repair = 0;
	double canonical_repair_remove_label_ratio = 0;

	double time_prepare = 0;
	double time_repair = 0;

};


canonical_repair_after_CT_info canonical_repair_after_CT(graph_hash_of_mixed_weighted_CT_v2_case_info& case_info) {

	canonical_repair_after_CT_info return_info;

	auto begin = std::chrono::high_resolution_clock::now();

	auto& isIntree = case_info.isIntree;
	auto& L = case_info.two_hop_case_info.L;
	auto& input_graph = case_info.core_graph;


	//graph_hash_of_mixed_weighted_print(input_graph);

	int N = isIntree.size();
	L_temp_595.resize(N);

	//case_info.two_hop_case_info.print_L();

	vector<pair<int, int>> sorted_vertices;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		sorted_vertices.push_back({ it->first, input_graph.degree(it->first) });
		//cout << "V " << it->first << " " << input_graph.degree(it->first) << endl;
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), compare_pair_second_large_to_small);
	unordered_map<int, int> vertexID_old_to_new;
	vertexID_new_to_old_595.resize(N, 0);
	for (int i = 0; i < N; i++) {
		vertexID_old_to_new[sorted_vertices[i].first] = i;
		vertexID_new_to_old_595[i] = sorted_vertices[i].first;
	}
	vector<pair<int, int>>().swap(sorted_vertices);

	//for (auto x : vertexID_new_to_old_595) {
	//	cout << "vertexID_new_to_old_595 " << x << endl;
	//}

	for (int i = 0; i < N; i++) {
		if (!isIntree[i]) {
			for (auto& l : L[i]) {

				//if (vertexID_old_to_new.count(l.vertex) == 0 || vertexID_old_to_new.count(l.parent_vertex) == 0) {
				//	cout << "here" << endl;
				//}

				l.vertex = vertexID_old_to_new[l.vertex];
				l.parent_vertex = vertexID_old_to_new[l.parent_vertex];
			}
			L_temp_595[vertexID_old_to_new[i]] = L[i];
			vector<two_hop_label_v1>().swap(L[i]);
		}
	}

	reduction_measures_2019R1_new_ID.resize(N, 0);
	reduction_measures_2019R2_new_ID.resize(N, 0);
	f_2019R1_new_ID.resize(N, 0);
	for (int i = 0; i < N; i++) {
		reduction_measures_2019R1_new_ID[vertexID_old_to_new[i]] = case_info.two_hop_case_info.reduction_measures_2019R1[i];
		reduction_measures_2019R2_new_ID[vertexID_old_to_new[i]] = case_info.two_hop_case_info.reduction_measures_2019R2[i];
		f_2019R1_new_ID[vertexID_old_to_new[i]] = vertexID_old_to_new[case_info.two_hop_case_info.f_2019R1[i]];
		sort(L_temp_595[i].begin(), L_temp_595[i].end(), compare_two_hop_label_small_to_large); // sort is necessary
	}
	graph_hash_of_mixed_weighted new_ID_g = graph_hash_of_mixed_weighted_update_vertexIDs(input_graph, vertexID_old_to_new);
	vector <vector<pair<int, double>>>().swap(adjs_new_IDs);
	adjs_new_IDs.resize(N);
	vector<pair<int, double>>().swap(min_adjs_new_IDs);
	min_adjs_new_IDs.resize(N);
	for (auto it = new_ID_g.hash_of_vectors.begin(); it != new_ID_g.hash_of_vectors.end(); it++) {
		adjs_new_IDs[it->first] = new_ID_g.adj_v_and_ec(it->first);
		min_adjs_new_IDs[it->first] = new_ID_g.min_adj(it->first);
	}

	return_info.time_prepare = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s

	auto begin2 = std::chrono::high_resolution_clock::now();
	canonical_repair_multi_threads(new_ID_g, return_info.label_size_before_canonical_repair, return_info.label_size_after_canonical_repair, 
		return_info.canonical_repair_remove_label_ratio, case_info.thread_num);
	return_info.time_repair = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin2).count() / 1e9; // s


	//cout << return_info.label_size_before_canonical_repair << endl;
	//cout << return_info.label_size_after_canonical_repair << endl;

	//cout << "print_L1:" << endl;
	//for (int i = 0; i < L_temp_595.size(); i++) {
	//	cout << "L[" << i << "]=";
	//	for (int j = 0; j < L_temp_595[i].size(); j++) {
	//		cout << "{" << L_temp_595[i][j].vertex << "," << L_temp_595[i][j].distance << "," << L_temp_595[i][j].parent_vertex << "}";
	//	}
	//	cout << "end" << endl;
	//}

	auto tempL = graph_hash_of_mixed_weighted_HL_PLL_v1_transform_labels_to_old_vertex_IDs(N, N, case_info.thread_num);

	//cout << "print_L2:" << endl;
	//for (int i = 0; i < L_temp_595.size(); i++) {
	//	cout << "L[" << i << "]=";
	//	for (int j = 0; j < L_temp_595[i].size(); j++) {
	//		cout << "{" << L_temp_595[i][j].vertex << "," << L_temp_595[i][j].distance << "," << L_temp_595[i][j].parent_vertex << "}";
	//	}
	//	cout << "end" << endl;
	//}

	for (int i = 0; i < N; i++) {
		if (!isIntree[i]) {
			L[i] = tempL[i];
			vector<two_hop_label_v1>().swap(tempL[i]);
		}
	}
	graph_hash_of_mixed_weighted_two_hop_clear_global_values();

	//case_info.two_hop_case_info.print_L();

	//cout << "end" << endl;

	return return_info;
}






