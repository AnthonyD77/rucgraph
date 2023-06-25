#pragma once

/*the following codes are for testing

---------------------------------------------------
a cpp file (try.cpp) for running the following test code:
----------------------------------------

#include <iostream>
#include <fstream>
using namespace std;

// header files in the Boost library: https://www.boost.org/
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

#include <graph_hash_of_mixed_weighted/HL/two_hop_v1/test_PLL_PSL_v1.h>


int main()
{
	test_PLL_PSL();
}

------------------------------------------------------------------------------------------
Commends for running the above cpp file on Linux:

g++ -std=c++17 -I/home/boost_1_75_0 -I/root/ysgraph try.cpp -lpthread -Ofast -o A
./A
rm A

(optional to put the above commends in run.sh, and then use the comment: sh run.sh)


*/
#include <graph_hash_of_mixed_weighted/HL/two_hop_v1/graph_hash_of_mixed_weighted_PLL_v1.h>
#include <graph_hash_of_mixed_weighted/HL/two_hop_v1/graph_hash_of_mixed_weighted_PSL_v1.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_generate_random_graph.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_read_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_save_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_shortest_paths.h>

void graph_hash_of_mixed_weighted_PLL_PSL_v1_check_correctness(graph_hash_of_mixed_weighted_two_hop_case_info_v1& case_info,
	graph_hash_of_mixed_weighted& instance_graph, int iteration_source_times, int iteration_terminal_times) {

	/*below is for checking whether the above labels are right (by randomly computing shortest paths)

	this function can only be used when 0 to n-1 is in the graph, i.e., the graph is an ideal graph

	*/

	boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(instance_graph.hash_of_vectors.size() - 1) };

	//graph_hash_of_mixed_weighted_print(instance_graph);

	for (int yy = 0; yy < iteration_source_times; yy++) {
		int source = dist(boost_random_time_seed);
		std::unordered_map<int, double> distances;
		std::unordered_map<int, int> predecessors;

		//source = 3; cout << "source = " << source << endl;

		graph_hash_of_mixed_weighted_shortest_paths_source_to_all(instance_graph, source, distances, predecessors);

		for (int xx = 0; xx < iteration_terminal_times; xx++) {

			int terminal = dist(boost_random_time_seed);

			//terminal = 4; cout << "terminal = " << terminal << endl;

			float dis = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance
			(case_info.L, case_info.reduction_measures_2019R2, case_info.reduction_measures_2019R1, case_info.f_2019R1, instance_graph, source, terminal);

			if (abs(dis - distances[terminal]) > 1e-4 && (dis < std::numeric_limits<float>::max() || distances[terminal] < std::numeric_limits<float>::max())) {
				cout << "source = " << source << endl;
				cout << "terminal = " << terminal << endl;
				cout << "source vector:" << endl;
				for (auto it = case_info.L[source].begin(); it != case_info.L[source].end(); it++) {
					cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
				}
				cout << endl;
				cout << "terminal vector:" << endl;
				for (auto it = case_info.L[terminal].begin(); it != case_info.L[terminal].end(); it++) {
					cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
				}
				cout << endl;

				cout << "dis = " << dis << endl;
				cout << "distances[terminal] = " << distances[terminal] << endl;
				cout << "abs(dis - distances[terminal]) > 1e-5!" << endl;
				getchar();
			}

			//cout << 0 << endl;
			//cout << source << " " << terminal << endl;
			//getchar();

			vector<pair<int, int>> path = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path(case_info.L,
				case_info.reduction_measures_2019R2, case_info.reduction_measures_2019R1, case_info.f_2019R1, instance_graph, source, terminal);

			float path_dis = 0;
			if (path.size() == 0) {
				if (source != terminal) { // disconnected
					path_dis = std::numeric_limits<float>::max();
				}
			}
			else {
				for (auto it = path.begin(); it != path.end(); it++) {
					path_dis = path_dis + graph_hash_of_mixed_weighted_edge_weight(instance_graph, it->first, it->second);
					if (path_dis > std::numeric_limits<float>::max()) {
						path_dis = std::numeric_limits<float>::max();
					}
				}
			}
			if (abs(dis - path_dis) > 1e-4 && (dis < std::numeric_limits<float>::max() || distances[terminal] < std::numeric_limits<float>::max())) {
				cout << "source = " << source << endl;
				cout << "terminal = " << terminal << endl;

				cout << "source vector:" << endl;
				for (auto it = case_info.L[source].begin(); it != case_info.L[source].end(); it++) {
					cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
				}
				cout << endl;
				cout << "terminal vector:" << endl;
				for (auto it = case_info.L[terminal].begin(); it != case_info.L[terminal].end(); it++) {
					cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
				}
				cout << endl;

				print_vector_pair_int(path);
				cout << "dis = " << dis << endl;
				cout << "path_dis = " << path_dis << endl;
				cout << "abs(dis - path_dis) > 1e-5!" << endl;
				getchar();
			}
		}

	}

}

void test_PLL_PSL() {

	/*parameters*/
	int iteration_graph_times = 1e2, iteration_source_times = 10, iteration_terminal_times = 10;
	int V = 100, E = 500, precision = 1, thread_num = 5;
	float ec_min = 0.1, ec_max = 1; // set ec_min=ec_max=1 for testing unweighted PLL_with_non_adj_reduction

	bool use_PLL = 1; // 1: PLL 0: PSL

	float avg_index_time = 0, avg_index_size_per_v = 0, avg_reduce_V_num_2019R1 = 0, avg_MG_num = 0;
	double avg_canonical_repair_remove_label_ratio = 0;


	bool weighted = true;
	if (ec_min == 1 && ec_max == 1) {
		weighted = false;
	}

	/*reduction method selection*/
	graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
	mm.use_2019R1 = 1;
	mm.use_2019R2 = 1;
	mm.use_enhanced2019R2 = 0;
	mm.use_non_adj_reduc_degree = 0;
	mm.max_degree_MG_enhanced2019R2 = 100;
	mm.max_labal_size = 6e9;
	mm.max_run_time_seconds = 1e9;
	mm.use_canonical_repair = true;

	/*iteration*/
	for (int i = 0; i < iteration_graph_times; i++) {
		cout << i << endl;

		/*input and output; below is for generating random new graph, or read saved graph*/
		int generate_new_graph = 1;
		std::unordered_set<int> generated_group_vertices;
		graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
		if (generate_new_graph == 1) {
			instance_graph = graph_hash_of_mixed_weighted_generate_random_graph(V, E, 0, 0, ec_min, ec_max, precision, boost_random_time_seed);
			graph_hash_of_mixed_weighted_save_graph_with_weight("simple_iterative_tests.txt", instance_graph, 0);
		}
		else {
			double lambda;
			graph_hash_of_mixed_weighted_read_graph_with_weight("simple_iterative_tests.txt", instance_graph, lambda);
		}
		//graph_hash_of_mixed_weighted_print(instance_graph);


		auto begin = std::chrono::high_resolution_clock::now();
		try {
			if (use_PLL) {
				graph_hash_of_mixed_weighted_PLL_v1(instance_graph, V + 1, weighted, thread_num, mm);
			}
			else {
				graph_hash_of_mixed_weighted_PSL_v1(instance_graph, V + 1, thread_num, mm);
			}
			if (0) {
				cout << "mm.time_initialization: " << mm.time_initialization << "s" << endl;
				cout << "mm.time_2019R1: " << mm.time_2019R1 << "s" << endl;
				cout << "mm.time_2019R2_or_enhanced_pre: " << mm.time_2019R2_or_enhanced_pre << "s" << endl;
				cout << "mm.time_2019R2_or_enhanced_fixlabels: " << mm.time_2019R2_or_enhanced_fixlabels << "s" << endl;
				cout << "mm.time_generate_labels: " << mm.time_generate_labels << "s" << endl;
				cout << "mm.time_canonical_repair1: " << mm.time_canonical_repair1 << "s" << endl;
				cout << "mm.time_canonical_repair2: " << mm.time_canonical_repair2 << "s" << endl;
				cout << "mm.time_update_old_IDs_in_labels: " << mm.time_update_old_IDs_in_labels << "s" << endl;
			}
		}
		catch (string s) {
			cout << s << endl;
			graph_hash_of_mixed_weighted_two_hop_clear_global_values();
			continue;
		}
		auto end = std::chrono::high_resolution_clock::now();
		double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
		avg_index_time = avg_index_time + runningtime / iteration_graph_times;

		avg_reduce_V_num_2019R1 = avg_reduce_V_num_2019R1 + (float)mm.reduce_V_num_2019R1 / iteration_graph_times;
		avg_MG_num = avg_MG_num + (float)mm.MG_num / iteration_graph_times;
		avg_canonical_repair_remove_label_ratio = avg_canonical_repair_remove_label_ratio + (double)mm.canonical_repair_remove_label_ratio / iteration_graph_times;

		/*debug*/
		if (0) {
			graph_hash_of_mixed_weighted_print(instance_graph);
			mm.print_L();
			mm.print_reduction_measures_2019R1();
			mm.print_reduction_measures_2019R2();
			mm.print_f_2019R1();
		}

		/*test canonical_repair proof*/
		if (1) {
			auto mm2 = mm;
			if (!use_PLL) { // use different method here
				graph_hash_of_mixed_weighted_PLL_v1(instance_graph, V + 1, weighted, 1, mm2); // single thread
			}
			else {
				graph_hash_of_mixed_weighted_PSL_v1(instance_graph, V + 1, 1, mm2); // single thread
			}

			auto& L1 = mm.L;
			auto& L2 = mm2.L;
			int size = L1.size();
			bool L1_is_l2 = true;
			for (int xx = 0; xx < size; xx++) {
				if (L1[xx].size() != L2[xx].size()) {
					if ((L1[xx].size() == 0 && L2[xx].size() == 1) || (L1[xx].size() == 1 && L2[xx].size() == 0)) {
					}
					else {
						L1_is_l2 == false;
						cout << "here" << endl;
						mm.print_L();
						mm2.print_L();
						getchar();
					}
				}
				else {
					int size2 = L1[xx].size();
					for (int yy = 0; yy < size2; yy++) {
						if (L1[xx][yy].vertex != L2[xx][yy].vertex) {
							L1_is_l2 == false;
							cout << "here" << endl;
							mm.print_L();
							mm2.print_L();
							getchar();
						}
					}
				}
			}
		}


		graph_hash_of_mixed_weighted_PLL_PSL_v1_check_correctness(mm, instance_graph, iteration_source_times, iteration_terminal_times);					  

		long long int index_size = 0;
		for (auto it = mm.L.begin(); it != mm.L.end(); it++) {
			index_size = index_size + (*it).size();
		}
		avg_index_size_per_v = avg_index_size_per_v + (float)index_size / V / iteration_graph_times;

		mm.clear_labels();
	}

	cout << "avg_index_time: " << avg_index_time << "s" << endl;
	cout << "avg_index_size_per_v: " << avg_index_size_per_v << endl;
	cout << "avg_reduce_V_num_2019R1: " << avg_reduce_V_num_2019R1 << endl;
	cout << "avg_MG_num: " << avg_MG_num << endl;
	cout << "avg_canonical_repair_remove_label_ratio: " << avg_canonical_repair_remove_label_ratio << endl;
}