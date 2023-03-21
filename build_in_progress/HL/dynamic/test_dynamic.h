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

#include <build_in_progress/HL/dynamic/test_dynamic.h>


int main()
{
	test_dynamic();
}

------------------------------------------------------------------------------------------
Commends for running the above cpp file on Linux:

g++ -std=c++17 -I/home/boost_1_75_0 -I/root/rucgraph try.cpp -lpthread -Ofast -o A
./A
rm A

(optional to put the above commends in run.sh, and then use the comment: sh run.sh)


*/
#include <build_in_progress/HL/dynamic/graph_hash_of_mixed_weighted_PLL_dynamic.h>
#include <build_in_progress/HL/dynamic/WeightIncreaseMaintenance.h>
#include <build_in_progress/HL/dynamic/WeightDecreaseMaintenance.h>
#include <build_in_progress/HL/sort_v/graph_hash_of_mixed_weighted_update_vertexIDs_by_degrees.h>
#include <graph_hash_of_mixed_weighted/two_graphs_operations/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID_2.h>
#include <graph_hash_of_mixed_weighted/random_graph/graph_hash_of_mixed_weighted_generate_random_graph.h>
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_read_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_save_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/common_algorithms/graph_hash_of_mixed_weighted_shortest_paths.h>





void graph_hash_of_mixed_weighted_PLL_PSL_v1_check_correctness_dynamic(graph_hash_of_mixed_weighted_two_hop_case_info_v1& case_info,
	graph_hash_of_mixed_weighted& instance_graph, int iteration_source_times, int iteration_terminal_times) {

	/*
	below is for checking whether the above labels are right (by randomly computing shortest paths)

	this function can only be used when 0 to n-1 is in the graph, i.e., the graph is an ideal graph
	*/

	boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(instance_graph.hash_of_vectors.size() - 1) };

	//graph_hash_of_mixed_weighted_print(instance_graph);

	for (int yy = 0; yy < iteration_source_times; yy++) {
		int source = dist(boost_random_time_seed);
		std::unordered_map<int, double> distances;
		std::unordered_map<int, int> predecessors;

		//source = 1; cout << "source = " << source << endl;

		graph_hash_of_mixed_weighted_shortest_paths_source_to_all(instance_graph, source, distances, predecessors);

		for (int xx = 0; xx < iteration_terminal_times; xx++) {

			int terminal = dist(boost_random_time_seed);

			//terminal = 6; cout << "terminal = " << terminal << endl;

			double dis = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance
			(case_info.L, case_info.reduction_measures_2019R2, case_info.reduction_measures_2019R1, case_info.f_2019R1, instance_graph, source, terminal);

			if (abs(dis - distances[terminal]) > 1e-2 && (dis < std::numeric_limits<weightTYPE>::max() || distances[terminal] < std::numeric_limits<double>::max())) {
				cout << "source = " << source << endl;
				cout << "terminal = " << terminal << endl;
				cout << "source vector:" << endl;
				for (auto it = case_info.L[source].begin(); it != case_info.L[source].end(); it++) {
					cout << "<" << it->vertex << "," << it->distance << ">";
				}
				cout << endl;
				cout << "terminal vector:" << endl;
				for (auto it = case_info.L[terminal].begin(); it != case_info.L[terminal].end(); it++) {
					cout << "<" << it->vertex << "," << it->distance << ">";
				}
				cout << endl;

				cout << "dis = " << dis << endl;
				cout << "distances[terminal] = " << distances[terminal] << endl;
				cout << "abs(dis - distances[terminal]) > 1e-5!" << endl;
				getchar();
			}
		}
	}
}

void graph_change_and_label_maintenance(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	int V, int weightIncrease_time, int weightDecrease_time, double weightChange_ratio) {

	while (weightIncrease_time + weightDecrease_time) {

		/*randomly select an edge*/
		pair<int, int> selected_edge;
		double selected_edge_weight;
		while (1) {
			boost::random::uniform_int_distribution<> dist_v1{ static_cast<int>(0), static_cast<int>(V - 1) };
			int v1 = dist_v1(boost_random_time_seed);
			if (instance_graph.degree(v1) > 0) {
				if (instance_graph.hash_of_vectors[v1].adj_vertices.size() > 0) {
					boost::random::uniform_int_distribution<> dist_v2{ static_cast<int>(0), static_cast<int>(instance_graph.hash_of_vectors[v1].adj_vertices.size() - 1) };
					int v2 = instance_graph.hash_of_vectors[v1].adj_vertices[dist_v2(boost_random_time_seed)].first;
					selected_edge = { v1,v2 };
					selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, v1, v2);
				}
				else {
					boost::random::uniform_int_distribution<> dist_v2{ static_cast<int>(0), static_cast<int>(instance_graph.hash_of_hashs[v1].size() - 1) };
					int r = dist_v2(boost_random_time_seed);
					auto it = instance_graph.hash_of_hashs[v1].begin();
					int id = 0;
					while (1) {
						if (id == r) {
							int v2 = it->first;
							selected_edge = { v1,v2 };
							selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, v1, v2);
							break;
						}
						else {
							it++, id++;
						}
					}
				}
				break;
			}
		}

		/*change weight*/
		if (weightIncrease_time >= weightDecrease_time) {
			weightIncrease_time--;

			/*debug*/
			//if (weightIncrease_time == 9) {
			//	selected_edge.first = 2, selected_edge.second = 4;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 8) {
			//	selected_edge.first = 6, selected_edge.second = 0;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 7) {
			//	selected_edge.first = 8, selected_edge.second = 4;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 6) {
			//	selected_edge.first = 0, selected_edge.second = 6;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 5) {
			//	selected_edge.first = 7, selected_edge.second = 0;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 4) {
			//	selected_edge.first = 1, selected_edge.second = 9;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 3) {
			//	selected_edge.first = 9, selected_edge.second = 1;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 2) {
			//	selected_edge.first = 0, selected_edge.second = 7;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 1) {
			//	selected_edge.first = 6, selected_edge.second = 3;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 0) {
			//	selected_edge.first = 0, selected_edge.second = 1;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}


			double new_ec = min(selected_edge_weight * (1 + weightChange_ratio), 1e6);
			graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight

			/*maintain labels*/
			WeightIncreaseMaintenance(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec);



			//cout << "1ec change " << selected_edge.first << " " << selected_edge.second << " " << selected_edge_weight * (1 + weightChange_ratio) << endl;
			//mm.print_L();
			//mm.print_PPR();
		}
		else {
			weightDecrease_time--;

			/*debug*/
			//if (weightDecrease_time == 9) {
			//	selected_edge.first = 1, selected_edge.second = 9;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightDecrease_time == 8) {
			//	selected_edge.first = 2, selected_edge.second = 4;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightDecrease_time == 7) {
			//	selected_edge.first = 4, selected_edge.second = 2;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightDecrease_time == 6) {
			//	selected_edge.first = 2, selected_edge.second = 4;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightDecrease_time == 5) {
			//	selected_edge.first = 0, selected_edge.second = 1;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightDecrease_time == 4) {
			//	selected_edge.first = 5, selected_edge.second = 2;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightDecrease_time == 3) {
			//	selected_edge.first = 9, selected_edge.second = 1;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightDecrease_time == 2) {
			//	selected_edge.first = 5, selected_edge.second = 1;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightDecrease_time == 1) {
			//	selected_edge.first = 7, selected_edge.second = 4;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightDecrease_time == 0) {
			//	selected_edge.first = 2, selected_edge.second = 5;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}

			double new_ec = max( selected_edge_weight * (1 - weightChange_ratio), 1e-2);
			graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight

			/*maintain labels*/
			WeightDecreaseMaintenance(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec);

			//cout << "2ec change " << selected_edge.first << " " << selected_edge.second << " " << selected_edge_weight * (1 - weightChange_ratio) << endl;
			//mm.print_L();
			//mm.print_PPR();
		}
	}

}

void test_dynamic() {

	/*parameters*/
	int iteration_graph_times = 1e5, iteration_source_times = 10, iteration_terminal_times = 10;
	int V = 100, E = 150, precision = 1, thread_num = 1;
	double ec_min = 1, ec_max = 10; // set ec_min=ec_max=1 for testing unweighted PLL_with_non_adj_reduction

	int weightIncrease_time = 100, weightDecrease_time = 100;
	double weightChange_ratio = 0.2;

	double avg_index_time = 0, avg_index_size_per_v = 0, avg_reduce_V_num_2019R1 = 0, avg_MG_num = 0;
	double avg_canonical_repair_remove_label_ratio = 0;

	/*reduction method selection*/
	graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
	mm.use_2019R1 = 0;
	mm.use_2019R2 = 0;
	mm.use_enhanced2019R2 = 0;
	mm.use_non_adj_reduc_degree = 0;
	mm.max_degree_MG_enhanced2019R2 = 100;
	mm.max_labal_size = 6e9;
	mm.max_run_time_seconds = 1e9;
	mm.use_canonical_repair = false; // canonical_repair needs to modify for PPR

	/*iteration*/
	for (int i = 0; i < iteration_graph_times; i++) {
		cout << i << endl;

		/*input and output; below is for generating random new graph, or read saved graph*/
		int generate_new_graph = 1;
		std::unordered_set<int> generated_group_vertices;
		graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
		if (generate_new_graph == 1) {
			instance_graph = graph_hash_of_mixed_weighted_generate_random_graph(V, E, 0, 0, ec_min, ec_max, precision, boost_random_time_seed);
			instance_graph = graph_hash_of_mixed_weighted_update_vertexIDs_by_degrees_large_to_small(instance_graph); // sort vertices
			graph_hash_of_mixed_weighted_save_graph_with_weight("simple_iterative_tests.txt", instance_graph, 0);
		}
		else {
			double lambda;
			graph_hash_of_mixed_weighted_read_graph_with_weight("simple_iterative_tests.txt", instance_graph, lambda);
		}
		//graph_hash_of_mixed_weighted_print(instance_graph);


		auto begin = std::chrono::high_resolution_clock::now();
		try {
			graph_hash_of_mixed_weighted_PLL_dynamic(instance_graph, V + 1, thread_num, mm);
			if (0) {
				cout << "mm.time_initialization: " << mm.time_initialization << "s" << endl;
				cout << "mm.time_2019R1: " << mm.time_2019R1 << "s" << endl;
				cout << "mm.time_2019R2_or_enhanced_pre: " << mm.time_2019R2_or_enhanced_pre << "s" << endl;
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
		avg_reduce_V_num_2019R1 = avg_reduce_V_num_2019R1 + (double)mm.reduce_V_num_2019R1 / iteration_graph_times;
		avg_MG_num = avg_MG_num + (double)mm.MG_num / iteration_graph_times;
		avg_canonical_repair_remove_label_ratio = avg_canonical_repair_remove_label_ratio + (double)mm.canonical_repair_remove_label_ratio / iteration_graph_times;


		//mm.print_L();
		//mm.print_PPR();
		/*debug*/
		if (0) {
			graph_hash_of_mixed_weighted_print(instance_graph);
			mm.print_L();
			mm.print_reduction_measures_2019R1();
			mm.print_reduction_measures_2019R2();
			mm.print_f_2019R1();
		}

		/*dynamic maintenance*/
		graph_change_and_label_maintenance(instance_graph, mm, V, weightIncrease_time, weightDecrease_time, weightChange_ratio);

		graph_hash_of_mixed_weighted_PLL_PSL_v1_check_correctness_dynamic(mm, instance_graph, iteration_source_times, iteration_terminal_times);

		long long int index_size = 0;
		for (auto it = mm.L.begin(); it != mm.L.end(); it++) {
			index_size = index_size + (*it).size();
		}
		avg_index_size_per_v = avg_index_size_per_v + (double)index_size / V / iteration_graph_times;

		mm.clear_labels();
	}

	cout << "avg_index_time: " << avg_index_time << "s" << endl;
	cout << "avg_index_size_per_v: " << avg_index_size_per_v << endl;
	cout << "avg_reduce_V_num_2019R1: " << avg_reduce_V_num_2019R1 << endl;
	cout << "avg_MG_num: " << avg_MG_num << endl;
	cout << "avg_canonical_repair_remove_label_ratio: " << avg_canonical_repair_remove_label_ratio << endl;
}


