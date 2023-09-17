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

g++ -std=c++17 -I/home/boost_1_75_0 -I/root/rucgraph try.cpp -lpthread -O3 -o A
./A
rm A

(optional to put the above commends in run.sh, and then use the comment: sh run.sh)


*/
#include <build_in_progress/HL/dynamic/PLL_dynamic.h>
#include <build_in_progress/HL/dynamic/WeightIncreaseMaintenance_improv.h>
#include <build_in_progress/HL/dynamic/WeightDecreaseMaintenance_improv.h>
#include <build_in_progress/HL/dynamic/WeightIncrease2021.h>
#include <build_in_progress/HL/dynamic/WeightDecrease2021.h>
#include <build_in_progress/HL/dynamic/WeightDecrease2014.h>
#include <build_in_progress/HL/dynamic/WeightIncrease2019_multiThread.h>
#include <build_in_progress/HL/sort_v/graph_hash_of_mixed_weighted_update_vertexIDs_by_degrees.h>
#include <graph_hash_of_mixed_weighted/two_graphs_operations/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID_2.h>
#include <graph_hash_of_mixed_weighted/random_graph/graph_hash_of_mixed_weighted_generate_random_graph.h>
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_read_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_save_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/common_algorithms/graph_hash_of_mixed_weighted_shortest_paths.h>
#include <build_in_progress/HL/dynamic/clean_labels.h>

/*check_correctness*/

void check_correctness_dynamic(graph_hash_of_mixed_weighted_two_hop_case_info_v1& case_info,
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

		//source = 6; cout << "source = " << source << endl;

		graph_hash_of_mixed_weighted_shortest_paths_source_to_all(instance_graph, source, distances, predecessors);

		for (int xx = 0; xx < iteration_terminal_times; xx++) {

			int terminal = dist(boost_random_time_seed);

			//terminal = 0; cout << "terminal = " << terminal << endl;

			double dis = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(case_info.L, source, terminal);

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
	int V, int weightIncrease_time, int weightDecrease_time, double weightChange_ratio, int thread_num, double& avg_maintain_time) {

	ThreadPool pool_dynamic(thread_num);
	std::vector<std::future<int>> results_dynamic;

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
			//	selected_edge.first = 1, selected_edge.second = 8;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 8) {
			//	selected_edge.first = 0, selected_edge.second = 6;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 7) {
			//	selected_edge.first = 3, selected_edge.second = 0;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 6) {
			//	selected_edge.first = 2, selected_edge.second = 1;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 5) {
			//	selected_edge.first = 2, selected_edge.second = 1;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 4) {
			//	selected_edge.first = 0, selected_edge.second = 3;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 3) {
			//	selected_edge.first = 5, selected_edge.second = 0;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 2) {
			//	selected_edge.first = 7, selected_edge.second = 5;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 1) {
			//	selected_edge.first = 0, selected_edge.second = 2;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightIncrease_time == 0) {
			//	selected_edge.first = 2, selected_edge.second = 0;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}

			double new_ec = selected_edge_weight * (1 + weightChange_ratio);
			if (new_ec > 1e6) {
				weightIncrease_time++;
				continue;
			}
			graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight

			auto begin = std::chrono::high_resolution_clock::now();

			/*maintain labels*/
			//WeightIncrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight);
			//WeightIncreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight);
			WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic, 1e-1);
			//WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, 1e2);
			//cout << "1ec change " << selected_edge.first << " " << selected_edge.second << " " << selected_edge_weight * (1 + weightChange_ratio) << endl;
			//mm.print_L();
			//mm.print_PPR();

			auto end = std::chrono::high_resolution_clock::now();
			avg_maintain_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
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
			//	selected_edge.first = 1, selected_edge.second = 2;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightDecrease_time == 3) {
			//	selected_edge.first = 3, selected_edge.second = 5;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightDecrease_time == 2) {
			//	selected_edge.first = 0, selected_edge.second = 1;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightDecrease_time == 1) {
			//	selected_edge.first = 0, selected_edge.second = 3;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}
			//if (weightDecrease_time == 0) {
			//	selected_edge.first = 1, selected_edge.second = 2;
			//	selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
			//}

			double new_ec = selected_edge_weight * (1 - weightChange_ratio);
			if (new_ec < 1e-2) {
				weightDecrease_time++;
				continue;
			}
			graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight

			auto begin = std::chrono::high_resolution_clock::now();

			/*maintain labels*/
			//WeightDecreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, new_ec);
			//WeightDecrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, new_ec);
			WeightDecrease2014(instance_graph, mm, selected_edge.first, selected_edge.second, new_ec);
			//cout << "2ec change " << selected_edge.first << " " << selected_edge.second << " " << selected_edge_weight * (1 - weightChange_ratio) << endl;
			//mm.print_L();
			//mm.print_PPR();

			auto end = std::chrono::high_resolution_clock::now();
			avg_maintain_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
		}
	}
}

void test_dynamic() {

	/*parameters*/
	int iteration_graph_times = 1e2, iteration_source_times = 10, iteration_terminal_times = 10;
	int V = 100, E = 500, precision = 1, thread_num = 10;
	double ec_min = 1, ec_max = 10;

	int weightIncrease_time = 30, weightDecrease_time = 0;
	double weightChange_ratio = 0.2;

	double avg_index_time = 0, avg_index_size_per_v = 0, avg_maintain_time = 0;

	/*iteration*/
	for (int i = 0; i < iteration_graph_times; i++) {
		cout << "iteration " << i << endl;
		//getchar();
		/*reduction method selection*/
		graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
		mm.max_labal_size = 6e8;
		mm.max_run_time_seconds = 1e2;

		/*input and output; below is for generating random new graph, or read saved graph*/
		int generate_new_graph = 1; // this value is for debug
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

		auto begin = std::chrono::high_resolution_clock::now();
		try {
			PLL_dynamic(instance_graph, V + 1, thread_num, mm);
			if (0) {
				graph_hash_of_mixed_weighted_print(instance_graph);
				mm.print_L();
				mm.print_PPR();
			}
			//mm.print_PPR();
			//binary_save_PPR("PPR.bin", mm.PPR);
			//binary_read_PPR("PPR.bin", mm.PPR);
			//PPR_type x(V + 1);
			//PPR_type(x).swap(mm.PPR);
			//mm.print_PPR();
		}
		catch (string s) {
			cout << s << endl;
			two_hop_clear_global_values();
			continue;
		}
		auto end = std::chrono::high_resolution_clock::now();
		double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
		avg_index_time = avg_index_time + runningtime / iteration_graph_times;

		/*dynamic maintenance*/
		ThreadPool pool_dynamic(thread_num);
		std::vector<std::future<int>> results_dynamic;
		clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic);
		initialize_global_values_dynamic(V, thread_num);
		graph_change_and_label_maintenance(instance_graph, mm, V, weightIncrease_time, weightDecrease_time, weightChange_ratio, thread_num, avg_maintain_time);
		check_correctness_dynamic(mm, instance_graph, iteration_source_times, iteration_terminal_times);

		long long int index_size = 0;
		for (auto it = mm.L.begin(); it != mm.L.end(); it++) {
			index_size = index_size + (*it).size();
		}
		avg_index_size_per_v = avg_index_size_per_v + (double)index_size / V / iteration_graph_times;
	}

	cout << "avg_index_time: " << avg_index_time << "s" << endl;
	cout << "avg_index_size_per_v: " << avg_index_size_per_v << endl;
	cout << "avg_maintain_time: " << (double)avg_maintain_time / iteration_graph_times << "s" << endl;
}









/*
compare_speed

the speed of WeightDecreaseMaintenance_improv increases with thread_num when thread_num<10, and barely changes when thread_num>10
*/

void compare_speed() {

	/*parameters*/
	int iteration_graph_times = 5e0, weightChange_time = 50; 
	double weightChange_ratio = 0.5;
	vector<int> thread_nums = { 1, 5, 10, 20 };

	int V = 5000, E = 20000, precision = 1;
	double ec_min = 1, ec_max = 10;
	
	bool use_WeightIncreaseMaintenance = 1, use_WeightIncreaseMaintenance_improv = 1, 
		use_WeightDecrease2014 = 1, use_WeightDecreaseMaintenance = 1, use_WeightDecreaseMaintenance_improv = 1;
	
	/*iteration*/
	cout << "V = " << V << "	E =" << E << endl;
	
	for (int i = 0; i < iteration_graph_times; i++) {
		cout << "iteration_graph_times = " << i << endl;

		/*reduction method selection*/
		graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;

		graph_hash_of_mixed_weighted instance_graph = graph_hash_of_mixed_weighted_generate_random_graph(V, E, 0, 0, ec_min, ec_max, precision, boost_random_time_seed);
		instance_graph = graph_hash_of_mixed_weighted_update_vertexIDs_by_degrees_large_to_small(instance_graph); // sort vertices
		PLL_dynamic(instance_graph, V + 1, 20, mm);

		for (int t_num : thread_nums) {
			cout << "thread_num = " << t_num << endl;

			double avg_time_WeightIncreaseMaintenance = 0, avg_time_WeightIncreaseMaintenance_improv = 0,
				avg_time_WeightDecrease2014 = 0, avg_time_WeightDecreaseMaintenance = 0, avg_time_WeightDecreaseMaintenance_improv = 0;

			initialize_global_values_dynamic(V, t_num);
			ThreadPool pool_dynamic(t_num);
			std::vector<std::future<int>> results_dynamic;

			/*selected_edge*/
			vector<pair<int, int>> selected_edges;
			boost::random::uniform_int_distribution<> dist_v1{ static_cast<int>(0), static_cast<int>(V - 1) };
			while (selected_edges.size() < weightChange_time) {
				int v1 = dist_v1(boost_random_time_seed);
				if (instance_graph.degree(v1) > 0) {
					if (instance_graph.hash_of_vectors[v1].adj_vertices.size() > 0) {
						boost::random::uniform_int_distribution<> dist_v2{ static_cast<int>(0), static_cast<int>(instance_graph.hash_of_vectors[v1].adj_vertices.size() - 1) };
						int v2 = instance_graph.hash_of_vectors[v1].adj_vertices[dist_v2(boost_random_time_seed)].first;
						selected_edges.push_back({ v1, v2 });
					}
					else {
						boost::random::uniform_int_distribution<> dist_v2{ static_cast<int>(0), static_cast<int>(instance_graph.hash_of_hashs[v1].size() - 1) };
						int r = dist_v2(boost_random_time_seed);
						auto it = instance_graph.hash_of_hashs[v1].begin();
						int id = 0;
						while (1) {
							if (id == r) {
								int v2 = it->first;
								selected_edges.push_back({ v1, v2 });
								break;
							}
							else {
								it++, id++;
							}
						}
					}
				}
			}

			if (use_WeightIncreaseMaintenance) {
				auto g = instance_graph;
				auto mm2 = mm;
				for (int j = 0; j < weightChange_time; j++) {
					double old_ec = graph_hash_of_mixed_weighted_edge_weight(g, selected_edges[j].first, selected_edges[j].second);
					double new_ec = min(old_ec * (1 + weightChange_ratio), 1e6);
					graph_hash_of_mixed_weighted_add_edge(g, selected_edges[j].first, selected_edges[j].second, new_ec);
					auto begin = std::chrono::high_resolution_clock::now();
					WeightIncrease2021(g, mm2, selected_edges[j].first, selected_edges[j].second, old_ec);
					auto end = std::chrono::high_resolution_clock::now();
					double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
					avg_time_WeightIncreaseMaintenance += runningtime / iteration_graph_times / weightChange_time;
				}
			}

			if (use_WeightIncreaseMaintenance_improv) {
				auto g = instance_graph;
				auto mm2 = mm;
				for (int j = 0; j < weightChange_time; j++) {
					double old_ec = graph_hash_of_mixed_weighted_edge_weight(g, selected_edges[j].first, selected_edges[j].second);
					double new_ec = min(old_ec * (1 + weightChange_ratio), 1e6);
					graph_hash_of_mixed_weighted_add_edge(g, selected_edges[j].first, selected_edges[j].second, new_ec);
					auto begin = std::chrono::high_resolution_clock::now();
					WeightIncreaseMaintenance_improv(g, mm2, selected_edges[j].first, selected_edges[j].second, old_ec);
					auto end = std::chrono::high_resolution_clock::now();
					double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
					avg_time_WeightIncreaseMaintenance_improv += runningtime / iteration_graph_times / weightChange_time;
				}
			}

			if (use_WeightDecrease2014) {
				auto g = instance_graph;
				auto mm2 = mm;
				for (int j = 0; j < weightChange_time; j++) {
					double new_ec = max(graph_hash_of_mixed_weighted_edge_weight(g, selected_edges[j].first, selected_edges[j].second) * (1 - weightChange_ratio), 1e-2);
					graph_hash_of_mixed_weighted_add_edge(g, selected_edges[j].first, selected_edges[j].second, new_ec);
					auto begin = std::chrono::high_resolution_clock::now();
					WeightDecrease2014(g, mm2, selected_edges[j].first, selected_edges[j].second, new_ec);
					auto end = std::chrono::high_resolution_clock::now();
					double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
					avg_time_WeightDecrease2014 += runningtime / iteration_graph_times / weightChange_time;
					//cout << "WeightDecreaseMaintenance " << runningtime << "s" << endl;
				}
			}

			if (use_WeightDecreaseMaintenance) {
				auto g = instance_graph;
				auto mm2 = mm;
				for (int j = 0; j < weightChange_time; j++) {
					double new_ec = max(graph_hash_of_mixed_weighted_edge_weight(g, selected_edges[j].first, selected_edges[j].second) * (1 - weightChange_ratio), 1e-2);
					graph_hash_of_mixed_weighted_add_edge(g, selected_edges[j].first, selected_edges[j].second, new_ec);
					auto begin = std::chrono::high_resolution_clock::now();
					WeightDecrease2021(g, mm2, selected_edges[j].first, selected_edges[j].second, new_ec);
					auto end = std::chrono::high_resolution_clock::now();
					double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
					avg_time_WeightDecreaseMaintenance += runningtime / iteration_graph_times / weightChange_time;
					//cout << "WeightDecreaseMaintenance " << runningtime << "s" << endl;
				}
			}

			if (use_WeightDecreaseMaintenance_improv) {
				auto g = instance_graph;
				auto mm2 = mm;
				for (int j = 0; j < weightChange_time; j++) {
					double new_ec = max(graph_hash_of_mixed_weighted_edge_weight(g, selected_edges[j].first, selected_edges[j].second) * (1 - weightChange_ratio), 1e-2);
					graph_hash_of_mixed_weighted_add_edge(g, selected_edges[j].first, selected_edges[j].second, new_ec);
					auto begin = std::chrono::high_resolution_clock::now();
					WeightDecreaseMaintenance_improv(g, mm2, selected_edges[j].first, selected_edges[j].second, new_ec);
					auto end = std::chrono::high_resolution_clock::now();
					double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
					avg_time_WeightDecreaseMaintenance_improv += runningtime / iteration_graph_times / weightChange_time;
					//cout << "WeightDecreaseMaintenance_improv " << runningtime << "s" << endl;
				}
			}



			cout << "avg_time_WeightIncreaseMaintenance = " << avg_time_WeightIncreaseMaintenance << "s" << endl
				<< "avg_time_WeightIncreaseMaintenance_improv = " << avg_time_WeightIncreaseMaintenance_improv << "s" << endl
				<< "avg_time_WeightDecrease2014 = " << avg_time_WeightDecrease2014 << "s" << endl
				<< "avg_time_WeightDecreaseMaintenance = " << avg_time_WeightDecreaseMaintenance << "s" << endl
				<< "avg_time_WeightDecreaseMaintenance_improv = " << avg_time_WeightDecreaseMaintenance_improv << "s" << endl;
			cout << "avg_time_WeightIncreaseMaintenance / avg_time_WeightIncreaseMaintenance_improv = " << avg_time_WeightIncreaseMaintenance / avg_time_WeightIncreaseMaintenance_improv << endl
				<< "avg_time_WeightDecreaseMaintenance / avg_time_WeightDecreaseMaintenance_improv = " << avg_time_WeightDecreaseMaintenance / avg_time_WeightDecreaseMaintenance_improv << endl << endl;
		}
	}
}








/*
example
*/
void PLL_PPR_example() {

	graph_hash_of_mixed_weighted instance_graph;
	graph_hash_of_mixed_weighted_add_edge(instance_graph, 0, 1, 4);
	graph_hash_of_mixed_weighted_add_edge(instance_graph, 0, 3, 5);
	graph_hash_of_mixed_weighted_add_edge(instance_graph, 0, 4, 2);
	graph_hash_of_mixed_weighted_add_edge(instance_graph, 0, 5, 9);
	graph_hash_of_mixed_weighted_add_edge(instance_graph, 1, 2, 6);
	graph_hash_of_mixed_weighted_add_edge(instance_graph, 1, 4, 7);
	graph_hash_of_mixed_weighted_add_edge(instance_graph, 1, 6, 5);
	graph_hash_of_mixed_weighted_add_edge(instance_graph, 2, 3, 1);
	graph_hash_of_mixed_weighted_add_edge(instance_graph, 2, 6, 4);
	graph_hash_of_mixed_weighted_add_edge(instance_graph, 3, 5, 15);
	graph_hash_of_mixed_weighted_add_edge(instance_graph, 4, 5, 3);

	graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
	PLL_dynamic(instance_graph, 7, 1, mm);

	mm.print_L();
	mm.print_PPR();
}









