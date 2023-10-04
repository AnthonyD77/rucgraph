#include <iostream>
#include <fstream>
using namespace std;

// header files in the Boost library: https://www.boost.org/
#include <boost/random.hpp>
#include <boost/range/algorithm.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

#include <build_in_progress/HL/dynamic/PLL_dynamic.h>
#include <build_in_progress/HL/dynamic/WeightIncreaseMaintenance_improv_multiThread.h>
#include <build_in_progress/HL/dynamic/WeightDecreaseMaintenance_improv_multiThread.h>
#include <build_in_progress/HL/dynamic/WeightIncrease2021_multiThread.h>
#include <build_in_progress/HL/dynamic/WeightDecrease2021_multiThread.h>
#include <build_in_progress/HL/dynamic/WeightDecrease2014_multiThread.h>
#include <build_in_progress/HL/dynamic/WeightIncrease2019_multiThread.h>
#include <build_in_progress/HL/sort_v/graph_hash_of_mixed_weighted_update_vertexIDs_by_degrees.h>
#include <graph_hash_of_mixed_weighted/two_graphs_operations/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID_2.h>
#include <graph_hash_of_mixed_weighted/random_graph/graph_hash_of_mixed_weighted_generate_random_graph.h>
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_read_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_save_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_binary_save_read.h>
#include <graph_hash_of_mixed_weighted/common_algorithms/graph_hash_of_mixed_weighted_shortest_paths.h>
#include <build_in_progress/HL/dynamic/clean_labels.h>
#include <text_mining/binary_save_read_vector_of_vectors.h>
#include <graph_v_of_v_idealID/graph_v_of_v_idealID_to_graph_hash_of_mixed_weighted.h>


void generate_L_PPR() {

	vector<string> data_names = { "astro", "condmat", "github", "google", "youtube", "skitter" };
	string path = "dynamicHL//";
	int thread_num = 50;
	graph_hash_of_mixed_weighted g;
	graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
	ThreadPool pool_dynamic(thread_num);
	std::vector<std::future<int>> results_dynamic;

	string save_file_name;
	ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);

	for (auto s : data_names) {
		g = graph_hash_of_mixed_weighted_binary_read(path + s + "_random.bin");
		PLL_dynamic(g, g.hash_of_vectors.size(), thread_num, mm);
		clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic, 80);
		binary_save_PPR(path + s + "_PPR_random.bin", mm.PPR);
		binary_save_vector_of_vectors(path + s + "_L_random.bin", mm.L);
		outputFile.open(path + s + "_L_random_generation.txt");
		outputFile << "time_generate_labels=" << mm.time_generate_labels << "s" << endl;
		outputFile.close();

		g = graph_hash_of_mixed_weighted_binary_read(path + s + "_unique.bin");
		PLL_dynamic(g, g.hash_of_vectors.size(), thread_num, mm);
		clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic, 80);
		binary_save_PPR(path + s + "_PPR_unique.bin", mm.PPR);
		binary_save_vector_of_vectors(path + s + "_L_unique.bin", mm.L);
		outputFile.open(path + s + "_L_unique_generation.txt");
		outputFile << "time_generate_labels=" << mm.time_generate_labels << "s" << endl;
		outputFile.close();
	}
}

void exp_element(string data_name, double weightChange_ratio, int change_times, double max_Maintain_time, int thread_num) {

	ThreadPool pool_dynamic(thread_num);
	std::vector<std::future<int>> results_dynamic;

	ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);

	string path = "dynamicHL//";
	graph_v_of_v_idealID instance_graph;
	graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
	vector<pair<int, int>> selected_edges;

	for (int i = 0; i < 2; i++) {

		string weight_type;
		if (i == 0) {
			weight_type = "unique";
		}
		else {
			weight_type = "random";
		}
		graph_hash_of_mixed_weighted instance_graph_initial_hash = graph_hash_of_mixed_weighted_binary_read(path + data_name + "_" + weight_type + ".bin");
		graph_v_of_v_idealID instance_graph_initial = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID_2(instance_graph_initial_hash, instance_graph_initial_hash.hash_of_vectors.size());
		graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
		binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm_initial.PPR);
		binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm_initial.L);
		outputFile.open("exp_" + data_name + "_T_" + to_string(thread_num) + "_changeRatio_" + to_string((int)(weightChange_ratio * 100)) + "_" + weight_type + ".csv");

		outputFile << "2014DE_time,2019IN_time,2021DE_time,2021DE_query_times,2021IN_time,2021IN_query_times,newDE_time,newDE_query_times,newIN_time,newIN_query_times,2014+2019_time,newDE2021IN_time,newDEIN_time," <<
			"L_bit_size_initial(1),PPR_bit_size_initial,L_bit_size_afterM1,PPR_bit_size_afterM1,L_bit_size_afterM2,PPR_bit_size_afterM2,L_bit_size_afterClean,PPR_bit_size_afterClean,cleanL_time,cleanPPR_time,rege_time" << endl;

		vector<double> _2014DE_time(change_times, 0), _2019IN_time(change_times, 0), _2021DE_time(change_times, 0), _2021DE_query_times(change_times, 0), _2021IN_time(change_times, 0), _2021IN_query_times(change_times, 0),
			_newDE_time(change_times, 0), _newDE_query_times(change_times, 0), _newIN_time(change_times, 0), _newIN_query_times(change_times, 0), _20142019_time(change_times, 0),
			_newDE2021IN_time(change_times, 0), _newDEIN_time(change_times, 0);
		double L_bit_size_initial = 0, PPR_bit_size_initial = 0, L_bit_size_afterM1 = 0, PPR_bit_size_afterM1 = 0, L_bit_size_afterM2 = 0, PPR_bit_size_afterM2 = 0,
			L_bit_size_afterClean = 0, PPR_bit_size_afterClean = 0, cleanL_time = 0, cleanPPR_time = 0, rege_time = 0;

		/*IN*/
		if (1) {
			instance_graph = instance_graph_initial;
			int V = instance_graph.size();
			vector<pair<int, int>>().swap(selected_edges);
			int left_change_times = change_times;
			while (left_change_times) {
				/*randomly select an edge*/
				pair<int, int> selected_edge;
				double selected_edge_weight;
				while (1) {
					boost::random::uniform_int_distribution<> dist_v1{ static_cast<int>(0), static_cast<int>(V - 1) };
					int v1 = dist_v1(boost_random_time_seed);
					if (instance_graph[v1].size() > 0) {
						boost::random::uniform_int_distribution<> dist_v2{ static_cast<int>(0), static_cast<int>(instance_graph[v1].size() - 1) };
						int rand = dist_v2(boost_random_time_seed);
						int v2 = instance_graph[v1][rand].first;
						selected_edge = { v1,v2 };
						selected_edge_weight = instance_graph[v1][rand].second;
						break;
					}
				}
				double new_ec = selected_edge_weight * (1 + weightChange_ratio);
				if (new_ec > 1e6) {
					continue;
				}
				else {
					left_change_times--;
					graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
					selected_edges.push_back(selected_edge);
				}
			}

			/*2019*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.size();
				initialize_global_values_dynamic(V, thread_num);

				if (//data_name == "google" || data_name == "youtube" || data_name == "skitter"
					0) {
					for (int i = 0; i < change_times; i++) {
						_2019IN_time[i] = INT_MAX;
					}
				}
				else {
					for (int k = 0; k < change_times; k++) {
						pair<int, int> selected_edge = selected_edges[k];
						double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
						double new_ec = selected_edge_weight * (1 + weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight						
						try {
							auto begin = std::chrono::high_resolution_clock::now();
							WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic, max_Maintain_time);
							//WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, max_Maintain_time);
							_2019IN_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
						}
						catch (string s) {
							_2019IN_time[k] = INT_MAX;
						}
					}
				}
			}

			/*2021*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.size();
				initialize_global_values_dynamic(V, thread_num);

				for (int k = 0; k < change_times; k++) {
					global_query_times = 0;
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					double new_ec = selected_edge_weight * (1 + weightChange_ratio);
					graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight						
					auto begin = std::chrono::high_resolution_clock::now();
					WeightIncrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic);
					_2021IN_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					_2021IN_query_times[k] = global_query_times;
				}
			}

			/*new*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.size();
				initialize_global_values_dynamic(V, thread_num);

				for (int k = 0; k < change_times; k++) {
					global_query_times = 0;
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					double new_ec = selected_edge_weight * (1 + weightChange_ratio);
					graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight						
					auto begin = std::chrono::high_resolution_clock::now();
					WeightIncreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic);
					_newIN_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					_newIN_query_times[k] = global_query_times;
				}
			}
		}

		/*DE*/
		if (1) {
			instance_graph = instance_graph_initial;
			int V = instance_graph.size();
			vector<pair<int, int>>().swap(selected_edges);
			int left_change_times = change_times;
			while (left_change_times) {
				/*randomly select an edge*/
				pair<int, int> selected_edge;
				double selected_edge_weight;
				while (1) {
					boost::random::uniform_int_distribution<> dist_v1{ static_cast<int>(0), static_cast<int>(V - 1) };
					int v1 = dist_v1(boost_random_time_seed);
					if (instance_graph[v1].size() > 0) {
						boost::random::uniform_int_distribution<> dist_v2{ static_cast<int>(0), static_cast<int>(instance_graph[v1].size() - 1) };
						int rand = dist_v2(boost_random_time_seed);
						int v2 = instance_graph[v1][rand].first;
						selected_edge = { v1,v2 };
						selected_edge_weight = instance_graph[v1][rand].second;
						break;
					}
				}
				double new_ec = selected_edge_weight * (1 - weightChange_ratio);
				if (new_ec < 1e-2) {
					continue;
				}
				else {
					left_change_times--;
					graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
					selected_edges.push_back(selected_edge);
				}
			}

			/*2014*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.size();
				initialize_global_values_dynamic(V, thread_num);

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					double new_ec = selected_edge_weight * (1 - weightChange_ratio);
					graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
					auto begin = std::chrono::high_resolution_clock::now();
					WeightDecrease2014(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic);
					_2014DE_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
				}
			}

			/*2021*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.size();
				initialize_global_values_dynamic(V, thread_num);

				for (int k = 0; k < change_times; k++) {
					global_query_times = 0;
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					double new_ec = selected_edge_weight * (1 - weightChange_ratio);
					graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight					
					auto begin = std::chrono::high_resolution_clock::now();
					WeightDecrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic);
					_2021DE_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					_2021DE_query_times[k] = global_query_times;
				}
			}

			/*new*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.size();
				initialize_global_values_dynamic(V, thread_num);

				for (int k = 0; k < change_times; k++) {
					global_query_times = 0;
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					double new_ec = selected_edge_weight * (1 - weightChange_ratio);
					graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight						
					auto begin = std::chrono::high_resolution_clock::now();
					WeightDecreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic);
					_newDE_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					_newDE_query_times[k] = global_query_times;
				}
			}
		}

		/*mixed*/
		if (1) {
			instance_graph = instance_graph_initial;
			int V = instance_graph.size();
			vector<pair<int, int>>().swap(selected_edges);
			int left_change_times = change_times;
			while (left_change_times) {
				/*randomly select an edge*/
				pair<int, int> selected_edge;
				double selected_edge_weight;
				while (1) {
					boost::random::uniform_int_distribution<> dist_v1{ static_cast<int>(0), static_cast<int>(V - 1) };
					int v1 = dist_v1(boost_random_time_seed);
					if (instance_graph[v1].size() > 0) {
						boost::random::uniform_int_distribution<> dist_v2{ static_cast<int>(0), static_cast<int>(instance_graph[v1].size() - 1) };
						int rand = dist_v2(boost_random_time_seed);
						int v2 = instance_graph[v1][rand].first;
						selected_edge = { v1,v2 };
						selected_edge_weight = instance_graph[v1][rand].second;
						break;
					}
				}

				if (left_change_times % 2 == 0) { // first increase
					double new_ec = selected_edge_weight * (1 + weightChange_ratio);
					if (new_ec > 1e6) {
						continue;
					}
					else {
						left_change_times--;
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						selected_edges.push_back(selected_edge);
					}
				}
				else { // then decrease
					double new_ec = selected_edge_weight * (1 - weightChange_ratio);
					if (new_ec < 1e-2) {
						continue;
					}
					else {
						left_change_times--;
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						selected_edges.push_back(selected_edge);
					}

				}
			}

			/*2014+2019*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.size();
				initialize_global_values_dynamic(V, thread_num);

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					if (k % 2 == 0) { // increase
						double new_ec = selected_edge_weight * (1 + weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						try {
							auto begin = std::chrono::high_resolution_clock::now();
							WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic, max_Maintain_time);
							//WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, max_Maintain_time);
							_20142019_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
						}
						catch (string s) {
							_20142019_time[k] = INT_MAX;
						}
					}
					else {
						double new_ec = selected_edge_weight * (1 - weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightDecrease2014(instance_graph, mm, selected_edge.first, selected_edge.second, new_ec, pool_dynamic, results_dynamic);
						_20142019_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					}
				}
			}

			/*newDE2021IN*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.size();
				initialize_global_values_dynamic(V, thread_num);

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					if (k % 2 == 0) { // increase
						double new_ec = selected_edge_weight * (1 + weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightIncrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic);
						//WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, max_Maintain_time);
						_newDE2021IN_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					}
					else {
						double new_ec = selected_edge_weight * (1 - weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightDecreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, new_ec, pool_dynamic, results_dynamic);
						_newDE2021IN_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					}
				}
			}

			/*new*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.size();
				initialize_global_values_dynamic(V, thread_num);

				L_bit_size_initial = mm.compute_L_bit_size();
				PPR_bit_size_initial = mm.compute_PPR_bit_size();

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					if (k % 2 == 0) { // increase
						double new_ec = selected_edge_weight * (1 + weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightIncreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic);
						_newDEIN_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					}
					else {
						double new_ec = selected_edge_weight * (1 - weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightDecreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, new_ec, pool_dynamic, results_dynamic);
						_newDEIN_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					}
				}

				L_bit_size_afterM1 = mm.compute_L_bit_size();
				PPR_bit_size_afterM1 = mm.compute_PPR_bit_size();

				int total_change_times = graph_v_of_v_idealID_total_edge_num(instance_graph) / 10;
				//if (data_name == "google" || data_name == "youtube" || data_name == "skitter") {
				//	total_change_times = 1e5;
				//}

				cout << "here" << endl;

				/*total_change_times-change_times changes*/
				auto gg = instance_graph;
				vector<pair<int, int>> edge_pool;
				for (int i = 0; i < V; i++) {
					for (auto adj : instance_graph[i]) {
						if (i < adj.first && adj.second >= 0.1 && adj.second <= 1e5) {
							edge_pool.push_back({ i, adj.first });
						}
					}
				}
				int left_change_times = total_change_times - change_times;
				while (left_change_times) {
					boost::range::random_shuffle(edge_pool);
					vector<pair<int, int>> new_edge_pool;
					for (auto e : edge_pool) {
						pair<int, int> selected_edge = e;
						double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
						selected_edges.push_back(selected_edge);
						if (left_change_times % 2 == 0) { // first increase
							double new_ec = selected_edge_weight * (1 + weightChange_ratio);
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
							if (new_ec >= 0.1 && new_ec <= 1e5) {
								new_edge_pool.push_back(selected_edge);
							}
						}
						else { // then decrease
							double new_ec = selected_edge_weight * (1 - weightChange_ratio);
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
							selected_edges.push_back(selected_edge);
							if (new_ec >= 0.1 && new_ec <= 1e5) {
								new_edge_pool.push_back(selected_edge);
							}
						}
						left_change_times--;
						if (left_change_times == 0) {
							break;
						}
					}
					edge_pool = new_edge_pool;
				}

				cout << "here1" << endl;

				instance_graph = gg;
				for (int k = change_times; k < total_change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					if (k % 2 == 0) { // increase
						double new_ec = selected_edge_weight * (1 + weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						WeightIncreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic);
					}
					else {
						double new_ec = selected_edge_weight * (1 - weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						WeightDecreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, new_ec, pool_dynamic, results_dynamic);
					}
				}

				L_bit_size_afterM2 = mm.compute_L_bit_size();
				PPR_bit_size_afterM2 = mm.compute_PPR_bit_size();

				ThreadPool pool_dynamic2(80);
				std::vector<std::future<int>> results_dynamic2;
				auto begin = std::chrono::high_resolution_clock::now();
				clean_L_dynamic(mm.L, mm.PPR, pool_dynamic2, results_dynamic2, thread_num);
				cleanL_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
				L_bit_size_afterClean = mm.compute_L_bit_size();

				begin = std::chrono::high_resolution_clock::now();
				clean_PPR(instance_graph, mm.L, mm.PPR, pool_dynamic2, results_dynamic2, thread_num);
				cleanPPR_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
				PPR_bit_size_afterClean = mm.compute_PPR_bit_size();


				mm.clear_labels();
				graph_hash_of_mixed_weighted g = graph_v_of_v_idealID_to_graph_hash_of_mixed_weighted(instance_graph);
				begin = std::chrono::high_resolution_clock::now();
				PLL_dynamic(g, instance_graph.size() + 1, thread_num, mm);
				clean_L_dynamic(mm.L, mm.PPR, pool_dynamic2, results_dynamic2, thread_num);
				rege_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
			}
		}

		double avg_2014DE_time = 0, avg_2019IN_time = 0, avg_2021DE_time = 0, avg_2021DE_query_times = 0, avg_2021IN_time = 0, avg_2021IN_query_times = 0,
			avg_newDE_time = 0, avg_newDE_query_times = 0, avg_newIN_time = 0, avg_newIN_query_times = 0, avg_20142019_time = 0, avg_newDE2021IN_time = 0, avg_newDEIN_time = 0;
		for (int k = 0; k < change_times; k++) {
			outputFile << _2014DE_time[k] << "," << _2019IN_time[k] << "," << _2021DE_time[k] << "," << _2021DE_query_times[k] << "," << _2021IN_time[k] << "," << _2021IN_query_times[k] << "," <<
				_newDE_time[k] << "," << _newDE_query_times[k] << "," << _newIN_time[k] << "," << _newIN_query_times[k] << "," << _20142019_time[k] << ","
				<< _newDE2021IN_time[k] << "," << _newDEIN_time[k] << "," <<
				L_bit_size_initial << "," << PPR_bit_size_initial / L_bit_size_initial << "," << L_bit_size_afterM1 / L_bit_size_initial << "," << PPR_bit_size_afterM1 / L_bit_size_initial << ","
				<< L_bit_size_afterM2 / L_bit_size_initial << "," << PPR_bit_size_afterM2 / L_bit_size_initial << ","
				<< L_bit_size_afterClean / L_bit_size_initial << "," << PPR_bit_size_afterClean / L_bit_size_initial << "," << cleanL_time << "," << cleanPPR_time << "," << rege_time << endl;
			avg_2014DE_time += _2014DE_time[k] / change_times;
			avg_2019IN_time += _2019IN_time[k] / change_times;
			avg_2021DE_time += _2021DE_time[k] / change_times;
			avg_2021DE_query_times += _2021DE_query_times[k] / change_times;
			avg_2021IN_time += _2021IN_time[k] / change_times;
			avg_2021IN_query_times += _2021IN_query_times[k] / change_times;
			avg_newDE_time += _newDE_time[k] / change_times;
			avg_newDE_query_times += _newDE_query_times[k] / change_times;
			avg_newIN_time += _newIN_time[k] / change_times;
			avg_newIN_query_times += _newIN_query_times[k] / change_times;
			avg_20142019_time += _20142019_time[k] / change_times;
			avg_newDE2021IN_time += _newDE2021IN_time[k] / change_times;
			avg_newDEIN_time += _newDEIN_time[k] / change_times;
		}
		outputFile << avg_2014DE_time << "," << avg_2019IN_time << "," << avg_2021DE_time << "," << avg_2021DE_query_times << "," << avg_2021IN_time << "," << avg_2021IN_query_times << "," <<
			avg_newDE_time << "," << avg_newDE_query_times << "," << avg_newIN_time << "," << avg_newIN_query_times << "," << avg_20142019_time << ","
			<< avg_newDE2021IN_time << "," << avg_newDEIN_time << "," <<
			L_bit_size_initial << "," << PPR_bit_size_initial / L_bit_size_initial << "," << L_bit_size_afterM1 / L_bit_size_initial << "," << PPR_bit_size_afterM1 / L_bit_size_initial << ","
			<< L_bit_size_afterM2 / L_bit_size_initial << "," << PPR_bit_size_afterM2 / L_bit_size_initial << ","
			<< L_bit_size_afterClean / L_bit_size_initial << "," << PPR_bit_size_afterClean / L_bit_size_initial << "," << cleanL_time << "," << cleanPPR_time << "," << rege_time << endl;

		outputFile.close(); // without this, multiple files cannot be successfully created
	}
}

void exp() {

	vector<string> data_names = { "astro", "condmat", "github", "google", "youtube", "skitter" };
	int change_times = 20, thread_num = 80;
	double max_Maintain_time = 600;

	/*weightChange_ratio 1*/
	if (1) {
		double weightChange_ratio = 0.8;
		for (auto data_name : data_names) {
			exp_element(data_name, weightChange_ratio, change_times, max_Maintain_time, thread_num);
		}
	}

	/*weightChange_ratio 2*/
	if (1) {
		double weightChange_ratio = 0.2;
		for (auto data_name : data_names) {
			exp_element(data_name, weightChange_ratio, change_times, max_Maintain_time, thread_num);
		}
	}
}

int main()
{
	cout << "Start running..." << endl;
	auto begin = std::chrono::high_resolution_clock::now();
	/*the two values below are for #include <graph_hash_of_mixed_weighted.h>*/
	graph_hash_of_mixed_weighted_turn_on_value = 1e3;
	graph_hash_of_mixed_weighted_turn_off_value = 1e1;

	exp();

	auto end = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	cout << "END    runningtime: " << runningtime << "s" << endl;
}