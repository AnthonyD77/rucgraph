#include <iostream>
#include <fstream>
using namespace std;

// header files in the Boost library: https://www.boost.org/
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

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
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_binary_save_read.h>
#include <graph_hash_of_mixed_weighted/common_algorithms/graph_hash_of_mixed_weighted_shortest_paths.h>
#include <build_in_progress/HL/dynamic/clean_labels.h>
#include <text_mining/binary_save_read_vector_of_vectors.h>

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
		clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic);
		binary_save_PPR(path + s + "_PPR_random.bin", mm.PPR);
		binary_save_vector_of_vectors(path + s + "_L_random.bin", mm.L);
		outputFile.open(path + s + "_L_random_generation.txt");
		outputFile << "time_generate_labels=" << mm.time_generate_labels << "s" << endl;
		outputFile.close();

		g = graph_hash_of_mixed_weighted_binary_read(path + s + "_unique.bin");
		PLL_dynamic(g, g.hash_of_vectors.size(), thread_num, mm);
		clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic);
		binary_save_PPR(path + s + "_PPR_unique.bin", mm.PPR);
		binary_save_vector_of_vectors(path + s + "_L_unique.bin", mm.L);
		outputFile.open(path + s + "_L_unique_generation.txt");
		outputFile << "time_generate_labels=" << mm.time_generate_labels << "s" << endl;
		outputFile.close();
	}
}

void exp_element(string data_name, double weightChange_ratio, int change_times, double max_Maintain_time) {

	ThreadPool pool_dynamic(80);
	std::vector<std::future<int>> results_dynamic;

	ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);

	string path = "dynamicHL//";
	graph_hash_of_mixed_weighted instance_graph;
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
		graph_hash_of_mixed_weighted instance_graph_initial = graph_hash_of_mixed_weighted_binary_read(path + data_name + "_" + weight_type + ".bin");
		graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
		binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm_initial.PPR);
		binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm_initial.L);
		outputFile.open("exp_" + data_name + "_changeRatio_" + to_string((int)(weightChange_ratio * 100)) + "_" + weight_type + ".csv");

		outputFile << "2014DE_time,2019IN_time,2021DE_time,2021DE_query_times,2021IN_time,2021IN_query_times,newDE_time,newDE_query_times,newIN_time,newIN_query_times,2014+2019_time,newDEIN_time," <<
			"L_bit_size_initial(1),PPR_bit_size_initial,L_bit_size_afterM,PPR_bit_size_afterM,L_bit_size_afterClean,PPR_bit_size_afterClean,clean_time" << endl;

		double _2014DE_time = 0, _2019IN_time = 0, _2021DE_time = 0, _2021DE_query_times = 0, _2021IN_time = 0, _2021IN_query_times = 0,
			_newDE_time = 0, _newDE_query_times = 0, _newIN_time = 0, _newIN_query_times = 0, _20142019_time = 0, _newDEIN_time = 0,
			L_bit_size_initial = 0, PPR_bit_size_initial = 0, L_bit_size_afterM = 0, PPR_bit_size_afterM = 0, L_bit_size_afterClean = 0, PPR_bit_size_afterClean = 0, clean_time = 0;

		/*IN*/
		if (1) {
			instance_graph = instance_graph_initial;
			int V = instance_graph.hash_of_vectors.size();
			vector<pair<int, int>>().swap(selected_edges);
			int left_change_times = change_times;
			while (left_change_times) {
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

				double new_ec = selected_edge_weight * (1 + weightChange_ratio);
				if (new_ec > 1e6) {
					continue;
				}
				else {
					left_change_times--;
					graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
					selected_edges.push_back(selected_edge);
				}
			}

			/*2019*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.hash_of_vectors.size();
				initialize_global_values_dynamic(V, 1);

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					double new_ec = selected_edge_weight * (1 + weightChange_ratio);
					graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight						
					try {
						auto begin = std::chrono::high_resolution_clock::now();
						WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic, max_Maintain_time);
						//WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, max_Maintain_time);
						_2019IN_time += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9 / change_times; // s
					}
					catch (string s) {
						_2019IN_time = INT_MAX;
						break;
					}
				}
			}

			/*2021*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.hash_of_vectors.size();
				initialize_global_values_dynamic(V, 1);
				global_query_times = 0;

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					double new_ec = selected_edge_weight * (1 + weightChange_ratio);
					graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight						
					auto begin = std::chrono::high_resolution_clock::now();
					WeightIncrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight);
					_2021IN_time += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9 / change_times; // s
				}

				_2021IN_query_times = global_query_times / change_times;
			}

			/*new*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.hash_of_vectors.size();
				initialize_global_values_dynamic(V, 1);
				global_query_times = 0;

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					double new_ec = selected_edge_weight * (1 + weightChange_ratio);
					graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight						
					auto begin = std::chrono::high_resolution_clock::now();
					WeightIncreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight);
					_newIN_time += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9 / change_times; // s
				}

				_newIN_query_times = global_query_times / change_times;
			}
		}

		/*DE*/
		if (1) {
			instance_graph = instance_graph_initial;
			int V = instance_graph.hash_of_vectors.size();
			vector<pair<int, int>>().swap(selected_edges);
			int left_change_times = change_times;
			while (left_change_times) {
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

				double new_ec = selected_edge_weight * (1 - weightChange_ratio);
				if (new_ec < 1e-2) {
					continue;
				}
				else {
					left_change_times--;
					graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
					selected_edges.push_back(selected_edge);
				}
			}

			/*2014*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.hash_of_vectors.size();
				initialize_global_values_dynamic(V, 1);

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					double new_ec = selected_edge_weight * (1 - weightChange_ratio);
					graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
					auto begin = std::chrono::high_resolution_clock::now();
					WeightDecrease2014(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight);
					_2014DE_time += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9 / change_times; // s
				}
			}

			/*2021*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.hash_of_vectors.size();
				initialize_global_values_dynamic(V, 1);
				global_query_times = 0;

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					double new_ec = selected_edge_weight * (1 - weightChange_ratio);
					graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight					
					auto begin = std::chrono::high_resolution_clock::now();
					WeightDecrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight);
					_2021DE_time += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9 / change_times; // s
				}

				_2021DE_query_times = global_query_times / change_times;
			}

			/*new*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.hash_of_vectors.size();
				initialize_global_values_dynamic(V, 1);
				global_query_times = 0;

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					double new_ec = selected_edge_weight * (1 - weightChange_ratio);
					graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight						
					auto begin = std::chrono::high_resolution_clock::now();
					WeightDecreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight);
					_newDE_time += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9 / change_times; // s
				}

				_newDE_query_times = global_query_times / change_times;
			}
		}

		/*mixed*/
		if (1) {
			instance_graph = instance_graph_initial;
			int V = instance_graph.hash_of_vectors.size();
			vector<pair<int, int>>().swap(selected_edges);
			int left_change_times = change_times;
			while (left_change_times) {
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

				if (left_change_times % 2 == 0) { // first increase
					double new_ec = selected_edge_weight * (1 + weightChange_ratio);
					if (new_ec > 1e6) {
						continue;
					}
					else {
						left_change_times--;
						graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
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
						graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						selected_edges.push_back(selected_edge);
					}

				}
			}

			/*2014+2019*/
			if (1) {
				if (_2019IN_time == INT_MAX) {
					_20142019_time = INT_MAX;
				}

				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.hash_of_vectors.size();
				initialize_global_values_dynamic(V, 1);

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					if (k % 2 == 0) { // increase
						double new_ec = selected_edge_weight * (1 + weightChange_ratio);
						graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						try {
							auto begin = std::chrono::high_resolution_clock::now();
							WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic, max_Maintain_time);
							//WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, max_Maintain_time);
							_20142019_time += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9 / change_times; // s
						}
						catch (string s) {
							_2019IN_time = INT_MAX;
							_20142019_time = INT_MAX;
							break;
						}
					}
					else {
						double new_ec = selected_edge_weight * (1 - weightChange_ratio);
						graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightDecrease2014(instance_graph, mm, selected_edge.first, selected_edge.second, new_ec);
						_20142019_time += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9 / change_times; // s
					}
				}
			}

			/*new*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				int V = instance_graph.hash_of_vectors.size();
				initialize_global_values_dynamic(V, 1);

				L_bit_size_initial = mm.compute_L_bit_size();
				PPR_bit_size_initial = mm.compute_PPR_bit_size();

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					if (k % 2 == 0) { // increase
						double new_ec = selected_edge_weight * (1 + weightChange_ratio);
						graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightIncreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight);
						_newDEIN_time += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9 / change_times; // s
					}
					else {
						double new_ec = selected_edge_weight * (1 - weightChange_ratio);
						graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightDecreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, new_ec);
						_newDEIN_time += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9 / change_times; // s
					}
				}

				L_bit_size_afterM = mm.compute_L_bit_size();
				PPR_bit_size_afterM = mm.compute_PPR_bit_size();

				ThreadPool pool_dynamic2(80);
				std::vector<std::future<int>> results_dynamic2;
				auto begin = std::chrono::high_resolution_clock::now();
				clean_L_dynamic(mm.L, mm.PPR, pool_dynamic2, results_dynamic2);
				clean_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
				L_bit_size_afterClean = mm.compute_L_bit_size();
				PPR_bit_size_afterClean = mm.compute_PPR_bit_size();
			}

		}

		outputFile << _2014DE_time << "," << _2019IN_time << "," << _2021DE_time << "," << _2021DE_query_times << "," << _2021IN_time << "," << _2021IN_query_times << "," <<
			_newDE_time << "," << _newDE_query_times << "," << _newIN_time << "," << _newIN_query_times << "," << _20142019_time << "," << _newDEIN_time << "," <<
			L_bit_size_initial << "," << PPR_bit_size_initial / L_bit_size_initial << "," << L_bit_size_afterM / L_bit_size_initial << ","
			<< PPR_bit_size_afterM / L_bit_size_initial << "," << L_bit_size_afterClean / L_bit_size_initial << "," << PPR_bit_size_afterClean / L_bit_size_initial << "," << clean_time << endl;

		outputFile.close(); // without this, multiple files cannot be successfully created
	}
}

void exp() {

	vector<string> data_names = { "astro", "condmat", "github", "google", "youtube", "skitter" };
	int change_times = 100;
	double max_Maintain_time = 600;

	/*weightChange_ratio 1*/
	if (1) {
		double weightChange_ratio = 0.8;
		for (auto data_name : data_names) {
			exp_element(data_name, weightChange_ratio, change_times, max_Maintain_time);
		}
	}

	/*weightChange_ratio 2*/
	if (1) {
		double weightChange_ratio = 0.2;
		for (auto data_name : data_names) {
			exp_element(data_name, weightChange_ratio, change_times, max_Maintain_time);
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