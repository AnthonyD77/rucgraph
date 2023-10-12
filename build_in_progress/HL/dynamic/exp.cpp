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

void exp_element1(string data_name, double weightChange_ratio, int change_times, double max_Maintain_time, int thread_num) {

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

	for (int type = 0; type < 2; type++) {

		string weight_type;
		if (type == 0) {
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
		string file_name = "exp_" + data_name + "_T_" + to_string(thread_num) + "_changeRatio_" + to_string((int)(weightChange_ratio * 100)) + "_" + weight_type + ".csv";
		cout << file_name << endl;
		outputFile.open(file_name);

		outputFile << "2014DE_time,2019IN_time,2021DE_time,2021DE_query_times,2021IN_time,2021IN_query_times,newDE_time,newDE_query_times,newIN_time,newIN_query_times," <<
			"2014+2019_time,2021DE2021IN_time,2021DEnewIN_time,newDE2021IN_time,newDEnewIN_time," <<
			"L_bit_size_initial(1),PPR_bit_size_initial,L_bit_size_afterM1,PPR_bit_size_afterM1,L_bit_size_afterClean1,PPR_bit_size_afterClean1,cleanL_time1,cleanPPR_time1,rege_time1," <<
			"L_bit_size_afterM2,PPR_bit_size_afterM2,L_bit_size_afterClean2,PPR_bit_size_afterClean2,cleanL_time2,cleanPPR_time2,rege_time2" << endl;

		int half_change_times = change_times / 2;
		vector<double> _2014DE_time(half_change_times, 0), _2019IN_time(half_change_times, 0), _2021DE_time(half_change_times, 0), _2021DE_query_times(half_change_times, 0), _2021IN_time(half_change_times, 0), _2021IN_query_times(half_change_times, 0),
			_newDE_time(half_change_times, 0), _newDE_query_times(half_change_times, 0), _newIN_time(half_change_times, 0), _newIN_query_times(half_change_times, 0),
			_20142019_time(half_change_times, 0), _2021DE2021IN_time(half_change_times, 0), _2021DEnewIN_time(half_change_times, 0), _newDE2021IN_time(half_change_times, 0), _newDEnewIN_time(half_change_times, 0);
		double L_bit_size_initial = 0, PPR_bit_size_initial = 0, L_bit_size_afterM1 = 0, PPR_bit_size_afterM1 = 0, L_bit_size_afterClean1 = 0, PPR_bit_size_afterClean1 = 0, cleanL_time1 = 0, cleanPPR_time1 = 0, rege_time1 = 0,
			L_bit_size_afterM2 = 0, PPR_bit_size_afterM2 = 0, L_bit_size_afterClean2 = 0, PPR_bit_size_afterClean2 = 0, cleanL_time2 = 0, cleanPPR_time2 = 0, rege_time2 = 0;

		/*mixed*/
		if (1) {
			instance_graph = instance_graph_initial;
			int V = instance_graph.size();
			vector<pair<int, int>>().swap(selected_edges);
			int left_change_times = change_times;
			while (left_change_times) {
				vector<pair<int, int>> edge_pool;
				for (int i = 0; i < V; i++) {
					for (auto adj : instance_graph[i]) {
						if (i < adj.first && adj.second >= 0.1 && adj.second <= 1e5) {
							edge_pool.push_back({ i, adj.first });
						}
					}
				}
				boost::range::random_shuffle(edge_pool);
				for (auto e : edge_pool) {
					pair<int, int> selected_edge = e;
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					selected_edges.push_back(selected_edge);
					if (left_change_times % 2 == 0) { // first increase
						double new_ec = selected_edge_weight * (1 + weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
					}
					else { // then decrease
						double new_ec = selected_edge_weight * (1 - weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						selected_edges.push_back(selected_edge);
					}
					left_change_times--;
					if (left_change_times == 0) {
						break;
					}
				}
			}

			cout << "step 1" << endl;

			/*2014+2019*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				initialize_global_values_dynamic(V, thread_num);

				for (int k = 0; k < change_times; k++) {
					cout << "k " << k << endl;
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					if (k % 2 == 0) { // increase
						cout << "step 1.1" << endl;
						double new_ec = selected_edge_weight * (1 + weightChange_ratio);
						auto mm_temp = mm;
						auto graph_temp = instance_graph;
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						try {
							auto begin = std::chrono::high_resolution_clock::now();
							WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic, max_Maintain_time);
							_2019IN_time[k / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
						}
						catch (string s) {
							instance_graph = graph_temp;
							mm = mm_temp; // WeightIncrease2019 may leave too many incorrect labels 
							_2019IN_time[k / 2] = INT_MAX;
						}
						cout << "step 1.2" << endl;
					}
					else {
						cout << "step 1.3" << endl;
						double new_ec = selected_edge_weight * (1 - weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightDecrease2014(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
						_2014DE_time[(k - 1) / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
						_20142019_time[(k - 1) / 2] = (_2019IN_time[(k - 1) / 2] + _2014DE_time[(k - 1) / 2]) / 2;
						cout << "step 1.4" << endl;
					}
				}
			}

			cout << "step 2" << endl;

			/*2021DE2021IN*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				initialize_global_values_dynamic(V, thread_num);

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					if (k % 2 == 0) { // increase
						global_query_times = 0;
						double new_ec = selected_edge_weight * (1 + weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightIncrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
						_2021IN_time[k / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
						_2021IN_query_times[k / 2] = global_query_times;
					}
					else {
						global_query_times = 0;
						double new_ec = selected_edge_weight * (1 - weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightDecrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
						_2021DE_time[(k - 1) / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
						_2021DE_query_times[(k - 1) / 2] = global_query_times;
						_2021DE2021IN_time[(k - 1) / 2] = (_2021IN_time[(k - 1) / 2] + _2021DE_time[(k - 1) / 2]) / 2;
					}
				}
			}

			cout << "step 3" << endl;

			/*new*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				initialize_global_values_dynamic(V, thread_num);

				L_bit_size_initial = mm.compute_L_bit_size();
				PPR_bit_size_initial = mm.compute_PPR_bit_size();

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					if (k % 2 == 0) { // increase
						global_query_times = 0;
						double new_ec = selected_edge_weight * (1 + weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightIncreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
						_newIN_time[k / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
						_newIN_query_times[k / 2] = global_query_times;
					}
					else {
						global_query_times = 0;
						double new_ec = selected_edge_weight * (1 - weightChange_ratio);
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightDecreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
						_newDE_time[(k - 1) / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
						_newDE_query_times[(k - 1) / 2] = global_query_times;
						_newDEnewIN_time[(k - 1) / 2] = (_newIN_time[(k - 1) / 2] + _newDE_time[(k - 1) / 2]) / 2;
						_2021DEnewIN_time[(k - 1) / 2] = (_newIN_time[(k - 1) / 2] + _2021DE_time[(k - 1) / 2]) / 2;
						_newDE2021IN_time[(k - 1) / 2] = (_2021IN_time[(k - 1) / 2] + _newDE_time[(k - 1) / 2]) / 2;
					}
				}

				cout << "step 4" << endl;
				mm_initial.clear_labels(); // to save RAM

				if (1) {

					int total_change_times = 1e4;

					/*total_change_times-change_times changes*/
					instance_graph = instance_graph_initial;
					vector<pair<int, int>>().swap(selected_edges);
					int left_change_times = total_change_times;
					while (left_change_times) {
						vector<pair<int, int>> edge_pool;
						for (int i = 0; i < V; i++) {
							for (auto adj : instance_graph[i]) {
								if (i < adj.first && adj.second >= 0.1 && adj.second <= 1e5) {
									edge_pool.push_back({ i, adj.first });
								}
							}
						}
						boost::range::random_shuffle(edge_pool);
						for (auto e : edge_pool) {
							pair<int, int> selected_edge = e;
							double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
							selected_edges.push_back(selected_edge);
							if (left_change_times % 2 == 0) { // first increase
								double new_ec = selected_edge_weight * (1 + weightChange_ratio);
								graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
							}
							else { // then decrease
								double new_ec = selected_edge_weight * (1 - weightChange_ratio);
								graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
								selected_edges.push_back(selected_edge);
							}
							left_change_times--;
							if (left_change_times == 0) {
								break;
							}
						}
					}

					instance_graph = instance_graph_initial;
					mm.clear_labels();
					binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm.PPR);
					binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm.L);
					initialize_global_values_dynamic(V, thread_num);

					for (int k = 0; k < total_change_times; k++) {
						pair<int, int> selected_edge = selected_edges[k];
						double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
						if (k % 2 == 0) { // increase
							double new_ec = selected_edge_weight * (1 + weightChange_ratio);
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
							WeightIncreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
						}
						else {
							double new_ec = selected_edge_weight * (1 - weightChange_ratio);
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
							WeightDecreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
						}
					}

					L_bit_size_afterM1 = mm.compute_L_bit_size();
					PPR_bit_size_afterM1 = mm.compute_PPR_bit_size();

					auto begin = std::chrono::high_resolution_clock::now();
					clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic, thread_num);
					cleanL_time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					L_bit_size_afterClean1 = mm.compute_L_bit_size();

					begin = std::chrono::high_resolution_clock::now();
					clean_PPR(instance_graph, mm.L, mm.PPR, pool_dynamic, results_dynamic, thread_num);
					cleanPPR_time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					PPR_bit_size_afterClean1 = mm.compute_PPR_bit_size();

					mm.clear_labels();
					begin = std::chrono::high_resolution_clock::now();
					graph_hash_of_mixed_weighted g = graph_v_of_v_idealID_to_graph_hash_of_mixed_weighted(instance_graph);
					PLL_dynamic(g, instance_graph.size() + 1, thread_num, mm);
					clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic, thread_num);
					rege_time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
				}

				cout << "step 5" << endl;

				/*skitter needs >1T RAM for 5e4, needs to swap() L and PPR*/
				if (0) {

					int total_change_times = 5e4;

					/*total_change_times-change_times changes*/
					instance_graph = instance_graph_initial;
					vector<pair<int, int>>().swap(selected_edges);
					int left_change_times = total_change_times;
					while (left_change_times) {
						vector<pair<int, int>> edge_pool;
						for (int i = 0; i < V; i++) {
							for (auto adj : instance_graph[i]) {
								if (i < adj.first && adj.second >= 0.1 && adj.second <= 1e5) {
									edge_pool.push_back({ i, adj.first });
								}
							}
						}
						boost::range::random_shuffle(edge_pool);
						for (auto e : edge_pool) {
							pair<int, int> selected_edge = e;
							double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
							selected_edges.push_back(selected_edge);
							if (left_change_times % 2 == 0) { // first increase
								double new_ec = selected_edge_weight * (1 + weightChange_ratio);
								graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
							}
							else { // then decrease
								double new_ec = selected_edge_weight * (1 - weightChange_ratio);
								graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
								selected_edges.push_back(selected_edge);
							}
							left_change_times--;
							if (left_change_times == 0) {
								break;
							}
						}
					}

					instance_graph = instance_graph_initial;
					mm.clear_labels();
					binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm.PPR);
					binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm.L);
					initialize_global_values_dynamic(V, thread_num);

					for (int k = 0; k < total_change_times; k++) {
						pair<int, int> selected_edge = selected_edges[k];
						double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
						if (k % 2 == 0) { // increase
							double new_ec = selected_edge_weight * (1 + weightChange_ratio);
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
							WeightIncreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
						}
						else {
							double new_ec = selected_edge_weight * (1 - weightChange_ratio);
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
							WeightDecreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
						}
					}

					L_bit_size_afterM2 = mm.compute_L_bit_size();
					PPR_bit_size_afterM2 = mm.compute_PPR_bit_size();

					auto begin = std::chrono::high_resolution_clock::now();
					clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic, thread_num);
					cleanL_time2 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					L_bit_size_afterClean2 = mm.compute_L_bit_size();

					begin = std::chrono::high_resolution_clock::now();
					clean_PPR(instance_graph, mm.L, mm.PPR, pool_dynamic, results_dynamic, thread_num);
					cleanPPR_time2 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					PPR_bit_size_afterClean2 = mm.compute_PPR_bit_size();

					mm.clear_labels();
					begin = std::chrono::high_resolution_clock::now();
					graph_hash_of_mixed_weighted g = graph_v_of_v_idealID_to_graph_hash_of_mixed_weighted(instance_graph);
					PLL_dynamic(g, instance_graph.size() + 1, thread_num, mm);
					clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic, thread_num);
					clean_PPR(instance_graph, mm.L, mm.PPR, pool_dynamic, results_dynamic, thread_num);
					rege_time2 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
				}

				cout << "step 6" << endl;
			}
		}

		double avg_2014DE_time = 0, avg_2019IN_time = 0, avg_2021DE_time = 0, avg_2021DE_query_times = 0, avg_2021IN_time = 0, avg_2021IN_query_times = 0,
			avg_newDE_time = 0, avg_newDE_query_times = 0, avg_newIN_time = 0, avg_newIN_query_times = 0,
			avg_20142019_time = 0, avg_2021DE2021IN_time = 0, avg_2021DEnewIN_time = 0, avg_newDE2021IN_time = 0, avg_newDEnewIN_time = 0;
		for (int k = 0; k < half_change_times; k++) {
			outputFile << _2014DE_time[k] << "," << _2019IN_time[k] << "," << _2021DE_time[k] << "," << _2021DE_query_times[k] << "," << _2021IN_time[k] << "," << _2021IN_query_times[k] << "," <<
				_newDE_time[k] << "," << _newDE_query_times[k] << "," << _newIN_time[k] << "," << _newIN_query_times[k] << ","
				<< _20142019_time[k] << "," << _2021DE2021IN_time[k] << "," << _2021DEnewIN_time[k] << "," << _newDE2021IN_time[k] << "," << _newDEnewIN_time[k] << "," <<
				L_bit_size_initial << "," << PPR_bit_size_initial / L_bit_size_initial << "," << L_bit_size_afterM1 / L_bit_size_initial << "," << PPR_bit_size_afterM1 / L_bit_size_initial << "," <<
				L_bit_size_afterClean1 / L_bit_size_initial << "," << PPR_bit_size_afterClean1 / L_bit_size_initial << "," << cleanL_time1 << "," << cleanPPR_time1 << "," << rege_time1 << "," <<
				L_bit_size_afterM2 / L_bit_size_initial << "," << PPR_bit_size_afterM2 / L_bit_size_initial << "," <<
				L_bit_size_afterClean2 / L_bit_size_initial << "," << PPR_bit_size_afterClean2 / L_bit_size_initial << "," << cleanL_time2 << "," << cleanPPR_time2 << "," << rege_time2 << endl;
			avg_2014DE_time += _2014DE_time[k] / half_change_times;
			avg_2019IN_time += _2019IN_time[k] / half_change_times;
			avg_2021DE_time += _2021DE_time[k] / half_change_times;
			avg_2021DE_query_times += _2021DE_query_times[k] / half_change_times;
			avg_2021IN_time += _2021IN_time[k] / half_change_times;
			avg_2021IN_query_times += _2021IN_query_times[k] / half_change_times;
			avg_newDE_time += _newDE_time[k] / half_change_times;
			avg_newDE_query_times += _newDE_query_times[k] / half_change_times;
			avg_newIN_time += _newIN_time[k] / half_change_times;
			avg_newIN_query_times += _newIN_query_times[k] / half_change_times;
			avg_20142019_time += _20142019_time[k] / half_change_times;
			avg_2021DE2021IN_time += _2021DE2021IN_time[k] / half_change_times;
			avg_2021DEnewIN_time += _2021DEnewIN_time[k] / half_change_times;
			avg_newDE2021IN_time += _newDE2021IN_time[k] / half_change_times;
			avg_newDEnewIN_time += _newDEnewIN_time[k] / half_change_times;
		}
		outputFile << avg_2014DE_time << "," << avg_2019IN_time << "," << avg_2021DE_time << "," << avg_2021DE_query_times << "," << avg_2021IN_time << "," << avg_2021IN_query_times << "," <<
			avg_newDE_time << "," << avg_newDE_query_times << "," << avg_newIN_time << "," << avg_newIN_query_times << "," <<
			avg_20142019_time << "," << avg_2021DE2021IN_time << "," << avg_2021DEnewIN_time << "," << avg_newDE2021IN_time << "," << avg_newDEnewIN_time << "," <<
			L_bit_size_initial << "," << PPR_bit_size_initial / L_bit_size_initial << "," << L_bit_size_afterM1 / L_bit_size_initial << "," << PPR_bit_size_afterM1 / L_bit_size_initial << "," <<
			L_bit_size_afterClean1 / L_bit_size_initial << "," << PPR_bit_size_afterClean1 / L_bit_size_initial << "," << cleanL_time1 << "," << cleanPPR_time1 << "," << rege_time1 << "," <<
			L_bit_size_afterM2 / L_bit_size_initial << "," << PPR_bit_size_afterM2 / L_bit_size_initial << "," <<
			L_bit_size_afterClean2 / L_bit_size_initial << "," << PPR_bit_size_afterClean2 / L_bit_size_initial << "," << cleanL_time2 << "," << cleanPPR_time2 << "," << rege_time2 << endl;

		outputFile.close(); // without this, multiple files cannot be successfully created
	}
}

class _edge {
public:
	int v1, v2;
	double ec;
};

void exp_element2(string data_name, int change_times, double max_Maintain_time, int thread_num) {

	ThreadPool pool_dynamic(thread_num);
	std::vector<std::future<int>> results_dynamic;

	ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);

	string path = "dynamicHL//";
	graph_v_of_v_idealID instance_graph;
	graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;

	for (int type = 0; type < 2; type++) {

		string weight_type;
		if (type == 0) {
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
		string file_name = "exp_" + data_name + "_T_" + to_string(thread_num) + "_DeleteInsert_" + weight_type + ".csv";
		cout << file_name << endl;
		outputFile.open(file_name);

		outputFile << "2014+2019_time,2021DE2021IN_time,newDEnewIN_time" << endl;

		vector<double> _20142019_time(change_times, 0), _2021DE2021IN_time(change_times, 0), _newDEnewIN_time(change_times, 0);

		/*mixed*/
		if (1) {
			double dummy_weight = 1e6;

			instance_graph = instance_graph_initial;
			int V = instance_graph.size();
			vector<vector<_edge>> changed_edges(change_times);
			if (1) {
				vector<int> vertex_pool;
				for (int i = 0; i < V; i++) {
					int degree = instance_graph[i].size();
					if (degree >= 2 && degree <= 5) {
						vertex_pool.push_back(i);
					}
				}
				boost::range::random_shuffle(vertex_pool);
				if (vertex_pool.size() < change_times) {
					cout << "vertex_pool.size() < left_change_times!" << endl;
					exit(1);
				}
				for (int i = 0; i < change_times; i++) {
					for (auto adj : instance_graph[vertex_pool[i]]) {
						changed_edges[i].push_back({ vertex_pool[i] , adj.first, adj.second });
					}
				}
			}

			cout << "step 1" << endl;

			/*2014+2019*/
			if (0) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				initialize_global_values_dynamic(V, thread_num);
				for (int k = 0; k < change_times; k++) {
					cout << "k " << k << endl;
					auto& _edges = changed_edges[k];
					auto mm_temp = mm;
					auto graph_temp = instance_graph;
					cout << "step 1.1" << endl;
					try {
						auto begin = std::chrono::high_resolution_clock::now();
						for (auto e : _edges) {
							graph_v_of_v_idealID_add_edge(instance_graph, e.v1, e.v2, dummy_weight);
							WeightIncrease2019(instance_graph, mm, e.v1, e.v2, e.ec, dummy_weight, pool_dynamic, results_dynamic, max_Maintain_time);
						}
						for (auto e : _edges) {
							graph_v_of_v_idealID_add_edge(instance_graph, e.v1, e.v2, e.ec);
							WeightDecrease2014(instance_graph, mm, e.v1, e.v2, dummy_weight, e.ec, pool_dynamic, results_dynamic);
						}
						_20142019_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9;
					}
					catch (string s) {
						instance_graph = graph_temp;
						mm = mm_temp; // WeightIncrease2019 may leave too many incorrect labels 
						_20142019_time[k] = INT_MAX;
					}
					cout << "step 1.2" << endl;
				}
			}

			cout << "step 2" << endl;

			/*2021DE2021IN*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				initialize_global_values_dynamic(V, thread_num);
				for (int k = 0; k < change_times; k++) {
					cout << "k " << k << endl;
					auto& _edges = changed_edges[k];
					auto begin = std::chrono::high_resolution_clock::now();
					cout << "_edges.size(): " << _edges.size() << endl;
					for (auto e : _edges) {
						cout << "IN " << e.v1 << " " << e.v2 << endl;
						graph_v_of_v_idealID_add_edge(instance_graph, e.v1, e.v2, dummy_weight);
						WeightIncrease2021(instance_graph, mm, e.v1, e.v2, e.ec, dummy_weight, pool_dynamic, results_dynamic);
						cout << "INend" << endl;
					}
					for (auto e : _edges) {
						cout << "DE " << e.v1 << " " << e.v2 << endl;
						graph_v_of_v_idealID_add_edge(instance_graph, e.v1, e.v2, e.ec);
						WeightDecrease2021(instance_graph, mm, e.v1, e.v2, dummy_weight, e.ec, pool_dynamic, results_dynamic);
						cout << "DEend" << endl;
					}
					_2021DE2021IN_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9;
				}
			}

			cout << "step 3" << endl;

			/*new*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				initialize_global_values_dynamic(V, thread_num);
				for (int k = 0; k < change_times; k++) {
					auto& _edges = changed_edges[k];
					auto begin = std::chrono::high_resolution_clock::now();
					for (auto e : _edges) {
						graph_v_of_v_idealID_add_edge(instance_graph, e.v1, e.v2, dummy_weight);
						WeightIncreaseMaintenance_improv(instance_graph, mm, e.v1, e.v2, e.ec, dummy_weight, pool_dynamic, results_dynamic);
					}
					for (auto e : _edges) {
						graph_v_of_v_idealID_add_edge(instance_graph, e.v1, e.v2, e.ec);
						WeightDecreaseMaintenance_improv(instance_graph, mm, e.v1, e.v2, dummy_weight, e.ec, pool_dynamic, results_dynamic);
					}
					_newDEnewIN_time[k] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9;
				}
			}

			cout << "step 4" << endl;
		}

		double avg_20142019_time = 0, avg_2021DE2021IN_time = 0, avg_newDEnewIN_time = 0;
		for (int k = 0; k < change_times; k++) {
			outputFile << _20142019_time[k] << "," << _2021DE2021IN_time[k] << "," << _newDEnewIN_time[k] << "," << endl;
			avg_20142019_time += _20142019_time[k] / change_times;
			avg_2021DE2021IN_time += _2021DE2021IN_time[k] / change_times;
			avg_newDEnewIN_time += _newDEnewIN_time[k] / change_times;
		}
		outputFile << avg_20142019_time << "," << avg_2021DE2021IN_time << "," << avg_newDEnewIN_time << endl;

		outputFile.close(); // without this, multiple files cannot be successfully created
	}
}

void exp_element3(string data_name, int change_times, double max_Maintain_time, int thread_num) {

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

	for (int type = 0; type < 2; type++) {

		string weight_type;
		if (type == 0) {
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
		string file_name = "exp_" + data_name + "_T_" + to_string(thread_num) + "_DeleteInsert_" + weight_type + ".csv";
		cout << file_name << endl;
		outputFile.open(file_name);

		outputFile << "2014+2019_time,2021DE2021IN_time,newDEnewIN_time" << endl;

		int half_change_times = change_times / 2;
		vector<double> _2014DE_time(half_change_times, 0), _2019IN_time(half_change_times, 0), _2021DE_time(half_change_times, 0), _2021DE_query_times(half_change_times, 0), _2021IN_time(half_change_times, 0), _2021IN_query_times(half_change_times, 0),
			_newDE_time(half_change_times, 0), _newDE_query_times(half_change_times, 0), _newIN_time(half_change_times, 0), _newIN_query_times(half_change_times, 0),
			_20142019_time(half_change_times, 0), _2021DE2021IN_time(half_change_times, 0), _2021DEnewIN_time(half_change_times, 0), _newDE2021IN_time(half_change_times, 0), _newDEnewIN_time(half_change_times, 0);

		/*mixed*/
		if (1) {
			double dummy_ec = 1e5, de_ec = 1;

			instance_graph = instance_graph_initial;
			int V = instance_graph.size();
			vector<pair<int, int>>().swap(selected_edges);
			int left_change_times = change_times;
			vector<pair<int, int>> edge_pool;
			for (int i = 0; i < V; i++) {
				for (auto adj : instance_graph[i]) {
					if (i < adj.first) {
						edge_pool.push_back({ i, adj.first });
					}
				}
			}
			boost::range::random_shuffle(edge_pool);
			boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(V - 1) };
			while (left_change_times) {
				if (left_change_times % 2 == 0) { // first increase
					selected_edges.push_back(edge_pool.back());
					edge_pool.pop_back();
				}
				else { // then decrease					
					while (1) {
						int v1 = dist(boost_random_time_seed), v2 = dist(boost_random_time_seed);
						if (graph_v_of_v_idealID_contain_edge(instance_graph, v1, v2) && v1 == v2) {
							continue;
						}
						selected_edges.push_back({ v1, v2 });
						break;
					}
				}
				left_change_times--;
			}

			cout << "step 1" << endl;

			/*2014+2019*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				initialize_global_values_dynamic(V, thread_num);

				for (int k = 0; k < change_times; k++) {
					cout << "k " << k << endl;
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					if (k % 2 == 0) { // increase
						cout << "step 1.1" << endl;
						double new_ec = dummy_ec;
						auto mm_temp = mm;
						auto graph_temp = instance_graph;
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						try {
							auto begin = std::chrono::high_resolution_clock::now();
							WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic, max_Maintain_time);
							_2019IN_time[k / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
						}
						catch (string s) {
							instance_graph = graph_temp;
							mm = mm_temp; // WeightIncrease2019 may leave too many incorrect labels 
							_2019IN_time[k / 2] = INT_MAX;
						}
						cout << "step 1.2" << endl;
					}
					else {
						cout << "step 1.3" << endl;
						double new_ec = de_ec;
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightDecrease2014(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
						_2014DE_time[(k - 1) / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
						_20142019_time[(k - 1) / 2] = (_2019IN_time[(k - 1) / 2] + _2014DE_time[(k - 1) / 2]) / 2;
						cout << "step 1.4" << endl;
					}
				}
			}

			cout << "step 2" << endl;

			/*2021DE2021IN*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				initialize_global_values_dynamic(V, thread_num);

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					if (k % 2 == 0) { // increase
						global_query_times = 0;
						double new_ec = dummy_ec;
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightIncrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
						_2021IN_time[k / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
						_2021IN_query_times[k / 2] = global_query_times;
					}
					else {
						global_query_times = 0;
						double new_ec = de_ec;
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightDecrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
						_2021DE_time[(k - 1) / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
						_2021DE_query_times[(k - 1) / 2] = global_query_times;
						_2021DE2021IN_time[(k - 1) / 2] = (_2021IN_time[(k - 1) / 2] + _2021DE_time[(k - 1) / 2]) / 2;
					}
				}
			}

			cout << "step 3" << endl;

			/*new*/
			if (1) {
				instance_graph = instance_graph_initial;
				mm = mm_initial;
				initialize_global_values_dynamic(V, thread_num);

				for (int k = 0; k < change_times; k++) {
					pair<int, int> selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					if (k % 2 == 0) { // increase
						global_query_times = 0;
						double new_ec = dummy_ec;
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightIncreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
						_newIN_time[k / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
						_newIN_query_times[k / 2] = global_query_times;
					}
					else {
						global_query_times = 0;
						double new_ec = de_ec;
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						auto begin = std::chrono::high_resolution_clock::now();
						WeightDecreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
						_newDE_time[(k - 1) / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
						_newDE_query_times[(k - 1) / 2] = global_query_times;
						_newDEnewIN_time[(k - 1) / 2] = (_newIN_time[(k - 1) / 2] + _newDE_time[(k - 1) / 2]) / 2;
						_2021DEnewIN_time[(k - 1) / 2] = (_newIN_time[(k - 1) / 2] + _2021DE_time[(k - 1) / 2]) / 2;
						_newDE2021IN_time[(k - 1) / 2] = (_2021IN_time[(k - 1) / 2] + _newDE_time[(k - 1) / 2]) / 2;
					}
				}
			}

			cout << "step 4" << endl;
		}

		double avg_20142019_time = 0, avg_2021DE2021IN_time = 0, avg_newDEnewIN_time = 0;
		for (int k = 0; k < half_change_times; k++) {
			outputFile << _20142019_time[k] << "," << _2021DE2021IN_time[k] << "," << _newDEnewIN_time[k] << endl;
			avg_20142019_time += _20142019_time[k] / half_change_times;
			avg_2021DE2021IN_time += _2021DE2021IN_time[k] / half_change_times;
			avg_newDEnewIN_time += _newDEnewIN_time[k] / half_change_times;
		}
		outputFile << avg_20142019_time << "," << avg_2021DE2021IN_time << "," << avg_newDEnewIN_time << endl;

		outputFile.close(); // without this, multiple files cannot be successfully created
	}
}



void exp() {

	vector<string> data_names = { "astro", "condmat", "github", "google", "youtube", "skitter" };
	int change_times = 100, thread_num = 80;
	double max_Maintain_time = 100;

	/*weightChange_ratio 1*/
	if (1) {
		double weightChange_ratio = 0.8;
		for (auto data_name : data_names) {
			exp_element1(data_name, weightChange_ratio, change_times, max_Maintain_time, thread_num);
		}
	}

	/*weightChange_ratio 2*/
	if (1) {
		double weightChange_ratio = 0.2;
		for (auto data_name : data_names) {
			exp_element1(data_name, weightChange_ratio, change_times, max_Maintain_time, thread_num);
		}
	}

	/*DeleteInsert*/
	if (1) {
		for (auto data_name : data_names) {
			exp_element3(data_name, change_times, max_Maintain_time, thread_num);
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
	//srand(time(NULL)); //  seed random number generator

	exp();

	auto end = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	cout << "END    runningtime: " << runningtime << "s" << endl;
}