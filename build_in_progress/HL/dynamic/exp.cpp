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
#include <build_in_progress/HL/dynamic/WeightIncrease2019.h>
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

	vector<string> data_names = { "google" };
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

		g = graph_hash_of_mixed_weighted_binary_read(path + s + "_Jacard.bin");
		PLL_dynamic(g, g.hash_of_vectors.size(), thread_num, mm);
		clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic);
		binary_save_PPR(path + s + "_PPR_Jacard.bin", mm.PPR);
		binary_save_vector_of_vectors(path + s + "_L_Jacard.bin", mm.L);
		outputFile.open(path + s + "_L_Jacard_generation.txt");
		outputFile << "time_generate_labels=" << mm.time_generate_labels << "s" << endl;
		outputFile.close();
	}
}

void exp_element(string data_name, double weightChange_ratio, int change_times, int multi_thread_num) {

	string save_file_name;
	ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);

	string path = "dynamicHL//";
	graph_hash_of_mixed_weighted instance_graph;
	graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;

	for (int i = 0; i < 4; i++) {

		string save_name, weight_type = "Jacard";
		int thread_num = 1;
		if (i % 2 == 1) {
			thread_num = multi_thread_num;
		}
		if (i < 2) {
			weight_type = "random";
		}
		save_name = "exp_" + data_name + "_t_" + to_string(thread_num) + "_changeRatio_" + to_string((int)(weightChange_ratio * 100)) + "_" + weight_type + ".csv";
		ThreadPool pool_dynamic(thread_num);
		std::vector<std::future<int>> results_dynamic;
		outputFile.open(save_name);
		outputFile << "2014DE_avg_time,2019IN_avg_time,L_bit_size_initial,PPR_bit_size_initial,L_bit_size_2014+2019," <<
			"2021DE_avg_time,2021IN_avg_time,L_bit_size_2021,PPR_bit_size_2021,"
			"newDE_avg_time,newIN_avg_time,L_bit_size_new,PPR_bit_size_new,clean_time,L_bit_size_new_clean,PPR_bit_size_new_clean" << endl;

		/*record edge changes: IN and DE mixed*/
		vector<pair<int, int>> selected_edges;
		if (i % 2 == 0) {
			instance_graph = graph_hash_of_mixed_weighted_binary_read(path + data_name + "_" + weight_type + ".bin");
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
		}

		for (int j = 0; j < 3; j++) {

			instance_graph = graph_hash_of_mixed_weighted_binary_read(path + data_name + "_" + weight_type + ".bin");
			binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm.PPR);
			binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm.L);

			double time_IN = 0, time_DE = 0, time_clean = 0;
			long long int L_bit_size_1 = 0, L_bit_size_2 = 0, L_bit_size_3 = 0, PPR_bit_size_1 = 0, PPR_bit_size_2 = 0, PPR_bit_size_3 = 0;

			L_bit_size_1 = mm.compute_L_bit_size();
			PPR_bit_size_1 = mm.compute_PPR_bit_size();

			for (int k = 0; k < change_times; k++) {
				pair<int, int> selected_edge = selected_edges[k];
				double selected_edge_weight = graph_hash_of_mixed_weighted_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
				if (k % 2 == 0) { // increase
					double new_ec = selected_edge_weight * (1 + weightChange_ratio);
					graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
					auto begin = std::chrono::high_resolution_clock::now();
					if (j == 0) {
						WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic);
					}
					else if (j == 1) {
						WeightIncrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic);
					}
					else {
						WeightIncreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, pool_dynamic, results_dynamic);
					}
					time_IN += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
				}
				else {
					double new_ec = selected_edge_weight * (1 - weightChange_ratio);
					graph_hash_of_mixed_weighted_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
					auto begin = std::chrono::high_resolution_clock::now();
					if (j == 0) {
						WeightDecrease2014(instance_graph, mm, selected_edge.first, selected_edge.second, new_ec, pool_dynamic, results_dynamic);
					}
					else if (j == 1) {
						WeightDecrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, new_ec, pool_dynamic, results_dynamic);
					}
					else {
						WeightDecreaseMaintenance_improv(instance_graph, mm, selected_edge.first, selected_edge.second, new_ec, pool_dynamic, results_dynamic);
					}
					time_DE += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
				}
			}

			L_bit_size_2 = mm.compute_L_bit_size();
			PPR_bit_size_2 = mm.compute_PPR_bit_size();

			if (j == 2) {
				auto begin = std::chrono::high_resolution_clock::now();
				clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic);
				time_clean = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
				L_bit_size_3 = mm.compute_L_bit_size();
				PPR_bit_size_3 = mm.compute_PPR_bit_size();
			}

			outputFile << "2014DE_avg_time,2019IN_avg_time,L_bit_size_initial,PPR_bit_size_initial,L_bit_size_2014+2019," <<
				"2021DE_avg_time,2021IN_avg_time,L_bit_size_2021,PPR_bit_size_2021,"
				"newDE_avg_time,newIN_avg_time,L_bit_size_new,PPR_bit_size_new,clean_time,L_bit_size_new_clean,PPR_bit_size_new_clean" << endl;

			if (j == 0) {
				outputFile << (double)time_DE / change_times * 2 << "," << (double)time_IN / change_times * 2 << "," << L_bit_size_1
					<< "," << PPR_bit_size_1 << "," << L_bit_size_2 << flush;
			}
			else if (j == 1) {
				outputFile << (double)time_DE / change_times * 2 << "," << (double)time_IN / change_times * 2 << "," << L_bit_size_2
					<< "," << PPR_bit_size_2 << flush;
			}
			else {
				outputFile << (double)time_DE / change_times * 2 << "," << (double)time_IN / change_times * 2 << "," << L_bit_size_2
					<< "," << PPR_bit_size_2 << "," << time_clean << "," << L_bit_size_3 << "," << PPR_bit_size_3 << endl;
			}

		}

	}
}

void exp() {

	vector<string> data_names = { "google" };
	int change_times = 1e4;
	int multi_thread_num = 50;

	/*weightChange_ratio 1*/
	if (1) {
		double weightChange_ratio = 0.1;
		for (auto data_name : data_names) {
			exp_element(data_name, weightChange_ratio, change_times, multi_thread_num);
		}
	}

	/*weightChange_ratio 2*/
	if (0) {
		double weightChange_ratio = 0.2;
		for (auto data_name : data_names) {
			exp_element(data_name, weightChange_ratio, change_times, multi_thread_num);
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