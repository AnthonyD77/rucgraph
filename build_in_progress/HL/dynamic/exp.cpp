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

	vector<string> data_names = { "amazon", "book" };
	string path = "dynamicHL//";
	int thread_num = 50;
	ThreadPool pool_dynamic(thread_num);
	std::vector<std::future<int>> results_dynamic;

	string save_file_name;
	ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);

	for (auto s : data_names) {
		if (1) {
			graph_hash_of_mixed_weighted g;
			graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
			g = graph_hash_of_mixed_weighted_binary_read(path + s + "_random.bin");
			PLL_dynamic(g, g.hash_of_vectors.size(), thread_num, mm);
			clean_L_dynamic(mm.L, mm.PPR, 80);
			binary_save_PPR(path + s + "_PPR_random.bin", mm.PPR);
			binary_save_vector_of_vectors(path + s + "_L_random.bin", mm.L);
			outputFile.open(path + s + "_L_random_generation.txt");
			outputFile << "time_generate_labels=" << mm.time_generate_labels << "s" << endl;
			outputFile.close();
		}

		if (1) {
			graph_hash_of_mixed_weighted g;
			graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
			g = graph_hash_of_mixed_weighted_binary_read(path + s + "_Jaccard.bin");
			PLL_dynamic(g, g.hash_of_vectors.size(), thread_num, mm);
			clean_L_dynamic(mm.L, mm.PPR, 80);
			binary_save_PPR(path + s + "_PPR_Jaccard.bin", mm.PPR);
			binary_save_vector_of_vectors(path + s + "_L_Jaccard.bin", mm.L);
			outputFile.open(path + s + "_L_Jaccard_generation.txt");
			outputFile << "time_generate_labels=" << mm.time_generate_labels << "s" << endl;
			outputFile.close();
		}

		if (1) {
			graph_hash_of_mixed_weighted g;
			graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
			g = graph_hash_of_mixed_weighted_binary_read(path + s + "_unique.bin");
			PLL_dynamic(g, g.hash_of_vectors.size(), thread_num, mm);
			clean_L_dynamic(mm.L, mm.PPR, 80);
			binary_save_PPR(path + s + "_PPR_unique.bin", mm.PPR);
			binary_save_vector_of_vectors(path + s + "_L_unique.bin", mm.L);
			outputFile.open(path + s + "_L_unique_generation.txt");
			outputFile << "time_generate_labels=" << mm.time_generate_labels << "s" << endl;
			outputFile.close();
		}
	}
}

class _edge {
public:
	int v1, v2;
	double ec;
};

void exp_main(string data_name, double weightChange_ratio, int change_times, double max_Maintain_time, int thread_num) {

	ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);

	string path = "dynamicHL//";

	for (int type = 0; type < 2; type++) {

		//if (data_name == "skitter" && type == 0) {
		//	continue;
		//}

		graph_v_of_v_idealID instance_graph;
		vector<_edge> selected_edges;

		string weight_type;
		if (type == 0) {
			weight_type = "Jaccard";
		}
		else if (type == 1) {
			weight_type = "random";
		}
		else {
			weight_type = "unique";
		}
		graph_hash_of_mixed_weighted instance_graph_initial_hash = graph_hash_of_mixed_weighted_binary_read(path + data_name + "_" + weight_type + ".bin");
		graph_v_of_v_idealID instance_graph_initial = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID_2(instance_graph_initial_hash, instance_graph_initial_hash.hash_of_vectors.size());
		graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
		binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm_initial.PPR);
		binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm_initial.L);
		string file_name = "exp_" + data_name + "_T_" + to_string(thread_num) + "_changeRatio_" + to_string((int)(weightChange_ratio * 100)) + "_" + weight_type + ".csv";
		cout << file_name << endl;
		outputFile.open(file_name);

		outputFile << "2014DE_time,2021DE_time,2021DE_query_times,newDE_time,newDE_query_times,DE_ratio,2019IN_time,2021IN_time,2021IN_query_times,newIN_time,newIN_query_times,IN_ratio," <<
			"2014+2019_time,2021DE2021IN_time,2021DEnewIN_time,newDE2021IN_time,newDEnewIN_time," <<
			"L_bit_size_0(1),PPR_size_0,L_size_1,PPR_size_1,L_size_1clean,PPR_size_1clean,cleanL_time1,cleanPPR_time1,rege_time1,L_size_2,PPR_size_2,L_size_2clean,PPR_size_2clean,cleanL_time2,cleanPPR_time2,rege_time2" << endl;

		int half_change_times = change_times / 2;
		vector<double> _2014DE_time(half_change_times, 0), _2019IN_time(half_change_times, 0), _2021DE_time(half_change_times, 0), _2021DE_query_times(half_change_times, 0), _2021IN_time(half_change_times, 0), _2021IN_query_times(half_change_times, 0),
			_newDE_time(half_change_times, 0), _newDE_query_times(half_change_times, 0), _newIN_time(half_change_times, 0), _newIN_query_times(half_change_times, 0),
			_20142019_time(half_change_times, 0), _2021DE2021IN_time(half_change_times, 0), _2021DEnewIN_time(half_change_times, 0), _newDE2021IN_time(half_change_times, 0), _newDEnewIN_time(half_change_times, 0);
		double L_size_0 = 0, PPR_size_0 = 0, L_size_1 = 0, PPR_size_1 = 0, L_size_1clean = 0, PPR_size_1clean = 0, cleanL_time1 = 0, cleanPPR_time1 = 0, rege_time1 = 0,
			L_size_2 = 0, PPR_size_2 = 0, L_size_2clean = 0, PPR_size_2clean = 0, cleanL_time2 = 0, cleanPPR_time2 = 0, rege_time2 = 0;

		/*mixed*/
		if (1) {
			double precision = std::pow(10, 3);
			int div = 50;

			instance_graph = instance_graph_initial;
			int V = instance_graph.size();
			vector<_edge>().swap(selected_edges);
			int left_change_times = change_times;
			while (left_change_times) {
				vector<pair<int, int>> edge_pool;
				for (int i = 0; i < V; i++) {
					for (auto adj : instance_graph[i]) {
						if (i < adj.first && adj.second >= 0.1 && adj.second <= 1e4 && instance_graph[i].size() > 1 && instance_graph[adj.first].size() > 1) {
							edge_pool.push_back({ i, adj.first });
						}
					}
				}
				boost::range::random_shuffle(edge_pool);
				for (auto e : edge_pool) {
					pair<int, int> selected_edge = e;
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					if (weightChange_ratio == 0) {
						if (left_change_times % 2 == 0) { // first increase
							boost::random::uniform_int_distribution<> dist{ static_cast<int>(selected_edge_weight * 2 * precision), static_cast<int>(selected_edge_weight * 10 * precision) };
							double new_ec = dist(boost_random_time_seed) / precision;
							selected_edges.push_back({ selected_edge.first, selected_edge.second, new_ec });
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						}
						else { // then decrease
							boost::random::uniform_int_distribution<> dist{ static_cast<int>(selected_edge_weight * 0.1 * precision), static_cast<int>(selected_edge_weight * 0.5 * precision) };
							double new_ec = dist(boost_random_time_seed) / precision;
							selected_edges.push_back({ selected_edge.first, selected_edge.second, new_ec });
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						}
					}
					else {
						if (left_change_times % 2 == 0) { // first increase
							double new_ec = selected_edge_weight * (1 + weightChange_ratio);
							selected_edges.push_back({ selected_edge.first, selected_edge.second, new_ec });
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
						}
						else { // then decrease
							double new_ec = selected_edge_weight * (1 - weightChange_ratio);
							selected_edges.push_back({ selected_edge.first, selected_edge.second, new_ec });
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
						}
					}
					left_change_times--;
					if (left_change_times == 0) {
						break;
					}
				}
			}

			cout << "step 1" << endl;

			if (data_name == "hyves" && type == 0) {
				div = 10;
			}
			else if (data_name == "skitter" && type == 0) {
				div = 6;
			}

			/*2014+2019*/
			if (1) {
				for (int j = 0; j < change_times / div; j++) {

					cout << "initialize L PPR" << endl;
					instance_graph = instance_graph_initial;
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm = mm_initial;
					initialize_global_values_dynamic(V, thread_num);

					for (int q = 0; q < div; q++) {
						int k = j * div + q;
						cout << "2014+2019 k " << k << endl;

						ThreadPool pool_dynamic(thread_num);
						std::vector<std::future<int>> results_dynamic;

						auto selected_edge = selected_edges[k];
						double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.v1, selected_edge.v2);
						if (k % 2 == 0) { // increase
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec); // increase weight
							try {
								auto begin = std::chrono::high_resolution_clock::now();
								WeightIncrease2019(instance_graph, mm, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic, max_Maintain_time);
								_2019IN_time[k / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
							}
							catch (string s) {
								instance_graph = instance_graph_initial;
								mm = mm_initial;
								initialize_global_values_dynamic(V, thread_num);
								_2019IN_time[k / 2] = INT_MAX;
							}
						}
						else {
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec); // decrease weight
							auto begin = std::chrono::high_resolution_clock::now();
							WeightDecrease2014(instance_graph, mm, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic);
							_2014DE_time[(k - 1) / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
							_20142019_time[(k - 1) / 2] = (_2019IN_time[(k - 1) / 2] + _2014DE_time[(k - 1) / 2]) / 2;
						}
					}
				}
			}

			cout << "step 2" << endl;

			/*2021DE2021IN*/
			if (1) {
				for (int j = 0; j < change_times / div; j++) {

					cout << "initialize L PPR" << endl;
					instance_graph = instance_graph_initial;
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm = mm_initial;
					initialize_global_values_dynamic(V, thread_num);

					for (int q = 0; q < div; q++) {
						int k = j * div + q;
						cout << "2021DE2021IN k " << k << endl;

						ThreadPool pool_dynamic(thread_num);
						std::vector<std::future<int>> results_dynamic;

						auto selected_edge = selected_edges[k];
						double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.v1, selected_edge.v2);
						if (k % 2 == 0) { // increase
							global_query_times = 0;
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec); // increase weight
							try {
								auto begin = std::chrono::high_resolution_clock::now();
								WeightIncrease2021(instance_graph, mm, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic);
								_2021IN_time[k / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
								_2021IN_query_times[k / 2] = global_query_times;
							}
							catch (string s) {
								instance_graph = instance_graph_initial;
								mm = mm_initial;
								initialize_global_values_dynamic(V, thread_num);
								_2021IN_time[k / 2] = INT_MAX;
								_2021IN_query_times[k / 2] = global_query_times;
							}
						}
						else {
							global_query_times = 0;
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec); // decrease weight
							try {
								auto begin = std::chrono::high_resolution_clock::now();
								WeightDecrease2021(instance_graph, mm, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic);
								_2021DE_time[(k - 1) / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
								_2021DE_query_times[(k - 1) / 2] = global_query_times;
								_2021DE2021IN_time[(k - 1) / 2] = (_2021IN_time[(k - 1) / 2] + _2021DE_time[(k - 1) / 2]) / 2;
							}
							catch (string s) {
								instance_graph = instance_graph_initial;
								mm = mm_initial;
								initialize_global_values_dynamic(V, thread_num);
								_2021DE_time[(k - 1) / 2] = INT_MAX; // s
								_2021DE_query_times[(k - 1) / 2] = global_query_times;
								_2021DE2021IN_time[(k - 1) / 2] = (_2021IN_time[(k - 1) / 2] + _2021DE_time[(k - 1) / 2]) / 2;
							}
						}
					}
				}
			}

			div = 50;

			cout << "step 3" << endl;

			/*new*/
			if (1) {
				for (int j = 0; j < change_times / div; j++) {

					cout << "initialize L PPR" << endl;
					instance_graph = instance_graph_initial;
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm = mm_initial;
					initialize_global_values_dynamic(V, thread_num);

					for (int q = 0; q < div; q++) {
						int k = j * div + q;
						cout << "new k " << k << endl;

						ThreadPool pool_dynamic(thread_num);
						std::vector<std::future<int>> results_dynamic;

						auto selected_edge = selected_edges[k];
						double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.v1, selected_edge.v2);
						if (k % 2 == 0) { // increase
							global_query_times = 0;
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec); // increase weight
							auto begin = std::chrono::high_resolution_clock::now();
							WeightIncreaseMaintenance_improv(instance_graph, mm, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic);
							_newIN_time[k / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
							_newIN_query_times[k / 2] = global_query_times;
						}
						else {
							global_query_times = 0;
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec);
							auto begin = std::chrono::high_resolution_clock::now();
							WeightDecreaseMaintenance_improv(instance_graph, mm, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic);
							_newDE_time[(k - 1) / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
							_newDE_query_times[(k - 1) / 2] = global_query_times;
							_newDEnewIN_time[(k - 1) / 2] = (_newIN_time[(k - 1) / 2] + _newDE_time[(k - 1) / 2]) / 2;
							_2021DEnewIN_time[(k - 1) / 2] = (_newIN_time[(k - 1) / 2] + _2021DE_time[(k - 1) / 2]) / 2;
							_newDE2021IN_time[(k - 1) / 2] = (_2021IN_time[(k - 1) / 2] + _newDE_time[(k - 1) / 2]) / 2;
						}
					}
				}

				cout << "step 4" << endl;

				ofstream outputFile2;
				outputFile2.precision(6);
				outputFile2.setf(ios::fixed);
				outputFile2.setf(ios::showpoint);
				outputFile2.open("temp_" + file_name);
				outputFile2 << "2014DE_time,2021DE_time,2021DE_query_times,newDE_time,newDE_query_times,DE_ratio,2019IN_time,2021IN_time,2021IN_query_times,newIN_time,newIN_query_times,IN_ratio," <<
					"2014+2019_time,2021DE2021IN_time,2021DEnewIN_time,newDE2021IN_time,newDEnewIN_time," <<
					"L_bit_size_initial(1),PPR_bit_size_initial,L_bit_size_afterM1,PPR_bit_size_afterM1,L_bit_size_afterClean1,PPR_bit_size_afterClean1,cleanL_time1,cleanPPR_time1,rege_time1" << endl;
				double avg_2014DE_time = 0, avg_2019IN_time = 0, avg_2021DE_time = 0, avg_2021DE_query_times = 0, avg_2021IN_time = 0, avg_2021IN_query_times = 0, avg_DEratio = 0, avg_INratio = 0,
					avg_newDE_time = 0, avg_newDE_query_times = 0, avg_newIN_time = 0, avg_newIN_query_times = 0,
					avg_20142019_time = 0, avg_2021DE2021IN_time = 0, avg_2021DEnewIN_time = 0, avg_newDE2021IN_time = 0, avg_newDEnewIN_time = 0;
				for (int k = 0; k < half_change_times; k++) {
					outputFile2 << _2014DE_time[k] << "," << _2021DE_time[k] << "," << _2021DE_query_times[k] << "," << _newDE_time[k] << "," << _newDE_query_times[k] << "," << _newDE_time[k] / _2021DE_time[k] << "," <<
						_2019IN_time[k] << "," << _2021IN_time[k] << "," << _2021IN_query_times[k] << "," << _newIN_time[k] << "," << _newIN_query_times[k] << "," << _newIN_time[k] / _2021IN_time[k] << "," <<
						_20142019_time[k] << "," << _2021DE2021IN_time[k] << "," << _2021DEnewIN_time[k] << "," << _newDE2021IN_time[k] << "," << _newDEnewIN_time[k] << endl;
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
					avg_DEratio += _newDE_time[k] / _2021DE_time[k] / half_change_times;
					avg_INratio += _newIN_time[k] / _2021IN_time[k] / half_change_times;
				}
				outputFile2 << avg_2014DE_time << "," << avg_2021DE_time << "," << avg_2021DE_query_times << "," << avg_newDE_time << "," << avg_newDE_query_times << "," << avg_DEratio << "," <<
					avg_2019IN_time << "," << avg_2021IN_time << "," << avg_2021IN_query_times << "," << avg_newIN_time << "," << avg_newIN_query_times << "," << avg_INratio << "," <<
					avg_20142019_time << "," << avg_2021DE2021IN_time << "," << avg_2021DEnewIN_time << "," << avg_newDE2021IN_time << "," << avg_newDEnewIN_time << endl;
				outputFile2.close(); // without this, multiple files cannot be successfully created	
			}
		}

		double avg_2014DE_time = 0, avg_2019IN_time = 0, avg_2021DE_time = 0, avg_2021DE_query_times = 0, avg_2021IN_time = 0, avg_2021IN_query_times = 0, avg_DEratio = 0, avg_INratio = 0,
			avg_newDE_time = 0, avg_newDE_query_times = 0, avg_newIN_time = 0, avg_newIN_query_times = 0,
			avg_20142019_time = 0, avg_2021DE2021IN_time = 0, avg_2021DEnewIN_time = 0, avg_newDE2021IN_time = 0, avg_newDEnewIN_time = 0;
		for (int k = 0; k < half_change_times; k++) {
			outputFile << _2014DE_time[k] << "," << _2021DE_time[k] << "," << _2021DE_query_times[k] << "," << _newDE_time[k] << "," << _newDE_query_times[k] << "," << _newDE_time[k] / _2021DE_time[k] << "," <<
				_2019IN_time[k] << "," << _2021IN_time[k] << "," << _2021IN_query_times[k] << "," << _newIN_time[k] << "," << _newIN_query_times[k] << "," << _newIN_time[k] / _2021IN_time[k] << "," <<
				_20142019_time[k] << "," << _2021DE2021IN_time[k] << "," << _2021DEnewIN_time[k] << "," << _newDE2021IN_time[k] << "," << _newDEnewIN_time[k] << "," <<
				L_size_0 << "," << PPR_size_0 / L_size_0 << "," << L_size_1 / L_size_0 << "," << PPR_size_1 / L_size_0 << "," <<
				L_size_1clean / L_size_0 << "," << PPR_size_1clean / L_size_0 << "," << cleanL_time1 << "," << cleanPPR_time1 << "," << rege_time1 << "," <<
				L_size_2 / L_size_0 << "," << PPR_size_2 / L_size_0 << "," << L_size_2clean / L_size_0 << "," << PPR_size_2clean / L_size_0 << "," << cleanL_time2 << "," << cleanPPR_time2 << "," << rege_time2 << endl;

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
			avg_DEratio += _newDE_time[k] / _2021DE_time[k] / half_change_times;
			avg_INratio += _newIN_time[k] / _2021IN_time[k] / half_change_times;
		}
		outputFile << avg_2014DE_time << "," << avg_2021DE_time << "," << avg_2021DE_query_times << "," << avg_newDE_time << "," << avg_newDE_query_times << "," << avg_DEratio << "," <<
			avg_2019IN_time << "," << avg_2021IN_time << "," << avg_2021IN_query_times << "," << avg_newIN_time << "," << avg_newIN_query_times << "," << avg_INratio << "," <<
			avg_20142019_time << "," << avg_2021DE2021IN_time << "," << avg_2021DEnewIN_time << "," << avg_newDE2021IN_time << "," << avg_newDEnewIN_time << "," <<
			L_size_0 << "," << PPR_size_0 / L_size_0 << "," << L_size_1 / L_size_0 << "," << PPR_size_1 / L_size_0 << "," <<
			L_size_1clean / L_size_0 << "," << PPR_size_1clean / L_size_0 << "," << cleanL_time1 << "," << cleanPPR_time1 << "," << rege_time1 << "," <<
			L_size_2 / L_size_0 << "," << PPR_size_2 / L_size_0 << "," << L_size_2clean / L_size_0 << "," << PPR_size_2clean / L_size_0 << "," << cleanL_time2 << "," << cleanPPR_time2 << "," << rege_time2 << endl;

		outputFile.close(); // without this, multiple files cannot be successfully created
	}
}

void exp_clean(string data_name, double weightChange_ratio, int change_times, double max_Maintain_time, int thread_num) {

	ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);

	string path = "dynamicHL//";

	for (int type = 0; type < 2; type++) {

		vector<_edge> selected_edges;

		string weight_type;
		if (type == 0) {
			weight_type = "Jaccard";
		}
		else if (type == 1) {
			weight_type = "random";
		}
		graph_hash_of_mixed_weighted instance_graph_initial_hash = graph_hash_of_mixed_weighted_binary_read(path + data_name + "_" + weight_type + ".bin");
		graph_v_of_v_idealID instance_graph = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID_2(instance_graph_initial_hash, instance_graph_initial_hash.hash_of_vectors.size());
		string file_name = "clean_exp_" + data_name + "_T_" + to_string(thread_num) + "_changeRatio_" + to_string((int)(weightChange_ratio * 100)) + "_" + weight_type + ".csv";
		cout << file_name << endl;
		outputFile.open(file_name);

		outputFile << "L_bit_size_0(1),PPR_size_0,L_size_1,PPR_size_1,L_size_1clean,PPR_size_1clean,cleanL_time1,cleanPPR_time1," <<
			"rege_time1,L_size_2,PPR_size_2,L_size_2clean,PPR_size_2clean,cleanL_time2,cleanPPR_time2,rege_time2,"
			"temp_size_1,temp_size_2,temp_size_3,temp_size_4,temp_size_5,temp_size_6,temp_size_7,temp_size_8,temp_size_9,"
			"temp_size_10,temp_size_11,temp_size_12,temp_size_13,temp_size_14,temp_size_15,temp_size_16,temp_size_17,temp_size_18" << endl;

		double L_size_0 = 0, PPR_size_0 = 0, L_size_1 = 0, PPR_size_1 = 0, L_size_1clean = 0, PPR_size_1clean = 0, cleanL_time1 = 0, cleanPPR_time1 = 0, rege_time1 = 0,
			L_size_2 = 0, PPR_size_2 = 0, L_size_2clean = 0, PPR_size_2clean = 0, cleanL_time2 = 0, cleanPPR_time2 = 0, rege_time2 = 0;


		double temp_size_1 = 0, temp_size_2 = 0, temp_size_3 = 0, temp_size_4 = 0, temp_size_5 = 0, temp_size_6 = 0, temp_size_7 = 0, temp_size_8 = 0,
			temp_size_9 = 0, temp_size_10 = 0, temp_size_11 = 0, temp_size_12 = 0, temp_size_13 = 0, temp_size_14 = 0, temp_size_15 = 0, temp_size_16 = 0, temp_size_17 = 0, temp_size_18 = 0;

		/*mixed*/
		if (1) {
			double precision = std::pow(10, 3);
			int div = 500;
			int V = instance_graph.size();
			initialize_global_values_dynamic(V, thread_num);

			int total_change_times1 = 1e3, total_change_times2 = 4e3;

			cout << "step 1" << endl;

			/*changes 1*/
			if (1) {

				/*total_change_times-change_times changes*/
				auto gg = instance_graph;
				vector<_edge>().swap(selected_edges);
				int left_change_times = total_change_times1;
				while (left_change_times) {
					vector<pair<int, int>> edge_pool;
					for (int i = 0; i < V; i++) {
						for (auto adj : instance_graph[i]) {
							if (i < adj.first && adj.second >= 0.1 && adj.second <= 1e4) {
								edge_pool.push_back({ i, adj.first });
							}
						}
					}
					boost::range::random_shuffle(edge_pool);
					for (auto e : edge_pool) {
						pair<int, int> selected_edge = e;
						double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
						if (weightChange_ratio == 0) {
							if (left_change_times % 2 == 0) { // first increase
								boost::random::uniform_int_distribution<> dist{ static_cast<int>(selected_edge_weight * 2 * precision), static_cast<int>(selected_edge_weight * 10 * precision) };
								double new_ec = dist(boost_random_time_seed) / precision;
								selected_edges.push_back({ selected_edge.first, selected_edge.second, new_ec });
								graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
							}
							else { // then decrease
								boost::random::uniform_int_distribution<> dist{ static_cast<int>(selected_edge_weight * 0.1 * precision), static_cast<int>(selected_edge_weight * 0.5 * precision) };
								double new_ec = dist(boost_random_time_seed) / precision;
								selected_edges.push_back({ selected_edge.first, selected_edge.second, new_ec });
								graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
							}
						}
						else {
							if (left_change_times % 2 == 0) { // first increase
								double new_ec = selected_edge_weight * (1 + weightChange_ratio);
								selected_edges.push_back({ selected_edge.first, selected_edge.second, new_ec });
								graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
							}
							else { // then decrease
								double new_ec = selected_edge_weight * (1 - weightChange_ratio);
								selected_edges.push_back({ selected_edge.first, selected_edge.second, new_ec });
								graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
							}
						}
						left_change_times--;
						if (left_change_times == 0) {
							break;
						}
					}
				}
				instance_graph = gg;

				cout << "step 2" << endl;

				for (int j = 0; j < total_change_times1 / div; j++) {
					ThreadPool pool_dynamic(thread_num);
					std::vector<std::future<int>> results_dynamic;

					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
					if (j == 0) {
						binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm_initial.PPR);
						binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm_initial.L);
						L_size_0 = mm_initial.compute_L_bit_size();
						PPR_size_0 = mm_initial.compute_PPR_bit_size();
						//outputFile << "L_size_0 " << L_size_0 << " PPR_size_0 " << PPR_size_0 << endl;
					}
					else {
						binary_read_vector_of_vectors("temp_L.bin", mm_initial.L);
						binary_read_PPR("temp_PPR.bin", mm_initial.PPR);
						std::remove("temp_L.bin");
						std::remove("temp_PPR.bin");
					}

					for (int q = 0; q < div; q++) {
						int k = j * div + q;
						cout << "new large k " << k << endl;
						auto selected_edge = selected_edges[k];
						double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.v1, selected_edge.v2);
						if (k % 2 == 0) { // increase
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec); // increase weight
							WeightIncreaseMaintenance_improv(instance_graph, mm_initial, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic);
						}
						else {
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec);
							WeightDecreaseMaintenance_improv(instance_graph, mm_initial, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic);
						}
					}

					binary_save_vector_of_vectors("temp_L.bin", mm_initial.L);
					binary_save_PPR("temp_PPR.bin", mm_initial.PPR);

					if (j == total_change_times1 / div - 1) {
						L_size_1 = mm_initial.compute_L_bit_size();
						PPR_size_1 = mm_initial.compute_PPR_bit_size();
						//outputFile << "L_size_1 " << L_size_1 << " PPR_size_1 " << PPR_size_1 << endl;
					}
				}

				cout << "step 3" << endl;
			}

			/*clean L*/
			if (1) {
				graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
				binary_read_vector_of_vectors("temp_L.bin", mm_initial.L);
				std::remove("temp_L.bin");

				temp_size_1 = mm_initial.compute_L_bit_size();

				cout << "step 4" << endl;

				auto begin = std::chrono::high_resolution_clock::now();
				clean_L_dynamic(mm_initial.L, mm_initial.PPR, thread_num);
				cleanL_time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s

				temp_size_2 = mm_initial.compute_L_bit_size();

				binary_save_vector_of_vectors("temp_L.bin", mm_initial.L);

				cout << "step 5" << endl;
			}

			/*re-ge PPR*/
			if (1) {
				graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
				binary_read_vector_of_vectors("temp_L.bin", mm_initial.L);
				std::remove("temp_L.bin");

				temp_size_3 = mm_initial.compute_L_bit_size();

				cout << "step 6" << endl;

				auto begin = std::chrono::high_resolution_clock::now();
				clean_PPR(instance_graph, mm_initial.L, mm_initial.PPR, thread_num);
				cleanPPR_time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s

				temp_size_4 = mm_initial.compute_PPR_bit_size();

				cout << "step 7" << endl;
			}

			/*re-ge*/
			if (1) {
				double time_rege = 0, time_cleanL = 0, time_clean_PPR = 0;

				if (1) {
					graph_hash_of_mixed_weighted g = graph_v_of_v_idealID_to_graph_hash_of_mixed_weighted(instance_graph);
					auto begin = std::chrono::high_resolution_clock::now();
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
					PLL_dynamic(g, instance_graph.size() + 1, thread_num, mm_initial);
					time_rege = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					binary_save_vector_of_vectors("temp_L.bin", mm_initial.L);

					temp_size_5 = mm_initial.compute_L_bit_size();
					temp_size_6 = mm_initial.compute_PPR_bit_size();
				}

				cout << "step 8" << endl;

				if (1) {
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
					binary_read_vector_of_vectors("temp_L.bin", mm_initial.L);
					std::remove("temp_L.bin");

					temp_size_7 = mm_initial.compute_L_bit_size();

					auto begin = std::chrono::high_resolution_clock::now();
					clean_L_dynamic(mm_initial.L, mm_initial.PPR, thread_num);
					time_cleanL = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					L_size_1clean = mm_initial.compute_L_bit_size();
					binary_save_vector_of_vectors("temp_L.bin", mm_initial.L);
					//outputFile << "L_size_1clean " << L_size_1clean << endl;
				}

				cout << "step 9" << endl;

				if (1) {
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
					binary_read_vector_of_vectors("temp_L.bin", mm_initial.L);
					std::remove("temp_L.bin");

					temp_size_8 = mm_initial.compute_L_bit_size();

					auto begin = std::chrono::high_resolution_clock::now();
					clean_PPR(instance_graph, mm_initial.L, mm_initial.PPR, thread_num);
					time_clean_PPR = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					PPR_size_1clean = mm_initial.compute_PPR_bit_size();
					binary_save_vector_of_vectors("temp_L.bin", mm_initial.L);
					binary_save_PPR("temp_PPR.bin", mm_initial.PPR);
					//outputFile << "PPR_size_1clean " << PPR_size_1clean << endl;
				}

				/*
				just parallel PLL_dynamic generates redundant 2-hop labels, and clean_L_dynamic is required to remove redundant 2-hop labels, and
				clean_PPR is further required to produce the correcponding PPR;

				the strategy of parallelizing PLL with clean_L_dynamic: K. Lakhotia, R. Kannan, Q. Dong, and V. Prasanna, ¡°Planting trees for scalable
				and efficient canonical hub labeling,¡± Proc. VLDB Endow. 13 (2019).
				*/

				rege_time1 = time_rege + time_cleanL + time_clean_PPR; // s

				initialize_global_values_dynamic(V, thread_num); // Qid_595 needs to be initialized

				cout << "step 10" << endl;
			}

			/*changes 2*/
			if (1) {

				/*total_change_times-change_times changes*/
				auto gg = instance_graph;
				vector<_edge>().swap(selected_edges);
				int left_change_times = total_change_times2;
				while (left_change_times) {
					vector<pair<int, int>> edge_pool;
					for (int i = 0; i < V; i++) {
						for (auto adj : instance_graph[i]) {
							if (i < adj.first && adj.second >= 0.1 && adj.second <= 1e4) {
								edge_pool.push_back({ i, adj.first });
							}
						}
					}
					boost::range::random_shuffle(edge_pool);
					for (auto e : edge_pool) {
						pair<int, int> selected_edge = e;
						double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
						if (weightChange_ratio == 0) {
							if (left_change_times % 2 == 0) { // first increase
								boost::random::uniform_int_distribution<> dist{ static_cast<int>(selected_edge_weight * 2 * precision), static_cast<int>(selected_edge_weight * 10 * precision) };
								double new_ec = dist(boost_random_time_seed) / precision;
								selected_edges.push_back({ selected_edge.first, selected_edge.second, new_ec });
								graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
							}
							else { // then decrease
								boost::random::uniform_int_distribution<> dist{ static_cast<int>(selected_edge_weight * 0.1 * precision), static_cast<int>(selected_edge_weight * 0.5 * precision) };
								double new_ec = dist(boost_random_time_seed) / precision;
								selected_edges.push_back({ selected_edge.first, selected_edge.second, new_ec });
								graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
							}
						}
						else {
							if (left_change_times % 2 == 0) { // first increase
								double new_ec = selected_edge_weight * (1 + weightChange_ratio);
								selected_edges.push_back({ selected_edge.first, selected_edge.second, new_ec });
								graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
							}
							else { // then decrease
								double new_ec = selected_edge_weight * (1 - weightChange_ratio);
								selected_edges.push_back({ selected_edge.first, selected_edge.second, new_ec });
								graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
							}
						}
						left_change_times--;
						if (left_change_times == 0) {
							break;
						}
					}
				}
				instance_graph = gg;

				for (int j = 0; j < total_change_times2 / div; j++) {
					ThreadPool pool_dynamic(thread_num);
					std::vector<std::future<int>> results_dynamic;

					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
					binary_read_vector_of_vectors("temp_L.bin", mm_initial.L);
					binary_read_PPR("temp_PPR.bin", mm_initial.PPR);
					std::remove("temp_L.bin");
					std::remove("temp_PPR.bin");
					if (j == 0) {
						temp_size_9 = mm_initial.compute_L_bit_size();
						temp_size_10 = mm_initial.compute_PPR_bit_size();
						//outputFile << "temp_size_9 " << temp_size_9 << " temp_size_10 " << temp_size_10 << endl;
					}

					for (int q = 0; q < div; q++) {
						int k = j * div + q;
						cout << "new large k " << k << endl;
						auto selected_edge = selected_edges[k];
						double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.v1, selected_edge.v2);
						if (k % 2 == 0) { // increase
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec); // increase weight
							WeightIncreaseMaintenance_improv(instance_graph, mm_initial, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic);
						}
						else {
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec);
							WeightDecreaseMaintenance_improv(instance_graph, mm_initial, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic);
						}
					}

					binary_save_vector_of_vectors("temp_L.bin", mm_initial.L);
					binary_save_PPR("temp_PPR.bin", mm_initial.PPR);

					if (j == total_change_times2 / div - 1) {
						L_size_2 = mm_initial.compute_L_bit_size();
						PPR_size_2 = mm_initial.compute_PPR_bit_size();
						//outputFile << "L_size_2 " << L_size_2 << " PPR_size_2 " << PPR_size_2 << endl;
					}
				}
			}

			/*clean L*/
			if (1) {
				graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
				binary_read_vector_of_vectors("temp_L.bin", mm_initial.L);
				std::remove("temp_L.bin");

				temp_size_11 = mm_initial.compute_L_bit_size();

				cout << "step 11" << endl;

				auto begin = std::chrono::high_resolution_clock::now();
				clean_L_dynamic(mm_initial.L, mm_initial.PPR, thread_num);
				cleanL_time2 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s

				temp_size_12 = mm_initial.compute_L_bit_size();

				binary_save_vector_of_vectors("temp_L.bin", mm_initial.L);

				cout << "step 12" << endl;
			}

			/*re-ge PPR*/
			if (1) {
				graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
				binary_read_vector_of_vectors("temp_L.bin", mm_initial.L);
				std::remove("temp_L.bin");

				temp_size_13 = mm_initial.compute_L_bit_size();

				cout << "step 13" << endl;

				auto begin = std::chrono::high_resolution_clock::now();
				clean_PPR(instance_graph, mm_initial.L, mm_initial.PPR, thread_num);
				cleanPPR_time2 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s

				temp_size_14 = mm_initial.compute_PPR_bit_size();

				cout << "step 14" << endl;
			}

			/*re-ge*/
			if (1) {

				double time_rege = 0, time_cleanL = 0, time_clean_PPR = 0;

				if (1) {
					graph_hash_of_mixed_weighted g = graph_v_of_v_idealID_to_graph_hash_of_mixed_weighted(instance_graph);
					auto begin = std::chrono::high_resolution_clock::now();
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
					PLL_dynamic(g, instance_graph.size() + 1, thread_num, mm_initial);
					time_rege = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					binary_save_vector_of_vectors("temp_L.bin", mm_initial.L);

					temp_size_15 = mm_initial.compute_L_bit_size();
					temp_size_16 = mm_initial.compute_PPR_bit_size();
				}

				cout << "step 15" << endl;

				if (1) {
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
					binary_read_vector_of_vectors("temp_L.bin", mm_initial.L);
					std::remove("temp_L.bin");

					temp_size_17 = mm_initial.compute_L_bit_size();

					auto begin = std::chrono::high_resolution_clock::now();
					clean_L_dynamic(mm_initial.L, mm_initial.PPR, thread_num);
					time_cleanL = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					L_size_2clean = mm_initial.compute_L_bit_size();
					binary_save_vector_of_vectors("temp_L.bin", mm_initial.L);
				}

				cout << "step 16" << endl;

				if (1) {
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
					binary_read_vector_of_vectors("temp_L.bin", mm_initial.L);
					std::remove("temp_L.bin");

					temp_size_18 = mm_initial.compute_L_bit_size();

					auto begin = std::chrono::high_resolution_clock::now();
					clean_PPR(instance_graph, mm_initial.L, mm_initial.PPR, thread_num);
					time_clean_PPR = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					PPR_size_2clean = mm_initial.compute_PPR_bit_size();
				}

				/*
				just parallel PLL_dynamic generates redundant 2-hop labels, and clean_L_dynamic is required to remove redundant 2-hop labels, and
				clean_PPR is further required to produce the correcponding PPR;

				the strategy of parallelizing PLL with clean_L_dynamic: K. Lakhotia, R. Kannan, Q. Dong, and V. Prasanna, ¡°Planting trees for scalable
				and efficient canonical hub labeling,¡± Proc. VLDB Endow. 13 (2019).
				*/

				rege_time2 = time_rege + time_cleanL + time_clean_PPR; // s

				cout << "step 17" << endl;
			}

			std::remove("temp_L.bin");
			std::remove("temp_PPR.bin");
		}

		outputFile << L_size_0 << "," << PPR_size_0 / L_size_0 << "," << L_size_1 / L_size_0 << "," << PPR_size_1 / L_size_0 << "," <<
			L_size_1clean / L_size_0 << "," << PPR_size_1clean / L_size_0 << "," << cleanL_time1 << "," << cleanPPR_time1 << "," << rege_time1 << "," <<
			L_size_2 / L_size_0 << "," << PPR_size_2 / L_size_0 << "," << L_size_2clean / L_size_0 << "," << PPR_size_2clean / L_size_0 << "," <<
			cleanL_time2 << "," << cleanPPR_time2 << "," << rege_time2 << "," <<
			temp_size_1 / L_size_0 << "," << temp_size_2 / L_size_0 << "," << temp_size_3 / L_size_0 << "," << temp_size_4 / L_size_0 << "," << temp_size_5 / L_size_0 << "," << temp_size_6 / L_size_0 << "," << temp_size_7 / L_size_0 << "," <<
			temp_size_8 / L_size_0 << "," << temp_size_9 / L_size_0 << "," << temp_size_10 / L_size_0 << "," << temp_size_11 / L_size_0 << "," << temp_size_12 / L_size_0 << "," << temp_size_13 / L_size_0 << "," << temp_size_14 / L_size_0 << "," <<
			temp_size_15 / L_size_0 << "," << temp_size_16 / L_size_0 << "," << temp_size_17 / L_size_0 << "," << temp_size_18 / L_size_0 << endl;

		outputFile.close(); // without this, multiple files cannot be successfully created
	}
}

void exp_insert_delete(string data_name, int change_times, double max_Maintain_time, int thread_num) {

	ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);

	string path = "dynamicHL//";

	for (int type = 0; type < 2; type++) {

		graph_v_of_v_idealID instance_graph;
		vector<pair<int, int>> selected_edges;

		string weight_type;
		if (type == 0) {
			weight_type = "Jaccard";
		}
		else if (type == 1) {
			weight_type = "random";
		}
		else {
			weight_type = "unique";
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
			double dummy_ec = 1e4, de_ec = 10;
			int div = 50;

			instance_graph = instance_graph_initial;
			int V = instance_graph.size();
			vector<pair<int, int>>().swap(selected_edges);
			int left_change_times = change_times;
			vector<pair<int, int>> edge_pool;
			for (int i = 0; i < V; i++) {
				for (auto adj : instance_graph[i]) {
					if (i < adj.first && instance_graph[i].size() > 1 && instance_graph[adj.first].size() > 1) {
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
						if (graph_v_of_v_idealID_contain_edge(instance_graph, v1, v2) || v1 == v2 || instance_graph[v1].size() < 5 || instance_graph[v2].size() < 5) {
							continue;
						}
						selected_edges.push_back({ v1, v2 });
						break;
					}
				}
				left_change_times--;
			}

			cout << "step 1" << endl;

			if (data_name == "hyves" && type == 0) {
				div = 10;
			}
			else if (data_name == "skitter" && type == 0) {
				div = 6;
			}

			/*2014+2019*/
			if (1) {
				for (int j = 0; j < change_times / div; j++) {

					cout << "initialize L PPR" << endl;
					instance_graph = instance_graph_initial;
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm = mm_initial;
					initialize_global_values_dynamic(V, thread_num);

					for (int q = 0; q < div; q++) {
						int k = j * div + q;
						cout << "2014+2019 k " << k << endl;

						ThreadPool pool_dynamic(thread_num);
						std::vector<std::future<int>> results_dynamic;

						pair<int, int> selected_edge = selected_edges[k];
						double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
						if (k % 2 == 0) { // increase
							double new_ec = dummy_ec;
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
							try {
								auto begin = std::chrono::high_resolution_clock::now();
								WeightIncrease2019(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic, max_Maintain_time);
								_2019IN_time[k / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
							}
							catch (string s) {
								instance_graph = instance_graph_initial;
								mm = mm_initial;
								initialize_global_values_dynamic(V, thread_num);
								_2019IN_time[k / 2] = INT_MAX;
							}
						}
						else {
							double new_ec = de_ec;
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
							auto begin = std::chrono::high_resolution_clock::now();
							WeightDecrease2014(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
							_2014DE_time[(k - 1) / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
							_20142019_time[(k - 1) / 2] = (_2019IN_time[(k - 1) / 2] + _2014DE_time[(k - 1) / 2]) / 2;
						}
					}
				}
			}

			cout << "step 2" << endl;

			/*2021DE2021IN*/
			if (1) {
				for (int j = 0; j < change_times / div; j++) {

					cout << "initialize L PPR" << endl;
					instance_graph = instance_graph_initial;
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm = mm_initial;
					initialize_global_values_dynamic(V, thread_num);

					for (int q = 0; q < div; q++) {
						int k = j * div + q;
						cout << "2021DE2021IN k " << k << endl;

						ThreadPool pool_dynamic(thread_num);
						std::vector<std::future<int>> results_dynamic;

						pair<int, int> selected_edge = selected_edges[k];
						double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
						if (k % 2 == 0) { // increase
							global_query_times = 0;
							double new_ec = dummy_ec;
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
							try {
								auto begin = std::chrono::high_resolution_clock::now();
								WeightIncrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
								_2021IN_time[k / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
								_2021IN_query_times[k / 2] = global_query_times;
							}
							catch (string s) {
								instance_graph = instance_graph_initial;
								mm = mm_initial;
								initialize_global_values_dynamic(V, thread_num);
								_2021IN_time[k / 2] = INT_MAX;
								_2021IN_query_times[k / 2] = global_query_times;
							}
						}
						else {
							global_query_times = 0;
							double new_ec = de_ec;
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
							try {
								auto begin = std::chrono::high_resolution_clock::now();
								WeightDecrease2021(instance_graph, mm, selected_edge.first, selected_edge.second, selected_edge_weight, new_ec, pool_dynamic, results_dynamic);
								_2021DE_time[(k - 1) / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
								_2021DE_query_times[(k - 1) / 2] = global_query_times;
								_2021DE2021IN_time[(k - 1) / 2] = (_2021IN_time[(k - 1) / 2] + _2021DE_time[(k - 1) / 2]) / 2;
							}
							catch (string s) {
								instance_graph = instance_graph_initial;
								mm = mm_initial;
								initialize_global_values_dynamic(V, thread_num);
								_2021DE_time[(k - 1) / 2] = INT_MAX; // s
								_2021DE_query_times[(k - 1) / 2] = global_query_times;
								_2021DE2021IN_time[(k - 1) / 2] = (_2021IN_time[(k - 1) / 2] + _2021DE_time[(k - 1) / 2]) / 2;
							}
						}
					}
				}
			}

			cout << "step 3" << endl;

			div = 50;

			/*new*/
			if (1) {
				for (int j = 0; j < change_times / div; j++) {

					cout << "initialize L PPR" << endl;
					instance_graph = instance_graph_initial;
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm = mm_initial;
					initialize_global_values_dynamic(V, thread_num);

					for (int q = 0; q < div; q++) {
						int k = j * div + q;
						cout << "new k " << k << endl;

						ThreadPool pool_dynamic(thread_num);
						std::vector<std::future<int>> results_dynamic;

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
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec);
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

	vector<string> data_names = { "condmat", "gnutella", "amazon", "book", "hyves", "skitter" };
	int change_times = 300, thread_num = 80;
	double max_Maintain_time = 100;

	if (1) {
		double weightChange_ratio = 0;
		for (auto data_name : data_names) {
			exp_main(data_name, weightChange_ratio, change_times, max_Maintain_time, thread_num);
			exp_clean(data_name, weightChange_ratio, change_times, max_Maintain_time, thread_num);
			exp_insert_delete(data_name, change_times, max_Maintain_time, thread_num);
		}
	}
}



void exp_case() {

	vector<string> data_names = { "condmat", "gnutella" };
	int change_times = 100, thread_num = 80;

	for (auto data_name : data_names) {

		ofstream outputFile;
		outputFile.precision(6);
		outputFile.setf(ios::fixed);
		outputFile.setf(ios::showpoint);

		string path = "dynamicHL//";

		graph_v_of_v_idealID instance_graph;
		vector<_edge> selected_edges;

		string weight_type = "Jaccard";

		graph_hash_of_mixed_weighted instance_graph_initial_hash = graph_hash_of_mixed_weighted_binary_read(path + data_name + "_" + weight_type + ".bin");
		graph_v_of_v_idealID instance_graph_initial = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID_2(instance_graph_initial_hash, instance_graph_initial_hash.hash_of_vectors.size());
		graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
		binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm_initial.PPR);
		binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm_initial.L);
		string file_name = "exp_case_" + data_name + ".csv";
		cout << file_name << endl;
		outputFile.open(file_name);

		outputFile << "change_v1,change_v2,change_w0,change_w1,change_time,query_v1,query_v2,query_time" << endl;

		double change_time = 0, query_v1 = 0, query_v2 = 0, query_time = 0;

		/*mixed*/
		if (1) {
			double precision = std::pow(10, 3);
			instance_graph = instance_graph_initial;
			int V = instance_graph.size();

			int left_change_times = change_times;
			while (left_change_times) {
				vector<pair<int, int>> edge_pool;
				for (int i = 0; i < V; i++) {
					for (auto adj : instance_graph[i]) {
						if (i < adj.first && adj.second >= 0.1 && adj.second <= 1e4 && instance_graph[i].size() > 1 && instance_graph[adj.first].size() > 1) {
							edge_pool.push_back({ i, adj.first });
						}
					}
				}
				boost::range::random_shuffle(edge_pool);
				for (auto e : edge_pool) {
					pair<int, int> selected_edge = e;
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.first, selected_edge.second);
					if (left_change_times % 2 == 0) { // first increase
						boost::random::uniform_int_distribution<> dist{ static_cast<int>(selected_edge_weight * 2 * precision), static_cast<int>(selected_edge_weight * 10 * precision) };
						double new_ec = dist(boost_random_time_seed) / precision;
						selected_edges.push_back({ selected_edge.first, selected_edge.second, new_ec });
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // increase weight
					}
					else { // then decrease
						boost::random::uniform_int_distribution<> dist{ static_cast<int>(selected_edge_weight * 0.1 * precision), static_cast<int>(selected_edge_weight * 0.5 * precision) };
						double new_ec = dist(boost_random_time_seed) / precision;
						selected_edges.push_back({ selected_edge.first, selected_edge.second, new_ec });
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.first, selected_edge.second, new_ec); // decrease weight
					}
					left_change_times--;
					if (left_change_times == 0) {
						break;
					}
				}
			}

			cout << "step 1" << endl;

			/*new*/
			if (1) {
				instance_graph = instance_graph_initial;
				initialize_global_values_dynamic(V, thread_num);

				for (int k = 0; k < change_times; k++) {

					ThreadPool pool_dynamic(thread_num);
					std::vector<std::future<int>> results_dynamic;

					auto selected_edge = selected_edges[k];
					double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.v1, selected_edge.v2);

					if (k % 2 == 0) { // increase
						auto begin = std::chrono::high_resolution_clock::now();
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec); // increase weight
						WeightIncreaseMaintenance_improv(instance_graph, mm_initial, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic);
						change_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					}
					else {
						auto begin = std::chrono::high_resolution_clock::now();
						graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec);
						WeightDecreaseMaintenance_improv(instance_graph, mm_initial, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic);
						change_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					}

					/*query*/
					if (1) {
						boost::random::uniform_int_distribution<> dist{ static_cast <int>(0), static_cast<int>(V - 1) };
						query_v1 = dist(boost_random_time_seed), query_v2 = dist(boost_random_time_seed);
						auto begin = std::chrono::high_resolution_clock::now();
						double dis = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(mm_initial.L, query_v1, query_v2);
						query_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
					}

					outputFile << selected_edge.v1 << "," << selected_edge.v2 << "," << selected_edge_weight << "," << selected_edge.ec
						<< "," << change_time << "," << query_v1 << "," << query_v2 << "," << query_time << endl;
				}
			}
		}
		outputFile.close(); // without this, multiple files cannot be successfully created
	}
}


void exp_check() {
	string data_name = "skitter";

	string path = "dynamicHL//";

	graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm_initial;
	binary_read_PPR(path + data_name + "_PPR_Jaccard.bin", mm_initial.PPR);

	double PPR_size_1 = mm_initial.compute_PPR_bit_size();

	binary_save_PPR("temp_PPR.bin", mm_initial.PPR);

	binary_read_PPR("temp_PPR.bin", mm_initial.PPR);

	double PPR_size_2 = mm_initial.compute_PPR_bit_size();

	cout << PPR_size_1 << "   " << PPR_size_2 << endl;
}



int main()
{
	cout << "Start running..." << endl;
	auto begin = std::chrono::high_resolution_clock::now();
	/*the two values below are for #include <graph_hash_of_mixed_weighted.h>*/
	graph_hash_of_mixed_weighted_turn_on_value = 1e3;
	graph_hash_of_mixed_weighted_turn_off_value = 1e1;
	srand(time(NULL)); //  seed random number generator

	exp();

	auto end = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	cout << "END    runningtime: " << runningtime << "s" << endl;
}