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

	vector<string> data_names = { //"astro", "condmat", "github", "google", "youtube", "skitter" 
	"imdb", "hyves", "amazon" };
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
			clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic, 80);
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
			clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic, 80);
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
			clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic, 80);
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

void exp_element1(string data_name, double weightChange_ratio, int change_times, double max_Maintain_time, int thread_num) {

	ThreadPool pool_dynamic(thread_num);
	std::vector<std::future<int>> results_dynamic;

	ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);

	string path = "dynamicHL//";

	for (int type = 0; type < 2; type++) {

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
		string file_name = "exp_" + data_name + "_T_" + to_string(thread_num) + "_changeRatio_" + to_string((int)(weightChange_ratio * 100)) + "_" + weight_type + ".csv";
		cout << file_name << endl;
		outputFile.open(file_name);

		outputFile << "2014DE_time,2021DE_time,2021DE_query_times,newDE_time,newDE_query_times,DE_ratio,2019IN_time,2021IN_time,2021IN_query_times,newIN_time,newIN_query_times,IN_ratio," <<
			"2014+2019_time,2021DE2021IN_time,2021DEnewIN_time,newDE2021IN_time,newDEnewIN_time," <<
			"L_bit_size_initial(1),PPR_bit_size_initial,L_bit_size_afterM1,PPR_bit_size_afterM1,L_bit_size_afterClean1,PPR_bit_size_afterClean1,cleanL_time1,cleanPPR_time1,rege_time1" << endl;

		int half_change_times = change_times / 2;
		vector<double> _2014DE_time(half_change_times, 0), _2019IN_time(half_change_times, 0), _2021DE_time(half_change_times, 0), _2021DE_query_times(half_change_times, 0), _2021IN_time(half_change_times, 0), _2021IN_query_times(half_change_times, 0),
			_newDE_time(half_change_times, 0), _newDE_query_times(half_change_times, 0), _newIN_time(half_change_times, 0), _newIN_query_times(half_change_times, 0),
			_20142019_time(half_change_times, 0), _2021DE2021IN_time(half_change_times, 0), _2021DEnewIN_time(half_change_times, 0), _newDE2021IN_time(half_change_times, 0), _newDEnewIN_time(half_change_times, 0);
		double L_bit_size_initial = 0, PPR_bit_size_initial = 0, L_bit_size_afterM1 = 0, PPR_bit_size_afterM1 = 0, L_bit_size_afterClean1 = 0, PPR_bit_size_afterClean1 = 0, cleanL_time1 = 0, cleanPPR_time1 = 0, rege_time1 = 0;

		/*mixed*/
		if (1) {
			double precision = std::pow(10, 3);
			int div = 10;

			instance_graph = instance_graph_initial;
			int V = instance_graph.size();
			vector<_edge>().swap(selected_edges);
			int left_change_times = change_times;
			while (left_change_times) {
				vector<pair<int, int>> edge_pool;
				for (int i = 0; i < V; i++) {
					for (auto adj : instance_graph[i]) {
						if (i < adj.first && adj.second >= 0.1 && adj.second <= 1e4 && instance_graph[i].size() > 5 && instance_graph[adj.first].size() > 5) {
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

			/*2014+2019*/
			if (1) {
				for (int j = 0; j < change_times / div; j++) {

					cout << "initialize L PPR" << endl;
					instance_graph = instance_graph_initial;
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
					binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm.PPR);
					binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm.L);
					initialize_global_values_dynamic(V, thread_num);

					for (int q = 0; q < div; q++) {
						int k = j * div + q;
						cout << "2014+2019 k " << k << endl;

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
								mm.clear_labels();
								binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm.PPR);
								binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm.L);
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
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
					binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm.PPR);
					binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm.L);
					initialize_global_values_dynamic(V, thread_num);

					for (int q = 0; q < div; q++) {
						int k = j * div + q;
						cout << "2021DE2021IN k " << k << endl;

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
								mm.clear_labels();
								binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm.PPR);
								binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm.L);
								_2021IN_time[k / 2] = INT_MAX;
								_2021IN_query_times[k / 2] = global_query_times;
							}
						}
						else {
							global_query_times = 0;
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec); // decrease weight
							auto begin = std::chrono::high_resolution_clock::now();
							WeightDecrease2021(instance_graph, mm, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic);
							_2021DE_time[(k - 1) / 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
							_2021DE_query_times[(k - 1) / 2] = global_query_times;
							_2021DE2021IN_time[(k - 1) / 2] = (_2021IN_time[(k - 1) / 2] + _2021DE_time[(k - 1) / 2]) / 2;
						}
					}
				}
			}

			cout << "step 3" << endl;

			/*new*/
			if (1) {
				for (int j = 0; j < change_times / div; j++) {

					cout << "initialize L PPR" << endl;
					instance_graph = instance_graph_initial;
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
					binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm.PPR);
					binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm.L);
					initialize_global_values_dynamic(V, thread_num);

					for (int q = 0; q < div; q++) {
						int k = j * div + q;
						cout << "new k " << k << endl;

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

				if (1) {
					int total_change_times = 1e4;

					/*total_change_times-change_times changes*/
					instance_graph = instance_graph_initial;
					vector<_edge>().swap(selected_edges);
					int left_change_times = total_change_times;
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

					instance_graph = instance_graph_initial;
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
					binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm.PPR);
					binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm.L);
					initialize_global_values_dynamic(V, thread_num);

					for (int k = 0; k < total_change_times; k++) {
						cout << "new large k " << k << endl;
						auto selected_edge = selected_edges[k];
						double selected_edge_weight = graph_v_of_v_idealID_edge_weight(instance_graph, selected_edge.v1, selected_edge.v2);
						if (k % 2 == 0) { // increase
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec); // increase weight
							WeightIncreaseMaintenance_improv(instance_graph, mm, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic);
						}
						else {
							graph_v_of_v_idealID_add_edge(instance_graph, selected_edge.v1, selected_edge.v2, selected_edge.ec);
							WeightDecreaseMaintenance_improv(instance_graph, mm, selected_edge.v1, selected_edge.v2, selected_edge_weight, selected_edge.ec, pool_dynamic, results_dynamic);
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
					graph_hash_of_mixed_weighted g = graph_v_of_v_idealID_to_graph_hash_of_mixed_weighted(instance_graph);
					begin = std::chrono::high_resolution_clock::now();
					PLL_dynamic(g, instance_graph.size() + 1, thread_num, mm);
					clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic, thread_num);
					clean_PPR(instance_graph, mm.L, mm.PPR, pool_dynamic, results_dynamic, thread_num);
					rege_time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
				}

				cout << "step 5" << endl;
			}
		}

		double avg_2014DE_time = 0, avg_2019IN_time = 0, avg_2021DE_time = 0, avg_2021DE_query_times = 0, avg_2021IN_time = 0, avg_2021IN_query_times = 0, avg_DEratio = 0, avg_INratio = 0,
			avg_newDE_time = 0, avg_newDE_query_times = 0, avg_newIN_time = 0, avg_newIN_query_times = 0,
			avg_20142019_time = 0, avg_2021DE2021IN_time = 0, avg_2021DEnewIN_time = 0, avg_newDE2021IN_time = 0, avg_newDEnewIN_time = 0;
		for (int k = 0; k < half_change_times; k++) {
			outputFile << _2014DE_time[k] << "," << _2021DE_time[k] << "," << _2021DE_query_times[k] << "," << _newDE_time[k] << "," << _newDE_query_times[k] << "," << _newDE_time[k] / _2021DE_time[k] << "," <<
				_2019IN_time[k] << "," << _2021IN_time[k] << "," << _2021IN_query_times[k] << "," << _newIN_time[k] << "," << _newIN_query_times[k] << "," << _newIN_time[k] / _2021IN_time[k] << "," <<
				_20142019_time[k] << "," << _2021DE2021IN_time[k] << "," << _2021DEnewIN_time[k] << "," << _newDE2021IN_time[k] << "," << _newDEnewIN_time[k] << "," <<
				L_bit_size_initial << "," << PPR_bit_size_initial / L_bit_size_initial << "," << L_bit_size_afterM1 / L_bit_size_initial << "," << PPR_bit_size_afterM1 / L_bit_size_initial << "," <<
				L_bit_size_afterClean1 / L_bit_size_initial << "," << PPR_bit_size_afterClean1 / L_bit_size_initial << "," << cleanL_time1 << "," << cleanPPR_time1 << "," << rege_time1 << endl;
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
			L_bit_size_initial << "," << PPR_bit_size_initial / L_bit_size_initial << "," << L_bit_size_afterM1 / L_bit_size_initial << "," << PPR_bit_size_afterM1 / L_bit_size_initial << "," <<
			L_bit_size_afterClean1 / L_bit_size_initial << "," << PPR_bit_size_afterClean1 / L_bit_size_initial << "," << cleanL_time1 << "," << cleanPPR_time1 << "," << rege_time1 << endl;

		outputFile.close(); // without this, multiple files cannot be successfully created
	}
}

void exp_element2(string data_name, int change_times, double max_Maintain_time, int thread_num) {

	ThreadPool pool_dynamic(thread_num);
	std::vector<std::future<int>> results_dynamic;

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
			int div = 10;

			instance_graph = instance_graph_initial;
			int V = instance_graph.size();
			vector<pair<int, int>>().swap(selected_edges);
			int left_change_times = change_times;
			vector<pair<int, int>> edge_pool;
			for (int i = 0; i < V; i++) {
				for (auto adj : instance_graph[i]) {
					if (i < adj.first && instance_graph[i].size() > 5 && instance_graph[adj.first].size() > 5) {
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

			/*2014+2019*/
			if (1) {
				for (int j = 0; j < change_times / div; j++) {

					cout << "initialize L PPR" << endl;
					instance_graph = instance_graph_initial;
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
					binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm.PPR);
					binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm.L);
					initialize_global_values_dynamic(V, thread_num);

					for (int q = 0; q < div; q++) {
						int k = j * div + q;
						cout << "2014+2019 k " << k << endl;

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
								mm.clear_labels();
								binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm.PPR);
								binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm.L);
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
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
					binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm.PPR);
					binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm.L);
					initialize_global_values_dynamic(V, thread_num);

					for (int q = 0; q < div; q++) {
						int k = j * div + q;
						cout << "2021DE2021IN k " << k << endl;

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
								mm.clear_labels();
								binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm.PPR);
								binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm.L);
								_2021IN_time[k / 2] = INT_MAX;
								_2021IN_query_times[k / 2] = global_query_times;
							}
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
			}

			cout << "step 3" << endl;

			/*new*/
			if (1) {
				for (int j = 0; j < change_times / div; j++) {

					cout << "initialize L PPR" << endl;
					instance_graph = instance_graph_initial;
					graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
					binary_read_PPR(path + data_name + "_PPR_" + weight_type + ".bin", mm.PPR);
					binary_read_vector_of_vectors(path + data_name + "_L_" + weight_type + ".bin", mm.L);
					initialize_global_values_dynamic(V, thread_num);

					for (int q = 0; q < div; q++) {
						int k = j * div + q;
						cout << "new k " << k << endl;

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

	vector<string> data_names = { //"astro", "condmat", "github", 
		"google", "youtube", "hyves", "skitter" };
	int change_times = 300, thread_num = 80;
	double max_Maintain_time = 100;

	if (1) {
		double weightChange_ratio = 0;
		for (auto data_name : data_names) {
			exp_element1(data_name, weightChange_ratio, change_times, max_Maintain_time, thread_num);
			exp_element2(data_name, change_times, max_Maintain_time, thread_num);
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
	srand(time(NULL)); //  seed random number generator

	exp();

	auto end = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	cout << "END    runningtime: " << runningtime << "s" << endl;
}