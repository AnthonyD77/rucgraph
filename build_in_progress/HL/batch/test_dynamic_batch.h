#pragma once
#include <build_in_progress/HL/dynamic/PLL_dynamic.h>
#include <build_in_progress/HL/dynamic/WeightIncreaseMaintenance_improv_multiThread.h>
#include <build_in_progress/HL/dynamic/WeightDecreaseMaintenance_improv_multiThread.h>
#include <build_in_progress/HL/dynamic/WeightIncrease2021_multiThread.h>
#include <build_in_progress/HL/dynamic/WeightDecrease2021_multiThread.h>
#include <build_in_progress/HL/sort_v/graph_hash_of_mixed_weighted_update_vertexIDs_by_degrees.h>
#include <graph_hash_of_mixed_weighted/two_graphs_operations/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID_2.h>
#include <graph_hash_of_mixed_weighted/random_graph/graph_hash_of_mixed_weighted_generate_random_graph.h>
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_read_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_save_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/common_algorithms/graph_hash_of_mixed_weighted_shortest_paths.h>
#include <build_in_progress/HL/dynamic/clean_labels.h>

#include <build_in_progress/HL/batch/WeightChangeMaintenance_batch_base.h>
#include <build_in_progress/HL/batch/WeightDecrease_batch.h>
#include <build_in_progress/HL/batch/WeightIncrease_batch.h>

void check_correctness_dynamic_batch(graph_hash_of_mixed_weighted_two_hop_case_info_v1& case_info,
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


void graph_changes_batch(graph_v_of_v_idealID& ideal_g, graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	int V, int weightIncrease_time, int weightDecrease_time, double weightChange_ratio, std::vector<graph_edge_change>& edge_changes) {

    std::vector<graph_edge_change>().swap(edge_changes);
	auto temp_g = ideal_g;
	
	//----------debug----------
	// weightIncrease_time = 0;
	// weightDecrease_time = 1;
	//-------------------------

	while (weightIncrease_time + weightDecrease_time) {
		int v2_id;
		/*randomly select an edge*/
		pair<int, int> selected_edge;
		double selected_edge_weight;
		while (1) {
			boost::random::uniform_int_distribution<> dist_v1{ static_cast<int>(0), static_cast<int>(V - 1) };
			int v1 = dist_v1(boost_random_time_seed);
			if (temp_g[v1].size() > 0) {
				boost::random::uniform_int_distribution<> dist_v2{ static_cast<int>(0), static_cast<int>(temp_g[v1].size() - 1) };
				int rand = dist_v2(boost_random_time_seed);

				#ifdef FIX_GRAPH
					// v1 = 0;
					// rand = 1;
				#endif 

				int v2 = temp_g[v1][rand].first;
				selected_edge = { v1,v2 };
				selected_edge_weight = temp_g[v1][rand].second;
				v2_id = rand;
				break;
			}
		}

		/*change weight*/
		double new_ec;
		if (weightIncrease_time > 0) { //weightIncrease_time >= weightDecrease_time
			weightIncrease_time--;
			new_ec = selected_edge_weight * (1 + weightChange_ratio);
			if (new_ec > 1e6) {
				weightIncrease_time++;
				continue;
			}
            edge_changes.push_back(graph_edge_change(selected_edge.first, selected_edge.second, v2_id, selected_edge_weight, new_ec));
		}
		else {
			weightDecrease_time--;
			new_ec = selected_edge_weight * (1 - weightChange_ratio);
			if (new_ec < 1e-2) {
				weightDecrease_time++;
				continue;
			}
            edge_changes.push_back(graph_edge_change(selected_edge.first, selected_edge.second, v2_id, selected_edge_weight, new_ec));
		}
		graph_v_of_v_idealID_add_edge(temp_g, selected_edge.first, selected_edge.second, new_ec);
	}
}

void graph_label_maintenance_batch(graph_v_of_v_idealID& ideal_g, graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	std::vector<graph_edge_change>& edge_changes, int thread_num, double& avg_maintain_time) {

	std::vector<graph_edge_change> edge_increases, edge_decreases;
	
	for (auto it : edge_changes) {
		if (it.new_weight > it.old_weight) {
			edge_increases.push_back(it);
		}
		else {
			edge_decreases.push_back(it);
		}
	}

	#ifdef CHECK_PROCESS
		cerr << "Starting" << endl;
	#endif
	
	WeightIncrease_batch_original(ideal_g, instance_graph, mm, edge_increases, thread_num, avg_maintain_time);
	// WeightIncrease_batch(ideal_g, instance_graph, mm, edge_increases, avg_maintain_time);
	#ifdef CHECK_PROCESS
		cerr << "Increase maintenance compeleted" << endl;
	#endif	
	// WeightDecrease_batch_original(ideal_g, instance_graph, mm, edge_decreases, thread_num, avg_maintain_time);
	WeightDecrease_batch(ideal_g, instance_graph, mm, edge_decreases, avg_maintain_time);
	#ifdef CHECK_PROCESS
		cerr << "Decrease maintenance compeleted" << endl;
	#endif	
}

// void graph_label_maintenance_batch1(graph_v_of_v_idealID& ideal_g, graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
// 	std::vector<graph_edge_change>& edge_changes, int thread_num, double& avg_maintain_time) {

// 	ThreadPool pool_dynamic(thread_num);
// 	std::vector<std::future<int>> results_dynamic;

//     for (auto it : edge_changes) {
// 		// it.old_weight = ideal_g[it.source][it.terminal_id].second;
//         if (it.new_weight > it.old_weight) {
//             graph_hash_of_mixed_weighted_add_edge(instance_graph, it.source, it.terminal, it.new_weight); // increase weight
// 			graph_v_of_v_idealID_add_edge(ideal_g, it.source, it.terminal, it.new_weight);

// 			auto begin = std::chrono::high_resolution_clock::now();

// 			/*maintain labels*/
// 			//WeightIncrease2021(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight, pool_dynamic, results_dynamic);
// 			WeightIncreaseMaintenance_improv(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight, pool_dynamic, results_dynamic);
// 			//WeightIncrease2019(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight, pool_dynamic, results_dynamic, 1e1);
// 			auto end = std::chrono::high_resolution_clock::now();
// 			avg_maintain_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
//         }
//         else {
//             graph_hash_of_mixed_weighted_add_edge(instance_graph, it.source, it.terminal, it.new_weight); // decrease weight
// 			graph_v_of_v_idealID_add_edge(ideal_g, it.source, it.terminal, it.new_weight);

// 			auto begin = std::chrono::high_resolution_clock::now();

// 			/*maintain labels*/
// 			WeightDecreaseMaintenance_improv(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight, pool_dynamic, results_dynamic);
// 			//WeightDecrease2021(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight, pool_dynamic, results_dynamic);
// 			//WeightDecrease2014(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight, pool_dynamic, results_dynamic);
// 			auto end = std::chrono::high_resolution_clock::now();
// 			avg_maintain_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
//         }
//     }
// }

void test_dynamic_batch() {

	srand(time(NULL)); //  seed random number generator

	/*parameters*/
	int iteration_graph_times = 1e4, iteration_source_times = 30, iteration_terminal_times = 30;
	int iteration_change_round = 30;
	int V = 100, E = 500, precision = 1, thread_num = 10;
	double ec_min = 1, ec_max = 10;

	int weightIncrease_time = 10, weightDecrease_time = 10;
	double weightChange_ratio = 0.2;

	double avg_index_time = 0, avg_index_size_per_v = 0, avg_maintain_time = 0;

    std::vector<graph_edge_change> edge_changes;

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

		#ifdef FIX_GRAPH
			generate_new_graph = 0;
		#endif

		graph_hash_of_mixed_weighted instance_graph;
		if (generate_new_graph == 1) {
			instance_graph = graph_hash_of_mixed_weighted_generate_random_graph(V, E, 0, 0, ec_min, ec_max, precision, boost_random_time_seed);
			instance_graph = graph_hash_of_mixed_weighted_update_vertexIDs_by_degrees_large_to_small(instance_graph); // sort vertices
			graph_hash_of_mixed_weighted_save_graph_with_weight("simple_iterative_tests.txt", instance_graph, 0);
		}
		else {
			double lambda;
			graph_hash_of_mixed_weighted_read_graph_with_weight("simple_iterative_tests.txt", instance_graph, lambda);
		}
		graph_v_of_v_idealID ideal_g = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID_2(instance_graph, V);

		auto begin = std::chrono::high_resolution_clock::now();
		try {
			PLL_dynamic(instance_graph, V + 1, thread_num, mm);	// ?
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
		clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic, thread_num);
		initialize_global_values_dynamic(V, thread_num);
		initialize_global_values_dynamic_batch(V, thread_num);

		for (int change_round = 0; change_round < iteration_change_round; change_round ++) {
			#ifdef CHECK_PROCESS
				cerr << "round " << change_round << endl;
			#endif
			graph_changes_batch(ideal_g, instance_graph, mm, V, weightIncrease_time, weightDecrease_time, weightChange_ratio, edge_changes);
			#ifdef CHECK_PROCESS
				cerr << "graph changes generated " << endl;
			#endif
			graph_label_maintenance_batch(ideal_g, instance_graph, mm, edge_changes, thread_num, avg_maintain_time);
			#ifdef CHECK_PROCESS
				cerr << "graph maintenance completed " << endl;
			#endif
		}
		// graph_change_and_label_maintenance(ideal_g, instance_graph, mm, V, weightIncrease_time, weightDecrease_time, weightChange_ratio, thread_num, avg_maintain_time);
		
        check_correctness_dynamic_batch(mm, instance_graph, iteration_source_times, iteration_terminal_times);

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
