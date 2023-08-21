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
		outputFile << "time_generate_labels=" << mm.time_generate_labels << endl;
		outputFile.close();

		g = graph_hash_of_mixed_weighted_binary_read(path + s + "_Jacard.bin");
		PLL_dynamic(g, g.hash_of_vectors.size(), thread_num, mm);
		clean_L_dynamic(mm.L, mm.PPR, pool_dynamic, results_dynamic);
		binary_save_PPR(path + s + "_PPR_Jacard.bin", mm.PPR);
		binary_save_vector_of_vectors(path + s + "_L_Jacard.bin", mm.L);
		outputFile.open(path + s + "_L_Jacard_generation.txt");
		outputFile << "time_generate_labels=" << mm.time_generate_labels << endl;
		outputFile.close();
	}
}



void exp() {



}


















int main()
{
	cout << "Start running..." << endl;
	auto begin = std::chrono::high_resolution_clock::now();
	/*the two values below are for #include <graph_hash_of_mixed_weighted.h>*/
	graph_hash_of_mixed_weighted_turn_on_value = 1e3;
	graph_hash_of_mixed_weighted_turn_off_value = 1e1;

	generate_L_PPR();

	auto end = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	cout << "END    runningtime: " << runningtime << "s" << endl;
}