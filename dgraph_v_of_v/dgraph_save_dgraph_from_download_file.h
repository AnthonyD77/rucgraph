#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <string>
#include <list>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <typeinfo>
#include <unordered_set>
#include <unordered_map>
#include <climits>
#include <math.h>
#include <thread>
#include <chrono>
#include <shared_mutex>
using namespace std;


/*header files in the Boost library: https://www.boost.org/ */
#include <boost/random.hpp>

//boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };
std::shared_mutex glock;

#include <tool_functions/ThreadPool.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_binary_save_read.h>
#include <build_in_progress/HL/dgraph/dgraph_two_hop_label.h>
#include <text_mining/parse_substring_between_pairs_of_delimiters.h>
#include <text_mining/parse_substring_between_two_unique_delimiters.h>
#include <text_mining/binary_save_read_vector.h>
#include <text_mining/binary_save_read_vector_of_vectors.h>
//#include <text_mining/list_all_files_in_a_directory.h>

using namespace std;

//t use to choose type(double or float)
template <typename weight_type>
void dgraph_hash_of_mixed_weighted_binary_save(graph_hash_of_mixed_weighted& save_graph, std::string save_path, weight_type t) {

	std::vector<std::vector<std::pair<int, weight_type>>> save_vectors;

	for (auto it = save_graph.hash_of_vectors.begin(); it != save_graph.hash_of_vectors.end(); it++) {
		int v = it->first;
		weight_type w_v = it->second.vertex_weight;

		std::vector<std::pair<int, weight_type>> save_vector;

		//save_vector = it->second.adj_vertices;
		for (auto& p : it->second.adj_vertices)
		{
			save_vector.push_back(p);
		}

		save_vector.push_back({ v, w_v });
		save_vector.push_back({ 0, NULL });

		save_vectors.push_back(save_vector);
	}

	binary_save_vector_of_vectors(save_path, save_vectors);
}

void dgraph_hash_of_mixed_weighted_add_edge(graph_hash_of_mixed_weighted& input_graph, int e1, int e2, double ec) {

	/*this function adds a weighted edge via binary search,
	and may add e1 and e2 into input_graph if they are new vertices;
	time complexity O(graph_hash_of_mixed_weighted_turn_on_value), which can be considered as O(1);
	this function can update edge weights*/

	/*add e1 and its adj*/
	auto search = input_graph.hash_of_vectors.find(e1);
	if (search == input_graph.hash_of_vectors.end()) { // e1 is a new vertex
		graph_hash_of_mixed_weighted_add_vertex(input_graph, e1, 0); // add e1; initial weight is 0
		search = input_graph.hash_of_vectors.find(e1);
	}

	graph_hash_of_mixed_weighted_binary_operations_insert(search->second.adj_vertices, e2, ec);

	/*add e2*/
	search = input_graph.hash_of_vectors.find(e2);
	if (search == input_graph.hash_of_vectors.end()) { // e2 is a new vertex
		graph_hash_of_mixed_weighted_add_vertex(input_graph, e2, 0); // add e2; initial weight is 0
	}
}

void update_random_ec(graph_hash_of_mixed_weighted& graph, int vertex_num) {

	boost::random::uniform_int_distribution<> dist{ 1, 100 }; // cannot use small decimal values due to large dummy cost and 1e-7 float precision

	vector<pair<pair<int, int>, two_hop_weight_type>> new_ec;

	for (auto it1 = graph.hash_of_vectors.begin(); it1 != graph.hash_of_vectors.end(); it1++) {
		int v1 = it1->first;
		if (v1 >= vertex_num) {
			continue;
		}
		std::vector<int> adj_v1 = graph.adj_v(v1);
		for (int i = 0; i < adj_v1.size(); i++) {
			int v2 = adj_v1[i];

			if (v2 < vertex_num) {
				two_hop_weight_type ec = (two_hop_weight_type)dist(boost_random_time_seed);
				new_ec.push_back({ {v1, v2} ,ec });
			}

		}
	}

	int size = new_ec.size();
	for (int i = 0; i < size; i++) {
		dgraph_hash_of_mixed_weighted_add_edge(graph, new_ec[i].first.first, new_ec[i].first.second, new_ec[i].second);
	}

}

void update_Jacard_ec_element(graph_hash_of_mixed_weighted* input_graph, graph_hash_of_mixed_weighted* reverse_graph, int vertex_num, int i,
	vector<pair<pair<int, int>, double>>* input_new_ec, vector<vector<pair<int, double>>>* input_adj_lists) {

	auto& graph = *input_graph;
	auto& new_ec = *input_new_ec;
	auto& adj_lists = *input_adj_lists;

	std::vector<int> adj_v1 = graph.adj_v(i);
	for (int x = 0; x < adj_v1.size(); x++) {
		int j = adj_v1[x];

		if (j < vertex_num) {

			int V_i_cap_V_j = 0;
			int V_i_non_dummy_num = adj_lists[i].size();
			int V_j_non_dummy_num = reverse_graph->hash_of_vectors[j].adj_vertices.size();

			/*merge join style compare*/
			auto vector_i_check_pointer = adj_lists[i].begin();
			auto vector_j_check_pointer = reverse_graph->hash_of_vectors[j].adj_vertices.begin();
			auto pointer_vector_i_end = adj_lists[i].end();
			auto pointer_vector_j_end = reverse_graph->hash_of_vectors[j].adj_vertices.end();
			while (vector_i_check_pointer != pointer_vector_i_end && vector_j_check_pointer != pointer_vector_j_end) {
				if (vector_i_check_pointer->first >= vertex_num || vector_j_check_pointer->first >= vertex_num) {
					break;
				}
				if (vector_i_check_pointer->first == vector_j_check_pointer->first) {
					V_i_cap_V_j++;
					vector_i_check_pointer++;
				}
				else if (vector_i_check_pointer->first > vector_j_check_pointer->first) {
					vector_j_check_pointer++;
				}
				else {
					vector_i_check_pointer++;
				}
			}
			int V_i_cup_V_j = V_i_non_dummy_num + V_j_non_dummy_num - V_i_cap_V_j;
			double ec = 1 - (double)V_i_cap_V_j / V_i_cup_V_j;

			glock.lock();
			new_ec.push_back({ {i, j} ,ec });
			glock.unlock();
		}

	}

}

graph_hash_of_mixed_weighted update_Jacard_ec(graph_hash_of_mixed_weighted graph, graph_hash_of_mixed_weighted& reverse_graph, int vertex_num) {

	ThreadPool pool(5);
	std::vector< std::future<int> > results; // return typename: xxx

	vector<pair<pair<int, int>, double>> new_ec;

	vector<vector<pair<int, double>>> adj_lists(graph.hash_of_vectors.size() + 1); // sorted lists
	for (auto it1 = graph.hash_of_vectors.begin(); it1 != graph.hash_of_vectors.end(); it1++) {
		adj_lists[it1->first] = graph.adj_v_and_ec(it1->first);
		/*if graph.hash_of_vectors.size() < number of vertices (some vertices are not added) then bugs here
		, this is due to the fact that some vertices are isolated, and could not be added by adding edges,*/
	}



	graph_hash_of_mixed_weighted* input_graph = &graph;
	graph_hash_of_mixed_weighted* reverse_input_graph = &reverse_graph;
	vector<pair<pair<int, int>, double>>* input_new_ec = &new_ec;
	vector<vector<pair<int, double>>>* input_adj_lists = &adj_lists;
	for (auto it1 = graph.hash_of_vectors.begin(); it1 != graph.hash_of_vectors.end(); it1++) {
		int i = it1->first;
		if (i >= vertex_num) {
			continue;
		}
		results.emplace_back(
			pool.enqueue([input_graph, vertex_num, i, input_new_ec, input_adj_lists, reverse_input_graph] { // pass const type value j to thread; [] can be empty
				update_Jacard_ec_element(input_graph, reverse_input_graph, vertex_num, i, input_new_ec, input_adj_lists);
				return 1; // return to results; the return type must be the same with results
				})
		);
	}
	for (auto&& result : results)
		result.get(); //all threads finish here

	int size = new_ec.size();
	for (int i = 0; i < size; i++) {
		dgraph_hash_of_mixed_weighted_add_edge(graph, new_ec[i].first.first, new_ec[i].first.second, new_ec[i].second);
	}

	return graph;
}

/*vertex_num is non_dummy v num; new_id_name includes both non-dummy and dummy*/
void save_data(graph_hash_of_mixed_weighted& graph, graph_hash_of_mixed_weighted& reverse_graph, int vertex_num, long long int edge_num, string save_name) {

	string save_file_name;
	ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);

	update_random_ec(graph, vertex_num);
	save_file_name = save_name + "_random.bin";
	dgraph_hash_of_mixed_weighted_binary_save(graph, save_file_name, (two_hop_weight_type)1);

	graph_hash_of_mixed_weighted Jacard_graph = update_Jacard_ec(graph, reverse_graph, vertex_num);
	save_file_name = save_name + "_Jacard.bin";
	dgraph_hash_of_mixed_weighted_binary_save(Jacard_graph, save_file_name, (two_hop_weight_type)1);


	save_file_name = save_name + ".txt";
	outputFile.open(save_file_name);
	outputFile << "Edge_terminal_1 Edge_terminal_2 Random_weight Jacard weight (V=" << vertex_num << "  E=" << edge_num << ")" << endl;
	for (auto it1 = graph.hash_of_vectors.begin(); it1 != graph.hash_of_vectors.end(); it1++) {
		int v1 = it1->first;
		std::vector<int> adj_v1 = graph.adj_v(v1);
		for (int i = 0; i < adj_v1.size(); i++) {
			int v2 = adj_v1[i];

			outputFile << v1 << " " << v2 << " " << (two_hop_weight_type)graph_hash_of_mixed_weighted_edge_weight(graph, v1, v2) << " " << (two_hop_weight_type)graph_hash_of_mixed_weighted_edge_weight(Jacard_graph, v1, v2) << endl;

		}
	}
	outputFile.close();
}

//for test
void save_epinions() {
	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "./out.soc-Epinions1";
	string out_folder = "./soc-Epinions1";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, "\t"); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_cityseer() {
	/*15min; very slow, due to updating Jacard dis*/

	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "../data/citeseer/out.citeseer";
	string out_folder = "../result/citeseer/citeseer";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, " "); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}
				else {
					graph_hash_of_mixed_weighted_add_vertex(graph_with_1_weight, new_id1, 0);
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	//cout << demmy_node_num << endl;
	int before_num = graph_with_1_weight.hash_of_vectors.size();
	if (old_id_to_new_id.size() != before_num)
	{
		for (int i = before_num;i < old_id_to_new_id.size();i++)
		{
			graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
			reverse_graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
		}
	}

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_amazon0601() {
	/*15min; very slow, due to updating Jacard dis*/

	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "../data/amazon0601/out.amazon0601";
	string out_folder = "../result/amazon0601/amazon0601";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, "	"); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}
				else {
					graph_hash_of_mixed_weighted_add_vertex(graph_with_1_weight, new_id1, 0);
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	//cout << demmy_node_num << endl;
	int before_num = graph_with_1_weight.hash_of_vectors.size();
	if (old_id_to_new_id.size() != before_num)
	{
		for (int i = before_num;i < old_id_to_new_id.size();i++)
		{
			graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
			reverse_graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
		}
	}

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_digg_friends() {
	/*15min; very slow, due to updating Jacard dis*/

	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "../data/digg-friends/out.digg-friends";
	string out_folder = "../result/digg-friends/digg-friends";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, " "); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}
				else {
					graph_hash_of_mixed_weighted_add_vertex(graph_with_1_weight, new_id1, 0);
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	//cout << demmy_node_num << endl;
	int before_num = graph_with_1_weight.hash_of_vectors.size();
	if (old_id_to_new_id.size() != before_num)
	{
		for (int i = before_num;i < old_id_to_new_id.size();i++)
		{
			graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
			reverse_graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
		}
	}

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_flickr_growth() {
	/*15min; very slow, due to updating Jacard dis*/

	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "../data/flickr-growth/out.flickr-growth";
	string out_folder = "../result/flickr-growth/flickr-growth";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, " "); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}
				else {
					graph_hash_of_mixed_weighted_add_vertex(graph_with_1_weight, new_id1, 0);
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	//cout << demmy_node_num << endl;
	int before_num = graph_with_1_weight.hash_of_vectors.size();
	if (old_id_to_new_id.size() != before_num)
	{
		for (int i = before_num;i < old_id_to_new_id.size();i++)
		{
			graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
			reverse_graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
		}
	}

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_higgs_twitter_social() {
	/*15min; very slow, due to updating Jacard dis*/

	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "../data/higgs-twitter-social/out.higgs-twitter-social";
	string out_folder = "../result/higgs-twitter-social/higgs-twitter-social";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, "	"); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}
				else {
					graph_hash_of_mixed_weighted_add_vertex(graph_with_1_weight, new_id1, 0);
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	//cout << demmy_node_num << endl;
	int before_num = graph_with_1_weight.hash_of_vectors.size();
	if (old_id_to_new_id.size() != before_num)
	{
		for (int i = before_num;i < old_id_to_new_id.size();i++)
		{
			graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
			reverse_graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
		}
	}

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_lasagne_yahoo() {
	/*15min; very slow, due to updating Jacard dis*/

	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "../data/lasagne-yahoo/out.lasagne-yahoo";
	string out_folder = "../result/lasagne-yahoo/lasagne-yahoo";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, " "); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}
				else {
					graph_hash_of_mixed_weighted_add_vertex(graph_with_1_weight, new_id1, 0);
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	//cout << demmy_node_num << endl;
	int before_num = graph_with_1_weight.hash_of_vectors.size();
	if (old_id_to_new_id.size() != before_num)
	{
		for (int i = before_num;i < old_id_to_new_id.size();i++)
		{
			graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
			reverse_graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
		}
	}

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_libimseti() {
	/*15min; very slow, due to updating Jacard dis*/

	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "../data/libimseti/out.libimseti";
	string out_folder = "../result/libimseti/libimseti";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, "	"); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}
				else {
					graph_hash_of_mixed_weighted_add_vertex(graph_with_1_weight, new_id1, 0);
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	//cout << demmy_node_num << endl;
	int before_num = graph_with_1_weight.hash_of_vectors.size();
	if (old_id_to_new_id.size() != before_num)
	{
		for (int i = before_num;i < old_id_to_new_id.size();i++)
		{
			graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
			reverse_graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
		}
	}

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_patentcite() {
	/*15min; very slow, due to updating Jacard dis*/

	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "../data/patentcite/out.patentcite";
	string out_folder = "../result/patentcite/patentcite";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, " "); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}
				else {
					graph_hash_of_mixed_weighted_add_vertex(graph_with_1_weight, new_id1, 0);
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	//cout << demmy_node_num << endl;
	int before_num = graph_with_1_weight.hash_of_vectors.size();
	if (old_id_to_new_id.size() != before_num)
	{
		for (int i = before_num;i < old_id_to_new_id.size();i++)
		{
			graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
			reverse_graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
		}
	}

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_prosper_loans() {
	/*15min; very slow, due to updating Jacard dis*/

	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "../data/prosper-loans/out.prosper-loans";
	string out_folder = "../result/prosper-loans/prosper-loans";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, " "); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}
				else {
					graph_hash_of_mixed_weighted_add_vertex(graph_with_1_weight, new_id1, 0);
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	//cout << demmy_node_num << endl;
	int before_num = graph_with_1_weight.hash_of_vectors.size();
	if (old_id_to_new_id.size() != before_num)
	{
		for (int i = before_num;i < old_id_to_new_id.size();i++)
		{
			graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
			reverse_graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
		}
	}

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_soc_LiveJournal1() {
	/*15min; very slow, due to updating Jacard dis*/

	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "../data/soc-LiveJournal1/out.soc-LiveJournal1-loans";
	string out_folder = "../result/soc-LiveJournal1/soc-LiveJournal1";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, "	"); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}
				else {
					graph_hash_of_mixed_weighted_add_vertex(graph_with_1_weight, new_id1, 0);
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	//cout << demmy_node_num << endl;
	int before_num = graph_with_1_weight.hash_of_vectors.size();
	if (old_id_to_new_id.size() != before_num)
	{
		for (int i = before_num;i < old_id_to_new_id.size();i++)
		{
			graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
			reverse_graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
		}
	}

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_soc_pokec_relationships() {
	/*15min; very slow, due to updating Jacard dis*/

	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "../data/soc-pokec-relationships/out.soc-pokec-relationships";
	string out_folder = "../result/soc-pokec-relationships/soc-pokec-relationships";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, "	"); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}
				else {
					graph_hash_of_mixed_weighted_add_vertex(graph_with_1_weight, new_id1, 0);
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	//cout << demmy_node_num << endl;
	int before_num = graph_with_1_weight.hash_of_vectors.size();
	if (old_id_to_new_id.size() != before_num)
	{
		for (int i = before_num;i < old_id_to_new_id.size();i++)
		{
			graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
			reverse_graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
		}
	}

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_web_Stanford() {
	/*15min; very slow, due to updating Jacard dis*/

	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "../data/web-Stanford/out.web-Stanford";
	string out_folder = "../result/web-Stanford/web-Stanford";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, "	"); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}
				else {
					graph_hash_of_mixed_weighted_add_vertex(graph_with_1_weight, new_id1, 0);
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	//cout << demmy_node_num << endl;
	int before_num = graph_with_1_weight.hash_of_vectors.size();
	if (old_id_to_new_id.size() != before_num)
	{
		for (int i = before_num;i < old_id_to_new_id.size();i++)
		{
			graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
			reverse_graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
		}
	}

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_youtube_links() {
	/*15min; very slow, due to updating Jacard dis*/

	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "../data/youtube-links/out.youtube-links";
	string out_folder = "../result/youtube-links/youtube-links";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, " "); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}
				else {
					graph_hash_of_mixed_weighted_add_vertex(graph_with_1_weight, new_id1, 0);
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	//cout << demmy_node_num << endl;
	int before_num = graph_with_1_weight.hash_of_vectors.size();
	if (old_id_to_new_id.size() != before_num)
	{
		for (int i = before_num;i < old_id_to_new_id.size();i++)
		{
			graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
			reverse_graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
		}
	}

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_zhishi_baidu_internallink() {
	/*15min; very slow, due to updating Jacard dis*/

	unordered_map<int, int> old_id_to_new_id;
	graph_hash_of_mixed_weighted graph_with_1_weight;
	graph_hash_of_mixed_weighted reverse_graph_with_1_weight;

	string from_file = "../data/zhishi-baidu-internallink/out.zhishi-baidu-internallink";
	string out_folder = "../result/zhishi-baidu-internallink/zhishi-baidu-internallink";
	string line_content;

	long long int edge_num = 0;

	ifstream myfile(from_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			if (count > 1) {
				std::vector<string> Parsed_content = parse_string(line_content, " "); // parse line_content

				int old_id1 = stoi(Parsed_content[0]);
				int new_id1 = 0, new_id2 = 0;
				if (old_id_to_new_id.count(old_id1) > 0) {
					new_id1 = old_id_to_new_id[old_id1];
				}
				else {
					new_id1 = old_id_to_new_id.size();
					old_id_to_new_id[old_id1] = new_id1;
				}

				int old_id2 = stoi(Parsed_content[1]);
				if (old_id_to_new_id.count(old_id2) > 0) {
					new_id2 = old_id_to_new_id[old_id2];
				}
				else {
					new_id2 = old_id_to_new_id.size();
					old_id_to_new_id[old_id2] = new_id2;
				}

				//new_id1 = old_id1;
				//new_id2 = old_id2;

				if (new_id1 != new_id2)
				{
					dgraph_hash_of_mixed_weighted_add_edge(graph_with_1_weight, new_id1, new_id2, 1);
					dgraph_hash_of_mixed_weighted_add_edge(reverse_graph_with_1_weight, new_id2, new_id1, 1);
					edge_num++;
				}
				else {
					graph_hash_of_mixed_weighted_add_vertex(graph_with_1_weight, new_id1, 0);
				}

			}
			count++;
		}
		myfile.close(); //close the file
	}
	else
	{
		std::cout << "Unable to open file " << from_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
	cout << "1" << endl;

	//cout << demmy_node_num << endl;
	int before_num = graph_with_1_weight.hash_of_vectors.size();
	if (old_id_to_new_id.size() != before_num)
	{
		for (int i = before_num;i < old_id_to_new_id.size();i++)
		{
			graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
			reverse_graph_with_1_weight.hash_of_vectors[i].vertex_weight = 0;
		}
	}

	save_data(graph_with_1_weight, reverse_graph_with_1_weight, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}