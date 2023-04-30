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
#include <dgraph_v_of_v/dgraph_v_of_v.h>
#include <dgraph_v_of_v/dgraph_binary_save_read_dgraph.h>
#include <build_in_progress/HL/dgraph/dgraph_two_hop_label.h>
#include <text_mining/parse_substring_between_pairs_of_delimiters.h>
#include <text_mining/parse_substring_between_two_unique_delimiters.h>
#include <text_mining/binary_save_read_vector.h>
#include <text_mining/binary_save_read_vector_of_vectors.h>
//#include <text_mining/list_all_files_in_a_directory.h>

using namespace std;

void update_random_ec(dgraph_v_of_v<two_hop_weight_type>& graph, int vertex_num) {

	boost::random::uniform_int_distribution<> dist{ 1, 100 }; // cannot use small decimal values due to large dummy cost and 1e-7 float precision

	for (int i = 0;i < graph.OUTs.size(); i++)
	{
		int adj_size = graph.OUTs[i].size();
		for (int j = 0;j < adj_size;j++)
		{
			two_hop_weight_type ec = (two_hop_weight_type)dist(boost_random_time_seed);
			graph.add_edge(i,graph.OUTs[i][j].first,ec);
		}
	}

	return;
}

void update_Jacard_ec_element(dgraph_v_of_v<two_hop_weight_type>* input_graph, vector<vector<pair<int, two_hop_weight_type>>>* input_new_ec, int s) {

	auto& graph = *input_graph;
	int adj_size = input_graph->OUTs[s].size();
	
	for (int i = 0;i < adj_size;i++)
	{
		int t = input_graph->OUTs[s][i].first;
		int s_cap_t = 0;

		auto s_pointer = input_graph->OUTs[s].begin();
		auto s_end = input_graph->OUTs[s].end();
		auto t_pointer = input_graph->INs[t].begin();
		auto t_end = input_graph->INs[t].end();
		while (s_pointer != s_end && t_pointer != t_end) {
			if (s_pointer->first == t_pointer->first)
			{
				s_cap_t++;
				s_pointer++;
			}
			else if (s_pointer->first > t_pointer->first)
			{
				t_pointer++;
			}
			else {
				s_pointer++;
			}
		}

		int s_cup_t = input_graph->INs[t].size() + adj_size - s_cap_t;
		two_hop_weight_type ec = 1.0 - (two_hop_weight_type)s_cap_t / (s_cup_t);

		(*input_new_ec)[s].push_back({t,ec});
	}

	return;
}

dgraph_v_of_v<two_hop_weight_type> update_Jacard_ec(dgraph_v_of_v<two_hop_weight_type>& graph, int vertex_num) {

	ThreadPool pool(5);
	std::vector< std::future<int> > results; // return typename: xxx

	dgraph_v_of_v<two_hop_weight_type>* input_graph = &graph;
	dgraph_v_of_v<two_hop_weight_type> output_graph;
	output_graph.INs.resize(graph.INs.size());
	output_graph.OUTs.resize(graph.INs.size());

	vector<vector<pair<int, two_hop_weight_type>>> input_new_ec(graph.INs.size());
	vector<vector<pair<int, two_hop_weight_type>>>* input_new_ec_pointer = &input_new_ec;
	
	for (int i = 0;i < graph.INs.size();i++)
	{
		results.emplace_back(
			pool.enqueue([input_graph, i, input_new_ec_pointer] { // pass const type value j to thread; [] can be empty
				update_Jacard_ec_element(input_graph, input_new_ec_pointer,  i);
				return 1; // return to results; the return type must be the same with results
				})
		);
	}

	for (auto&& result : results)
		result.get(); //all threads finish here

	for (int i = 0;i < input_new_ec.size(); i++)
	{
		int adj_size = input_new_ec[i].size();
		for (int j = 0;j < adj_size;j++)
		{
			output_graph.add_edge(i, input_new_ec[i][j].first, input_new_ec[i][j].second);
		}
	}

	return output_graph;
}

/*vertex_num is non_dummy v num; new_id_name includes both non-dummy and dummy*/
void save_data(dgraph_v_of_v<two_hop_weight_type>& graph, int vertex_num, long long int edge_num, string save_name) {

	string save_file_name;
	ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);

	update_random_ec(graph, vertex_num);
	save_file_name = save_name + "_random.bin";
	dgraph_binary_save_dgraph(save_file_name , graph);

	dgraph_v_of_v<two_hop_weight_type> Jacard_graph = update_Jacard_ec(graph , vertex_num);
	save_file_name = save_name + "_Jacard.bin";
	dgraph_binary_save_dgraph(save_file_name, Jacard_graph);


	save_file_name = save_name + ".txt";
	outputFile.open(save_file_name);
	outputFile << "Edge_terminal_1 Edge_terminal_2 Random_weight Jacard weight (V=" << vertex_num << "  E=" << edge_num << ")" << endl;
	for (int i = 0;i < graph.OUTs.size(); i++)
	{
		int adj_size = graph.OUTs[i].size();
		for (int j=0;j<adj_size;j++)
		{	
			outputFile << i << " " << graph.OUTs[i][j].first << " " << graph.OUTs[i][j].second << " " << Jacard_graph.OUTs[i][j].second << endl;
		}
	}
	
	outputFile.close();
}

//for test
void save_epinions() {
	unordered_map<int, int> old_id_to_new_id;
	dgraph_v_of_v<two_hop_weight_type> instance_graph;

	instance_graph.INs.resize(75879);
	instance_graph.OUTs.resize(75879);
	
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

				if (new_id1 != new_id2)
				{
					instance_graph.add_edge(new_id1,new_id2,1);
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

	save_data(instance_graph, old_id_to_new_id.size(), edge_num, out_folder);

	return;
}

void save_files()
{
	vector<string> files = { "soc-Epinions1" , "citeseer" };
	vector<string> splits = { "\t" ," "};
	vector<int> v_nums = { 75879 , 384413};

	for (int i=0;i<files.size();i++)
	{
		string filename = files[i];
		string split = splits[i];
		int v_num = v_nums[i];

		unordered_map<int, int> old_id_to_new_id;
		dgraph_v_of_v<two_hop_weight_type> instance_graph;

		instance_graph.INs.resize(v_num);
		instance_graph.OUTs.resize(v_num);

		string from_file = "./out." + filename;
		string out_folder = filename;
		string line_content;

		long long int edge_num = 0;

		ifstream myfile(from_file); // open the file
		if (myfile.is_open()) // if the file is opened successfully
		{
			int count = 0;
			while (getline(myfile, line_content)) // read file line by line
			{
				if (count > 1) {
					std::vector<string> Parsed_content = parse_string(line_content, split); // parse line_content

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

					if (new_id1 != new_id2)
					{
						instance_graph.add_edge(new_id1, new_id2, 1);
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

		save_data(instance_graph, old_id_to_new_id.size(), edge_num, out_folder);

	}

	return;
}

