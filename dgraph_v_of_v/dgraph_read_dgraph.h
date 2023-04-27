#pragma once
#include <text_mining/parse_string.h>
#include <dgraph_v_of_v/dgraph_v_of_v.h>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>

/**
 * @brief 
 * 
 * @param instance_name : file name of storage dgraph
 * @param input_graph   : dgraph to be read in
 */

template <typename weight_type>
dgraph_v_of_v<weight_type> dgraph_hash_of_mixed_weighted_binary_read(std::string save_path, dgraph_v_of_v<weight_type>& save_graph) {

    std::vector<std::vector<std::pair<int, weight_type>>> save_vectors;
    binary_read_vector_of_vectors(save_path, save_vectors);

    //dgraph_v_of_v<weight_type> save_graph;

    save_graph.INs.resize(save_vectors.size());
    save_graph.OUTs.resize(save_vectors.size());

    for (auto it = save_vectors.begin(); it != save_vectors.end(); it++) {

        int vector_size = (*it).size();

        int from_node = (*it)[vector_size - 2].first;
        //std::vector<std::pair<int, double>> adj_list((*it).begin(), (*it).begin() + vector_size - 2); // https://www.techiedelight.com/get-slice-sub-vector-from-vector-cpp/

        for (auto it2 = (*it).begin(); it2 != (*it).begin() + vector_size - 2; it2++)
        {
            int end_node = it2->first;
            weight_type w = it2->second;
            save_graph.add_edge(from_node, end_node, w);
        }
    }
    return save_graph;
}

template <typename weight_type>
void dgraph_read_dgraph(std::string instance_name, dgraph_v_of_v<weight_type> &input_graph)
{
    input_graph.clear();
    std::string line_content;
    std::ifstream myfile(instance_name);
    if (myfile.is_open())
    {
        while (getline(myfile, line_content)) // read file line by line
        {
            std::vector<std::string> Parsed_content = parse_string(line_content, " ");

            if (!Parsed_content[0].compare("input_graph"))
            // when it's equal, compare returns 0
            {
                int V = std::stoi(Parsed_content[2]);
                input_graph = dgraph_v_of_v<weight_type>(V);
            }
            else if (!Parsed_content[0].compare("Edge"))
            {
                int v1 = std::stoi(Parsed_content[1]);
                int v2 = std::stoi(Parsed_content[2]);
                weight_type ec = std::stod(Parsed_content[3]);
                input_graph.add_edge(v1, v2, ec);
            }
        }
        myfile.close(); // close the file
    }
    else
    {
        std::cout << "Unable to open file " << instance_name << std::endl
                  << "Please check the file location or file name." << std::endl;
        getchar();
        exit(1);
    }
}

template <typename weight_type>
void dgraph_read_dgraph_from_txt(std::string instance_name, dgraph_v_of_v<weight_type>& input_graph, int Jacard_or_Random)
{
    int vertex_num;
    input_graph.clear();
    std::string line_content;
    std::ifstream myfile(instance_name);
    if (myfile.is_open())
    {
        getline(myfile, line_content);
        int index1 = line_content.find(' ', 63);
        vertex_num = stoi(line_content.substr(63, index1-63));
        input_graph = dgraph_v_of_v<weight_type>(vertex_num);
        getline(myfile, line_content);
        while (getline(myfile, line_content)) // read file line by line
        {
            std::vector<std::string> Parsed_content = parse_string(line_content, " ");
            int v1 = std::stoi(Parsed_content[0]);
            int v2 = std::stoi(Parsed_content[1]);
            weight_type w;
            if (Jacard_or_Random) {
                w = std::stod(Parsed_content[3]);
            }
            else
                w = std::stod(Parsed_content[2]);

            input_graph.add_edge(v1, v2, w);
        }
        myfile.close(); // close the file
    }
    else
    {
        std::cout << "Unable to open file " << instance_name << std::endl
            << "Please check the file location or file name." << std::endl;
        getchar();
        exit(1);
    }
}

//测试binary读取的graph是否正确
bool test_dgraph_read_from_bin(dgraph_v_of_v<two_hop_weight_type>& input_graph, int Jacard_or_random ,int iteration_source_times)
{
    boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(input_graph.INs.size() - 1) };

    for (int yy = 0; yy < iteration_source_times; yy++) {
        int source = dist(boost_random_time_seed);

        int source_out_size = input_graph.OUTs[source].size();

        for (int xx = 0; xx < source_out_size; xx++) {
            int terminal = input_graph.OUTs[source][xx].first;
            two_hop_weight_type dis = input_graph.OUTs[source][xx].second;

            if (Jacard_or_random)
            {
                int s_cap_t = 0;
                int terminal_in_size = input_graph.INs[terminal].size();

                auto s_pointer = input_graph.OUTs[source].begin();
                auto s_end = input_graph.OUTs[source].end();

                auto t_pointer = input_graph.INs[terminal].begin();
                auto t_end = input_graph.INs[terminal].end();

                while (s_pointer != s_end && t_pointer != t_end)
                {
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

                int s_cup_t = terminal_in_size + source_out_size - s_cap_t;
                two_hop_weight_type ec = 1.0 - (two_hop_weight_type)s_cap_t / (s_cup_t);

                if (ec - dis < 1e-5)
                    continue;

                return false;
            }
            else {
                if (dis >= 0 && dis <= 100)
                    continue;

                return false;
            }            
        }
    }

    return true;
}

template <typename weight_type>
void dgraph_read_dgraph_from_bin(std::string instance_name, dgraph_v_of_v<weight_type>& input_graph)
{
    dgraph_hash_of_mixed_weighted_binary_read(instance_name, input_graph);

    //1 means test Jacard
    if (test_dgraph_read_from_bin(input_graph, 1, 100))
        cout << "by test " << instance_name << " is most likely right" << endl;
    else {
        cout << instance_name << " is not right" << endl;
        getchar();
    }
    
    return;
}