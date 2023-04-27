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
void dgraph_read_dgraph(std::string instance_name, dgraph_v_of_v<weight_type>& input_graph)
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
