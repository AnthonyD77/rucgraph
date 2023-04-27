#pragma once
#include <fstream>
#include <dgraph_v_of_v/dgraph_v_of_v.h>



/**
 * @brief 
 * 
 * @param instance_name : file name of storage dgraph
 * @param input_graph   : dgraph to be stored
 */
template <typename weight_type>
void dgraph_save_dgraph(std::string instance_name, dgraph_v_of_v<weight_type> &input_graph)
{

    std::ofstream outputFile;
    outputFile.precision(10);
    outputFile.setf(std::ios::fixed);
    outputFile.setf(std::ios::showpoint);
    outputFile.open(instance_name);

    outputFile << "SECTION Comments" << std::endl;
    outputFile << "Name \"" << instance_name << "\"" << std::endl;
    outputFile << "Creator \"save_dgraph\"" << std::endl;
    outputFile << "END" << std::endl;
    outputFile << std::endl;

    // print to file
    int V = input_graph.INs.size();
    outputFile << "input_graph |V|= " << V << " |E|= " << input_graph.edge_number() << std::endl;
    outputFile << std::endl;

    for (int i = 0; i < V; i++)
    {
        int out_size = input_graph.OUTs[i].size();
        for (int j = 0; j < out_size; j++)
        {
            outputFile << "Edge " << i << " " << input_graph.OUTs[i][j].first << " "
                       << input_graph.OUTs[i][j].second << '\n';
        }
    }
    outputFile << std::endl;

    outputFile << "EOF" << std::endl;
}