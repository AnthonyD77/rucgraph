#pragma once
using namespace std;
/*
    to include this file, please copy
    #include <build_in_progress/HL/dgraph/dHL_base.h>
*/

template <typename weight_type>
class d_2hoplabel {
public:
	int vertex; // std::vector<PLL_with_non_adj_reduction_sorted_label> should be sorted by vertex
	weight_type distance;
};

template <typename weight_type>
class d_2hoplabels {
public:
	std::vector<std::vector<d_2hoplabel<weight_type>>> INlabels, OUTlabels;
};

/*
example:
------------------------------------
#include <build_in_progress/HL/dgraph/dHL_base.h>
int main()
{
	d_2hoplabels<float> L; // ³õÊ¼»¯
}
----------------------
*/

void label_output_to_file(std::string instance_name, d_2hoplabels<float> L)
{
    std::ofstream outputFile;
    outputFile.open(instance_name);

    outputFile << "in label" << std::endl;

    int size1 = L.INlabels.size();
    for (int i = 0; i < size1; i++)
    {
        outputFile << i << ":\t";
        int size2 = L.INlabels[i].size();
        for (int j = 0; j < size2; j++)
        {
            outputFile << "(" << L.INlabels[i][j].vertex << "," << L.INlabels[i][j].distance << ") ";
        }
        outputFile << std::endl;
    }

    outputFile << std::endl << std::endl << std::endl << "out label" << std::endl;;

    size1 = L.OUTlabels.size();
    for (int i = 0; i < size1; i++)
    {
        outputFile << i << ":\t";
        int size2 = L.OUTlabels[i].size();
        for (int j = 0; j < size2; j++)
        {
            outputFile << "(" << L.OUTlabels[i][j].vertex << "," << L.OUTlabels[i][j].distance << ") "; 
        }
        outputFile << std::endl;
    }
}