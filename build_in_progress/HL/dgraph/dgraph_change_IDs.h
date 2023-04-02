#pragma once

#include <dgraph_v_of_v/dgraph_v_of_v.h>

template <typename weight_type>
void sort_graph(dgraph_v_of_v<weight_type>& graph)
{
    int N = graph.INs.size();
    std::vector<std::vector<std::pair<int, weight_type>>> New_INs(N);
    std::vector<std::vector<std::pair<int, weight_type>>> New_OUTs(N);

    vector<pair<int, int>> value(N);

    for (int i = 0; i < N; i++)
        //change to in_degree\outdegree
        value[i] = { graph.OUTs[i].size() + graph.INs[i].size() , i };

    //升序排列
    sort(value.begin(), value.end());

    vector<int> old_to_new(N);

    for (int j = 0; j < N; j++)
    {
        old_to_new[value[j].second] = N - 1 - j;
    }

    for (int i = 0; i < N; i++)
    {
        int in_size = graph.INs[i].size();
        std::vector<std::pair<int, weight_type>> v1 = graph.INs[i];

        int out_size = graph.OUTs[i].size();
        std::vector<std::pair<int, weight_type>> v2 = graph.OUTs[i];

        for (int j = 0; j < in_size; j++)
        {
            v1[j].first = old_to_new[v1[j].first];
        }
        New_INs[old_to_new[i]] = v1;


        for (int j = 0; j < out_size; j++)
        {
            v2[j].first = old_to_new[v2[j].first];
        }
        New_OUTs[old_to_new[i]] = v2;
    }

    graph.INs = New_INs;
    graph.OUTs = New_OUTs;

    return;
}
