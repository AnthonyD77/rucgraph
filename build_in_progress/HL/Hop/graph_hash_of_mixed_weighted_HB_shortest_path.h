#pragma once
#include "boost/assert/source_location.hpp"
#include <build_in_progress/HL/Hop/graph_hash_of_mixed_weighted_two_hop_labels_v1.h>
#include <build_in_progress/HL/Hop/graph_hash_of_mixed_weighted_HB_v1.h>
#include <vector>

bool print_info = 0;

/* this func get the all vertices distances from vertex source with hop constraint hop_cst  */
void graph_hash_of_mixed_weighted_HB_shortest_distance(graph_v_of_v_idealID &instance_graph, int source, int hop_cst, vector<double> &distance)
{
    /* instance_graph is the same with ideal_graph before reduction */
    int N = instance_graph.size();
    int v_k = source;

    vector<vector<pair<int, double>>> Q;
    Q.resize(hop_cst + 2);
    Q[0].push_back({v_k, 0});

    distance.resize(N);
    for (int i = 0; i < N; i++)
        distance[i] = (i == v_k) ? 0 : std::numeric_limits<double>::max();

    int h = 0;

    if (print_info)
        cout << "-----" << source << "-----" << endl;

    /* BFS */
    while (1)
    {
        if (h > hop_cst || Q[h].empty())
            break;

        for (auto it = Q[h].begin(); it != Q[h].end(); it++)
        {
            int v = it->first;
            double distance_v = it->second;
            if (distance[v] != 0 && distance[v] < distance_v)
                continue;
            distance[v] = distance_v;
            if (print_info)
                cout << "update " << v << " with " << distance_v << endl;

            int v_adj_size = instance_graph[v].size();
            for (int i = 0; i < v_adj_size; i++)
            {
                int adj_v = instance_graph[v][i].first;
                double ec = instance_graph[v][i].second;
                if (distance_v + ec < distance[adj_v])
                {
                    if (print_info)
                        cout << "push neighber " << adj_v << " with " << distance_v + ec << endl;
                    Q[h+1].push_back({adj_v, distance_v + ec});
                }
            }
        }
        h++;
    }

}