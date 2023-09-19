#pragma once

#include <text_mining/print_items.h>
#include <boost/random.hpp>
#include <build_in_progress/HL/Hop/graph_v_of_v_idealID_HB_two_hop_query.h>

boost::random::mt19937 boost_random_time_seed{static_cast<std::uint32_t>(std::time(0))};

/* 
    HB-HFS
    this func get the all vertices distances from vertex source with hop constraint hop_cst  
*/
void graph_v_of_v_HB_shortest_distance(graph_v_of_v_idealID &instance_graph, int source, int hop_cst, vector<double> &distance) {
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

    /* BFS */
    while (1) {
        if (h > hop_cst || Q[h].empty())
            break;

        for (auto it = Q[h].begin(); it != Q[h].end(); it++) {
            int v = it->first;
            double distance_v = it->second;
            if (distance[v] != 0 && distance[v] < distance_v)
                continue;
            distance[v] = distance_v;

            int v_adj_size = instance_graph[v].size();
            for (int i = 0; i < v_adj_size; i++) {
                int adj_v = instance_graph[v][i].first;
                double ec = instance_graph[v][i].second;
                if (distance_v + ec < distance[adj_v]) {
                    Q[h + 1].push_back({adj_v, distance_v + ec});
                }
            }
        }
        h++;
    }
}

bool debug = 0;
int source_debug = 0;
int terminal_debug = 0;
int hop_cst_debug = 0;
bool check_path = 0;

void graph_v_of_v_idealID_HB_v1_check_correctness(graph_v_of_v_idealID_two_hop_case_info_v1 &case_info, graph_v_of_v_idealID &instance_graph, int iteration_source_times, int iteration_terminal_times, int hop_cst) {
    /*below is for checking whether the above labels are right (by randomly computing shortest paths)

	this function can only be used when 0 to n-1 is in the graph, i.e., the graph is an ideal graph

	*/

    boost::random::uniform_int_distribution<> vertex_range{static_cast<int>(0), static_cast<int>(instance_graph.size() - 1)};
    boost::random::uniform_int_distribution<> hop_range{static_cast<int>(1), static_cast<int>(10)};

    for (int yy = 0; yy < iteration_source_times; yy++) {
        int source = vertex_range(boost_random_time_seed);
        std::vector<double> distances;
        distances.resize(instance_graph.size());
        std::vector<int> predecessors;
        predecessors.resize(instance_graph.size());

        hop_cst = hop_range(boost_random_time_seed);

        if (debug) {
            source = source_debug;
            hop_cst = hop_cst_debug;
        }

        graph_v_of_v_HB_shortest_distance(instance_graph, source, hop_cst, distances);

        for (int xx = 0; xx < iteration_terminal_times; xx++) {
            int terminal = vertex_range(boost_random_time_seed);

            if (debug)
                terminal = terminal_debug;

            double dis;
            if (case_info.use_2019R2 || case_info.use_enhanced2019R2 || case_info.use_non_adj_reduc_degree)
                dis = graph_v_of_v_idealID_two_hop_v1_extract_distance_st_no_R1(case_info.L, case_info.reduction_measures_2019R2, source, terminal, hop_cst, case_info.value_M);
            else
                dis = graph_v_of_v_idealID_two_hop_v1_extract_distance_no_reduc(case_info.L, source, terminal, hop_cst, case_info.value_M);

            if (abs(dis - distances[terminal]) > 1e-4 && (dis < std::numeric_limits<double>::max() || distances[terminal] < std::numeric_limits<double>::max())) {
                cout << "source = " << source << endl;
                cout << "terminal = " << terminal << endl;
                cout << "hop_cst = " << hop_cst << endl;
                cout << "source vector:" << endl;
                for (auto it = case_info.L[source].begin(); it != case_info.L[source].end(); it++) {
                    cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << "," << it->hop << ">";
                }
                cout << endl;
                cout << "terminal vector:" << endl;
                for (auto it = case_info.L[terminal].begin(); it != case_info.L[terminal].end(); it++) {
                    cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << "," << it->hop << ">";
                }
                cout << endl;

                cout << "dis = " << dis << endl;
                cout << "distances[terminal] = " << distances[terminal] << endl;
                cout << "abs(dis - distances[terminal]) > 1e-5!" << endl;
                getchar();
            }

            if (check_path) {
                vector<pair<int, int>> path = graph_v_of_v_idealID_two_hop_v1_extract_shortest_path_st_no_R1(case_info.L, case_info.reduction_measures_2019R2, source, terminal, hop_cst, case_info.value_M);

                double path_dis = 0;
                if (path.size() == 0) {
                    if (source != terminal) {
                        path_dis = std::numeric_limits<double>::max();
                    }
                } else {
                    for (auto it = path.begin(); it != path.end(); it++) {
                        path_dis = path_dis + graph_v_of_v_idealID_edge_weight(instance_graph, it->first, it->second);
                        if (path_dis > std::numeric_limits<double>::max()) {
                            path_dis = std::numeric_limits<double>::max();
                        }
                    }
                }
                if (abs(dis - path_dis) > 1e-4 && (dis < std::numeric_limits<double>::max() || distances[terminal] < std::numeric_limits<double>::max())) {
                    cout << "source = " << source << endl;
                    cout << "terminal = " << terminal << endl;
                    cout << "hop_cst = " << hop_cst << endl;

                    cout << "source vector:" << endl;
                    for (auto it = case_info.L[source].begin(); it != case_info.L[source].end(); it++) {
                        cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << "," << it->hop << ">";
                    }
                    cout << endl;
                    cout << "terminal vector:" << endl;
                    for (auto it = case_info.L[terminal].begin(); it != case_info.L[terminal].end(); it++) {
                        cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << "," << it->hop << ">";
                    }
                    cout << endl;

                    print_vector_pair_int(path);
                    cout << "dis = " << dis << endl;
                    cout << "path_dis = " << path_dis << endl;
                    cout << "abs(dis - path_dis) > 1e-5!" << endl;
                    getchar();
                }
            }
        }
    }
}