#pragma once
/*test function


main file content:
----------------------------------------------------------------

#include <iostream>
#include <string>
using namespace std;

// header files in the Boost library: https://www.boost.org/
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };
#include <build_in_progress/HL/two_hop_v2/test_CT_v2.h>

int main()
{
	test_CT();
}

---------------------------------------------------------------------------------------------

*/
#include <build_in_progress/HL/two_hop_v2_old/graph_hash_of_mixed_weighted_CT_v2.h>
#include <build_in_progress/HL/two_hop_v2_old/graph_hash_of_mixed_weighted_CT_extension.h>
#include <graph_hash_of_mixed_weighted/random_graph/graph_hash_of_mixed_weighted_generate_random_graph.h>
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_read_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/read_save/graph_hash_of_mixed_weighted_save_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/common_algorithms/graph_hash_of_mixed_weighted_shortest_paths.h>

pair<double, double> querying_element(graph_hash_of_mixed_weighted_CT_v2_case_info* case_info, int s, int t) {
    auto begin = std::chrono::high_resolution_clock::now();
    CT_extract_distance(*case_info, s, t);
    auto end = std::chrono::high_resolution_clock::now();
    double query_dis_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;  // s

    begin = std::chrono::high_resolution_clock::now();
    vector<pair<int, int>> path;
    CT_extract_path(*case_info, s, t, path);
    end = std::chrono::high_resolution_clock::now();
    double query_path_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;  // s

    return {query_dis_time, query_path_time};
}

pair<double, double> querying(graph_hash_of_mixed_weighted_CT_v2_case_info& case_info, vector<pair<int, int>>& query_list) {
    double query_dis_avg_time = 0, query_path_avg_time = 0;
    int query_times = query_list.size();

    ThreadPool pool(5);                                      // thread_num for querying
    std::vector<std::future<pair<double, double>>> results;  // return typename: xxx
    graph_hash_of_mixed_weighted_CT_v2_case_info* case_info_p = &case_info;
    for (int i = 0; i < query_times; i++) {
        int s = query_list[i].first, t = query_list[i].second;
        results.emplace_back(
            pool.enqueue([s, t, case_info_p] {  // pass const type value j to thread; [] can be empty
                return querying_element(case_info_p, s, t);
            }));
    }
    for (auto&& result : results) {
        pair<double, double> r = result.get();
        query_dis_avg_time += r.first;
        query_path_avg_time += r.second;
    }
    query_dis_avg_time = query_dis_avg_time / (double)query_times;
    query_path_avg_time = query_path_avg_time / (double)query_times;

    return {query_dis_avg_time, query_path_avg_time};
}


int debug = 0;
int ss = 9, tt = 8;

void test_CT() {

    /*parameters*/
    int iteration_graph_times = 1e2, iteration_source_times = 100, iteration_terminal_times = 100;
    int V = 10, E = 15, precision = 1;
    double ec_min = 0.1, ec_max = 1;  // since ec_min = 0.01, precision should be at least 2! Otherwise ec may be 0, and causes bugs in CT

    double avg_CT_time = 0, avg_PLL_time = 0;
    double avg_query_distance_time = 0, avg_query_path_time = 0;
    long long int avg_CT_index_bit_size = 0;

    int generate_new_graph = 1;
    int print_time_detail = 0;
    int check_correctness = 1;
    int print_L = 0;

    int use_CT_extension = 1;
    int case_info_d = 3;
    int case_info_use_P2H = 0;

    if (debug) {
        generate_new_graph = 0;
        iteration_graph_times = 1;
        iteration_source_times = 1;
        iteration_terminal_times = 1;
        print_L = 1;
    }

    /*iteration*/
    std::time_t now = std::time(0);
    boost::random::mt19937 gen{static_cast<std::uint32_t>(now)};
    for (int i = 0; i < iteration_graph_times; i++) {
        cout << ">>>>>iteration times: " << i << endl;

        std::unordered_set<int> generated_group_vertices;
        graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
        if (generate_new_graph == 1) {
            instance_graph = graph_hash_of_mixed_weighted_generate_random_graph(V, E, 0, 0, ec_min, ec_max, precision, boost_random_time_seed);
            graph_hash_of_mixed_weighted_save_graph_with_weight("simple_iterative_tests.txt", instance_graph, 0);
        }
        else {
            double lambda;
            graph_hash_of_mixed_weighted_read_graph_with_weight("simple_iterative_tests.txt", instance_graph, lambda);
        }

        graph_hash_of_mixed_weighted_CT_v2_case_info case_info;
        case_info.two_hop_case_info.use_2019R1 = 0;
        case_info.two_hop_case_info.use_2019R2 = 0;
        case_info.two_hop_case_info.use_enhanced2019R2 = 0;
        case_info.two_hop_case_info.use_non_adj_reduc_degree = 0;
        case_info.two_hop_case_info.use_canonical_repair = 0;
        case_info.d = case_info_d;
        case_info.use_PLL = 1;
        case_info.thread_num = 5;
        case_info.use_P2H = case_info_use_P2H;

        if (1) {
            auto begin = std::chrono::high_resolution_clock::now();

            if (!use_CT_extension) {
                CT_v2(instance_graph, V + 1, case_info);
            }
            else {
                CT_extension(instance_graph, V + 1, case_info);
            }

            auto end = std::chrono::high_resolution_clock::now();
            double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;  // s
            avg_CT_time += runningtime / iteration_graph_times;
            avg_CT_index_bit_size += case_info.compute_label_bit_size() / iteration_graph_times;
            avg_PLL_time += case_info.time_core_indexs / iteration_graph_times;

            if (print_L) {
                case_info.two_hop_case_info.print_L();
            }

            if (print_time_detail) {
                cout << "case_info.time1_initialization = " << case_info.time_initialization << "s" << endl;
                cout << "case_info.time2_tree_decomposition = " << case_info.time_tree_decomposition << "s" << endl;
                cout << "case_info.time3_tree_indexs = " << case_info.time_tree_indexs << "s" << endl;
                cout << "case_info.time_P2H = " << case_info.time_P2H << "s" << endl;
                cout << "case_info.time_extension = " << case_info.time_extension << "s" << endl;
                cout << "case_info.time4_lca = " << case_info.time_lca << "s" << endl;
                cout << "case_info.time5_core_indexs = " << case_info.time_core_indexs << "s" << endl;
                cout << "case_info.time6_post = " << case_info.time_post << "s" << endl;
            }
        }

        if (0) {
            vector<pair<int, int>> query_list;
            boost::random::uniform_int_distribution<> dist{static_cast<int>(0), static_cast<int>(V - 1)};
            for (int i = 0; i < 1e3; i++) {
                int source = dist(gen);
                int terminal = dist(gen);
                query_list.push_back({source, terminal});
            }
            querying(case_info, query_list);
        }

        /**
         * check correctness
         */

        if (check_correctness) {
            boost::random::uniform_int_distribution<> dist{static_cast<int>(0), static_cast<int>(V - 1)};

            for (int yy = 0; yy < iteration_source_times; yy++) {
                int source = dist(gen);

                if (debug) {
                    source = ss;
                }
                std::unordered_map<int, double> distances;
                std::unordered_map<int, int> predecessors;
                graph_hash_of_mixed_weighted_shortest_paths_source_to_all(instance_graph, source, distances, predecessors);

                for (int xx = 0; xx < iteration_terminal_times; xx++) {
                    int terminal = dist(gen);

                    if (debug) {
                        terminal = tt;
                    }

                    auto begin = std::chrono::high_resolution_clock::now();

                    double dis = CT_extract_distance(case_info, source, terminal);

                    auto end = std::chrono::high_resolution_clock::now();
                    avg_query_distance_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;

                    if (abs(dis - distances[terminal]) > 1e-5 && (dis < 1e10 || distances[terminal] < 1e10)) {
                        cout << "source = " << source << endl;
                        cout << "terminal = " << terminal << endl;

                        cout << "source vector:" << endl;
                        for (auto it = case_info.two_hop_case_info.L[source].begin();
                             it != case_info.two_hop_case_info.L[source].end(); it++) {
                            cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
                        }
                        cout << endl;
                        cout << "terminal vector:" << endl;
                        for (auto it = case_info.two_hop_case_info.L[terminal].begin();
                             it != case_info.two_hop_case_info.L[terminal].end(); it++) {
                            cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
                        }
                        cout << endl;

                        cout << "dis = " << dis << endl;
                        cout << "distances[terminal] = " << distances[terminal] << endl;
                        cout << "abs(dis - distances[terminal]) > 1e-5!" << endl;
                        getchar();
                    }

                    vector<pair<int, int>> path;

                    begin = std::chrono::high_resolution_clock::now();

                    if (!use_CT_extension) {
                        CT_extract_path(case_info, source, terminal, path);
                    }
                    else {
                        CT_extension_extract_path(instance_graph, case_info, source, terminal, path);
                    }
                    end = std::chrono::high_resolution_clock::now();

                    avg_query_path_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;

                    double path_dis = 0;
                    if (path.size() == 0) {
                        if (source != terminal) {  // disconnected
                            path_dis = std::numeric_limits<double>::max();
                        }
                    }
                    else {
                        for (auto it = path.begin(); it != path.end(); it++) {
                            path_dis = path_dis +
                                       graph_hash_of_mixed_weighted_edge_weight(instance_graph, it->first, it->second);
                            if (path_dis > std::numeric_limits<double>::max()) {
                                path_dis = std::numeric_limits<double>::max();
                            }
                        }
                    }
                    if (abs(dis - path_dis) > 1e-5 && (dis < 1e10 || distances[terminal] < 1e10)) {
                        cout << "source = " << source << endl;
                        cout << "terminal = " << terminal << endl;

                        cout << "source vector:" << endl;
                        for (auto it = case_info.two_hop_case_info.L[source].begin();
                             it != case_info.two_hop_case_info.L[source].end(); it++) {
                            cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
                        }
                        cout << endl;
                        cout << "terminal vector:" << endl;
                        for (auto it = case_info.two_hop_case_info.L[terminal].begin();
                             it != case_info.two_hop_case_info.L[terminal].end(); it++) {
                            cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
                        }
                        cout << endl;
                        print_vector_pair_int(path);
                        for (int i = 0; i < path.size(); i++) {
                            if (graph_hash_of_mixed_weighted_contain_edge(instance_graph, path[i].first,
                                                                          path[i].second) == false) {
                                cout << "not-existing edge: " << path[i].first << "," << path[i].second << endl;
                                cout << "isIntree[path[i].first]: " << case_info.isIntree[path[i].first] << endl;
                                cout << "isIntree[path[i].second]: " << case_info.isIntree[path[i].second] << endl;
                            }
                        }
                        cout << "dis = " << dis << endl;
                        cout << "distances[terminal] = " << distances[terminal] << endl;
                        cout << "path_dis = " << path_dis << endl;
                        cout << "abs(dis - path_dis) > 1e-5!" << endl;
                        getchar();
                    }
                }
            }
        }
    }

    avg_query_distance_time /= (iteration_source_times * iteration_terminal_times);
    avg_query_path_time /= (iteration_source_times * iteration_terminal_times);

    cout << "avg_CT_time: " << avg_CT_time << endl;
    cout << "avg_PLL_time: " << avg_PLL_time << "s" << endl;
    cout << "avg_CT_index_bit_size: " << avg_CT_index_bit_size << endl;
    cout << "avg_query_distance_time: " << avg_query_distance_time << endl;
    cout << "avg_query_path_time: " << avg_query_path_time << endl;
}