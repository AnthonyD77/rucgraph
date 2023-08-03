#pragma once

/*the following codes are for testing

---------------------------------------------------
a cpp file (try.cpp) for running the following test code:
----------------------------------------

#include <iostream>
#include <fstream>
using namespace std;

// header files in the Boost library: https://www.boost.org/
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

#include <build_in_progress/HL/two_hop_v1/test_PLL_PSL_v1.h>


int main()
{
	test_PLL_PSL();
}

------------------------------------------------------------------------------------------
Commends for running the above cpp file on Linux:

g++ -std=c++17 -I/home/boost_1_75_0 -I/root/rucgraph try.cpp -lpthread -Ofast -o A
./A
rm A

(optional to put the above commends in run.sh, and then use the comment: sh run.sh)


*/

#include "build_in_progress/HL/Hop/graph_hash_of_mixed_weighted_two_hop_labels_v1.h"
#include "graph_v_of_v_idealID/graph_v_of_v_idealID.h"
#include <build_in_progress/HL/Hop/graph_hash_of_mixed_weighted_HB_v1.h>
#include <build_in_progress/HL/Hop/graph_hash_of_mixed_weighted_HB_shortest_path.h>
#include <cstdio>
#include <cstdlib>
#include <graph_v_of_v_idealID/random_graph/graph_v_of_v_idealID_generate_random_connected_graph.h>
#include <graph_v_of_v_idealID/read_save/graph_v_of_v_idealID_read.h>
#include <graph_v_of_v_idealID/read_save/graph_v_of_v_idealID_save.h>
#include <graph_v_of_v_idealID/graph_v_of_v_idaelID_sort.h>
#include <graph_v_of_v_idealID/common_algorithms/graph_v_of_v_idealID_shortest_paths.h>

#include <graph_hash_of_mixed_weighted/common_algorithms/graph_hash_of_mixed_weighted_shortest_paths.h>

#include <boost/random.hpp>
#include <numeric>
#include <ostream>
boost::random::mt19937 boost_random_time_seed{static_cast<std::uint32_t>(std::time(0))};

bool debug = 0;
int source_debug = 0;
int terminal_debug = 0; 
bool check_path = 0;

void graph_hash_of_mixed_weighted_HB_v1_check_correctness(graph_hash_of_mixed_weighted_two_hop_case_info_v1 &case_info, graph_v_of_v_idealID &instance_graph, int iteration_source_times, int iteration_terminal_times, int hop_cst)
{

    /*below is for checking whether the above labels are right (by randomly computing shortest paths)

	this function can only be used when 0 to n-1 is in the graph, i.e., the graph is an ideal graph

	*/

    boost::random::uniform_int_distribution<> vertex_range{static_cast<int>(0), static_cast<int>(instance_graph.size() - 1)};

    for (int yy = 0; yy < iteration_source_times; yy++)
    {
        int source = vertex_range(boost_random_time_seed);
        std::vector<double> distances;
        distances.resize(instance_graph.size());
        std::vector<int> predecessors;
        predecessors.resize(instance_graph.size());

        if (debug)
            source = source_debug;

        graph_hash_of_mixed_weighted_HB_shortest_distance(instance_graph, source, hop_cst, distances);
        // graph_v_of_v_idealID_shortest_paths(instance_graph, source, distances, predecessors);

        for (int xx = 0; xx < iteration_terminal_times; xx++)
        {
            int terminal = vertex_range(boost_random_time_seed);

            if(debug)
                terminal = terminal_debug;
            
            double dis;
            if (case_info.use_2019R2 || case_info.use_enhanced2019R2 || case_info.use_non_adj_reduc_degree)
                dis = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(case_info.L, case_info.reduction_measures_2019R2, source, terminal, hop_cst);
            else
                dis = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(case_info.L, source, terminal, hop_cst);

            if (abs(dis - distances[terminal]) > 1e-4 && (dis < std::numeric_limits<double>::max() || distances[terminal] < std::numeric_limits<double>::max()))
            {
                cout << "source = " << source << endl;
                cout << "terminal = " << terminal << endl;
                cout << "source vector:" << endl;
                for (auto it = case_info.L[source].begin(); it != case_info.L[source].end(); it++)
                {
                    cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << "," << it->hop << ">";
                }
                cout << endl;
                cout << "terminal vector:" << endl;
                for (auto it = case_info.L[terminal].begin(); it != case_info.L[terminal].end(); it++)
                {
                    cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << "," << it->hop << ">";
                }
                cout << endl;

                cout << "dis = " << dis << endl;
                cout << "distances[terminal] = " << distances[terminal] << endl;
                cout << "abs(dis - distances[terminal]) > 1e-5!" << endl;
                getchar();
            }

            if (check_path)
            {
                vector<pair<int, int>> path = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(case_info.L, case_info.reduction_measures_2019R2, source, terminal, hop_cst);

                double path_dis = 0;
                if (path.size() == 0)
                {
                    if (source != terminal)
                    {
                        path_dis = std::numeric_limits<double>::max();
                    }
                }
                else
                {
                    for (auto it = path.begin(); it != path.end(); it++)
                    {
                        path_dis = path_dis + graph_v_of_v_idealID_edge_weight(instance_graph, it->first, it->second);
                        if (path_dis > std::numeric_limits<double>::max())
                        {
                            path_dis = std::numeric_limits<double>::max();
                        }
                    }
                }
                if (abs(dis - path_dis) > 1e-4 && (dis < std::numeric_limits<double>::max() || distances[terminal] < std::numeric_limits<double>::max()))
                {
                    cout << "source = " << source << endl;
                    cout << "terminal = " << terminal << endl;

                    cout << "source vector:" << endl;
                    for (auto it = case_info.L[source].begin(); it != case_info.L[source].end(); it++)
                    {
                        cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
                    }
                    cout << endl;
                    cout << "terminal vector:" << endl;
                    for (auto it = case_info.L[terminal].begin(); it != case_info.L[terminal].end(); it++)
                    {
                        cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
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

void test_HBPLL()
{
    /*parameters*/
    int iteration_graph_times = 1e2, iteration_source_times = 1000, iteration_terminal_times = 1000;
    int V = 1e2, E = 5e2, precision = 1, thread_num = 10;
    double ec_min = 0.1, ec_max = 1;
    bool weighted = (ec_min == 1 && ec_max == 1) ? false : true;
    graph_hash_of_mixed_weighted_two_hop_case_info_v1 mm;
    bool use_PLL = 1; // 1: PLL 0: PSL
    int query_hop_cst = 2;

    /*control varible*/
    bool generate_new_graph = 1;
    bool print_time_details_in_every_loop = 0;
    bool print_label_before_canonical_fix = 0;
    bool print_L = 0;
    bool check_correctness = 1;
    check_path = 1;

    debug = 0;
    if (debug)
    {
        source_debug = 1;
        terminal_debug = 0;
        iteration_graph_times = 1;
        generate_new_graph = 0;
        iteration_source_times = 1;
        iteration_terminal_times = 1;
        print_L = 1;
    } 

    /*hop bounded upper limit*/
    mm.upper_k = 10; // 0 means there is no limit
    mm.use_hbdij = 1;

    /*reduction method selection*/
        /* use_hb 目前不支持 R2 */
    mm.use_2019R2 = 0;
    mm.use_enhanced2019R2 = 0;
    mm.use_non_adj_reduc_degree = 0;
    mm.max_degree_MG_enhanced2019R2 = 100;
    mm.max_labal_size = 6e9;
    mm.max_run_time_seconds = 1e9;
    mm.print_label_before_canonical_fix = print_label_before_canonical_fix;
    mm.use_canonical_repair = true;

    /*result info*/
    double avg_index_time = 0, avg_index_size_per_v = 0, avg_reduce_V_num_2019R1 = 0, avg_MG_num = 0;
    double avg_canonical_repair_remove_label_ratio = 0;

    /*iteration*/
    for (int i = 0; i < iteration_graph_times; i++)
    {
        cout << ">>>iteration_graph_times: " << i << endl;

        /*input and output; below is for generating random new graph, or read saved graph*/

        std::unordered_set<int> generated_group_vertices;
        graph_v_of_v_idealID instance_graph;
        if (generate_new_graph == 1)
        {
            instance_graph = graph_v_of_v_idealID_generate_random_connected_graph(V, E, ec_min, ec_max, precision, boost_random_time_seed);
            instance_graph = graph_v_of_v_idealID_sort(instance_graph);
            graph_v_of_v_idealID_save("simple_iterative_tests_HBPLL.txt", instance_graph);
        }
        else
        {
            graph_v_of_v_idealID_read("simple_iterative_tests_HBPLL.txt", instance_graph);
        }

        auto begin = std::chrono::high_resolution_clock::now();
        try
        {
            if (use_PLL)
            {
                graph_hash_of_mixed_weighted_HB_v1(instance_graph, V, weighted, thread_num, mm);
            }
            else
            {
            }
            if (print_time_details_in_every_loop)
            {
                cout << "mm.time_initialization: " << mm.time_initialization << "s" << endl;
                cout << "mm.time_2019R1: " << mm.time_2019R1 << "s" << endl;
                cout << "mm.time_2019R2_or_enhanced_pre: " << mm.time_2019R2_or_enhanced_pre << "s" << endl;
                cout << "mm.time_2019R2_or_enhanced_fixlabels: " << mm.time_2019R2_or_enhanced_fixlabels << "s" << endl;
                cout << "mm.time_generate_labels: " << mm.time_generate_labels << "s" << endl;
                cout << "mm.time_canonical_repair1: " << mm.time_canonical_repair1 << "s" << endl;
                cout << "mm.time_canonical_repair2: " << mm.time_canonical_repair2 << "s" << endl;
                cout << "mm.time_update_old_IDs_in_labels: " << mm.time_update_old_IDs_in_labels << "s" << endl;
            }
        }
        catch (string s)
        {
            cout << s << endl;
            graph_hash_of_mixed_weighted_two_hop_clear_global_values();
            continue;
        }

        auto end = std::chrono::high_resolution_clock::now();

        double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
        avg_index_time = avg_index_time + runningtime / iteration_graph_times;
        avg_reduce_V_num_2019R1 = avg_reduce_V_num_2019R1 + (double)mm.reduce_V_num_2019R1 / iteration_graph_times;
        avg_MG_num = avg_MG_num + (double)mm.MG_num / iteration_graph_times;
        avg_canonical_repair_remove_label_ratio = avg_canonical_repair_remove_label_ratio + (double)mm.canonical_repair_remove_label_ratio / iteration_graph_times;

        if (print_L)
            mm.print_L();

        /*debug*/
        if (0)
        {
            // graph_hash_of_mixed_weighted_print(instance_graph);
            mm.print_L();
            mm.print_reduction_measures_2019R1();
            mm.print_reduction_measures_2019R2();
            mm.print_f_2019R1();
        }

        /*test canonical_repair proof*/
        if (0)
        {
            auto mm2 = mm;
            if (!use_PLL)
            {                                                                                   // use different method here
                // graph_hash_of_mixed_weighted_HB_v1(instance_graph, V + 1, weighted, 1, mm2); // single thread
            }
            else
            {
            }

            auto &L1 = mm.L;
            auto &L2 = mm2.L;
            int size = L1.size();
            bool L1_is_l2 = true;
            for (int xx = 0; xx < size; xx++)
            {
                if (L1[xx].size() != L2[xx].size())
                {
                    if ((L1[xx].size() == 0 && L2[xx].size() == 1) || (L1[xx].size() == 1 && L2[xx].size() == 0))
                    {
                    }
                    else
                    {
                        L1_is_l2 == false;
                        cout << "here" << endl;
                        mm.print_L();
                        mm2.print_L();
                        getchar();
                    }
                }
                else
                {
                    int size2 = L1[xx].size();
                    for (int yy = 0; yy < size2; yy++)
                    {
                        if (L1[xx][yy].vertex != L2[xx][yy].vertex)
                        {
                            L1_is_l2 == false;
                            cout << "here" << endl;
                            mm.print_L();
                            mm2.print_L();
                            getchar();
                        }
                    }
                }
            }
        }

        if (check_correctness)
            graph_hash_of_mixed_weighted_HB_v1_check_correctness(mm, instance_graph, iteration_source_times, iteration_terminal_times, query_hop_cst);

        long long int index_size = 0;
        for (auto it = mm.L.begin(); it != mm.L.end(); it++)
        {
            index_size = index_size + (*it).size();
        }
        avg_index_size_per_v = avg_index_size_per_v + (double)index_size / V / iteration_graph_times;
        
        graph_hash_of_mixed_weighted_two_hop_clear_global_values2();
        mm.clear_labels();
    }

    cout << "avg_index_time: " << avg_index_time << "s" << endl;
    cout << "avg_index_size_per_v: " << avg_index_size_per_v << endl;
    if (mm.use_2019R1)
        cout << "avg_reduce_V_num_2019R1: " << avg_reduce_V_num_2019R1 << endl;
    if (mm.use_2019R2 || mm.use_enhanced2019R2)
        cout << "avg_MG_num: " << avg_MG_num << endl;
    if (mm.use_canonical_repair)
        cout << "avg_canonical_repair_remove_label_ratio: " << avg_canonical_repair_remove_label_ratio << endl;
}
