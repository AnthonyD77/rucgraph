#pragma once

/*the following codes are for testing

---------------------------------------------------
a cpp file (try.cpp) for running the following test code:
----------------------------------------

#include <iostream>
#include <fstream>
using namespace std;

#include <build_in_progress/HL/Hop/test_HBPLL.h>

int main()
{
	test_HBPLL();
}

------------------------------------------------------------------------------------------
Commends for running the above cpp file on Linux:

g++ -std=c++17 -I/__PATH__/boost_1_75_0 -I/__PATH__/rucgraph run.cpp -lpthread -Ofast -o A
./A
rm A 

(optional to put the above commends in run.sh, and then use the comment: sh run.sh)

*/

#include <build_in_progress/HL/Hop/graph_v_of_v_idealID_HB_v1.h>
#include <build_in_progress/HL/Hop/graph_v_of_v_idealID_HB_check_corectness.h>
#include <graph_v_of_v_idealID/random_graph/graph_v_of_v_idealID_generate_random_connected_graph.h>
#include <graph_v_of_v_idealID/read_save/graph_v_of_v_idealID_read.h>
#include <graph_v_of_v_idealID/read_save/graph_v_of_v_idealID_save.h>
#include <graph_v_of_v_idealID/graph_v_of_v_idaelID_sort.h>

void test_HBPLL() {
    /*parameters*/
    int iteration_graph_times = 1e2, iteration_source_times = 1000, iteration_terminal_times = 1000;
    int V = 20, E = 50, precision = 1, thread_num = 10;
    double ec_min = 0.1, ec_max = 1;
    bool weighted = (ec_min == 1 && ec_max == 1) ? false : true;
    graph_v_of_v_idealID_two_hop_case_info_v1 mm;

    /*control varible*/
    bool generate_new_graph = 1;
    bool print_time_details_in_every_loop = 0;
    bool print_label_before_canonical_fix = 0;
    bool print_L = 0;
    bool check_correctness = 1;
    check_path = 1;

    debug = 0;
    if (debug) {
        source_debug = 15;
        terminal_debug = 16;
        hop_cst_debug = 1;
        iteration_graph_times = 1;
        generate_new_graph = 0;
        iteration_source_times = 1;
        iteration_terminal_times = 1;
        print_L = 1;
    }

    /*hop bounded info*/
    int query_hop_cst = 5;
    mm.use_M = 0;
    mm.upper_k = 10;  // 0 means there is no limit
    mm.value_M = mm.use_M ? ec_max * E : 0;

    /*reduction method selection*/
    mm.use_2019R2 = 0;
    mm.use_enhanced2019R2 = 0;
    mm.use_non_adj_reduc_degree = 1;
    mm.print_label_before_canonical_fix = print_label_before_canonical_fix;
    mm.use_canonical_repair = true;

    /*result info*/
    double avg_index_time = 0, avg_index_size_per_v = 0, avg_MG_num = 0;
    double avg_canonical_repair_remove_label_ratio = 0;

    /*iteration*/
    for (int i = 0; i < iteration_graph_times; i++) {
        cout << ">>>iteration_graph_times: " << i << endl;

        /*input and output; below is for generating random new graph, or read saved graph*/

        graph_v_of_v_idealID instance_graph;
        if (generate_new_graph == 1) {
            instance_graph = graph_v_of_v_idealID_generate_random_connected_graph(V, E, ec_min, ec_max, precision, boost_random_time_seed);
            instance_graph = graph_v_of_v_idealID_sort(instance_graph);
            graph_v_of_v_idealID_save("simple_iterative_tests_HBPLL.txt", instance_graph);
        } else {
            graph_v_of_v_idealID_read("simple_iterative_tests_HBPLL.txt", instance_graph);
        }

        auto begin = std::chrono::high_resolution_clock::now();
        try {
            graph_v_of_v_idealID_HB_v1(instance_graph, V, weighted, thread_num, mm);

            if (print_time_details_in_every_loop) {
                cout << "mm.time_initialization: " << mm.time_initialization << "s" << endl;
                cout << "mm.time_reduction: " << mm.time_reduction << "s" << endl;
                cout << "mm.time_generate_labels: " << mm.time_generate_labels << "s" << endl;
                cout << "mm.time_update_predecessors: " << mm.time_update_predecessors << "s" << endl;
                cout << "mm.time_canonical_repair: " << mm.time_canonical_repair << "s" << endl;
                cout << "mm.time_update_labels: " << mm.time_update_labels << "s" << endl;
            }
        } catch (string s) {
            cout << s << endl;
            graph_v_of_v_idealID_two_hop_clear_global_values();
            continue;
        }

        auto end = std::chrono::high_resolution_clock::now();

        double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;  // s
        avg_index_time = avg_index_time + runningtime / iteration_graph_times;
        avg_MG_num = avg_MG_num + (double)mm.MG_num / iteration_graph_times;
        avg_canonical_repair_remove_label_ratio = avg_canonical_repair_remove_label_ratio + (double)mm.canonical_repair_remove_label_ratio / iteration_graph_times;

        if (print_L)
            mm.print_L();

        if (check_correctness)
            graph_v_of_v_idealID_HB_v1_check_correctness(mm, instance_graph, iteration_source_times, iteration_terminal_times, query_hop_cst);

        long long int index_size = 0;
        for (auto it = mm.L.begin(); it != mm.L.end(); it++) {
            index_size = index_size + (*it).size();
        }
        avg_index_size_per_v = avg_index_size_per_v + (double)index_size / V / iteration_graph_times;

        graph_v_of_v_idealID_two_hop_clear_global_values2();
        mm.clear_labels();
    }

    cout << "avg_index_time: " << avg_index_time << "s" << endl;
    cout << "avg_index_size_per_v: " << avg_index_size_per_v << endl;
    if (mm.use_2019R2 || mm.use_enhanced2019R2 || mm.use_non_adj_reduc_degree)
        cout << "avg_MG_num: " << avg_MG_num << endl;
    if (mm.use_canonical_repair)
        cout << "avg_canonical_repair_remove_label_ratio: " << avg_canonical_repair_remove_label_ratio << endl;
}
