#pragma once

/*the following codes are for testing

---------------------------------------------------
a cpp file (try.cpp) for running the following test code:
----------------------------------------

#include <iostream>
#include <fstream>
using namespace std;
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

#include <dgraph_v_of_v/test_dgraph.h>

int main()
{
    test_dgraph();
}

------------------------------------------------------------------------------------------
Commends for running the above cpp file on Linux:

g++ -std=c++17 -I/home/boost_1_75_0 -I/home/dengs/dgraph try.cpp -lpthread -Ofast -o A
./A
rm A

(optional to put the above commends in run.sh, and then use the comment: sh run.sh)


*/

#include <dgraph_v_of_v/dgraph_v_of_v.h>
#include <build_in_progress/HL/dgraph/dgraph_generate_random_dgraph.h>
#include <build_in_progress/HL/dgraph/dgraph_save_dgraph_with_weight.h>
#include <build_in_progress/HL/dgraph/dgraph_read_dgraph_with_weight.h>
#include <build_in_progress/HL/dgraph/dgraph_change_IDs.h>
#include <build_in_progress/HL/dgraph/dgraph_PLL.h>
#include <build_in_progress/HL/dgraph/dgraph_PSL.h>
#include <build_in_progress/HL/dgraph/dHL_base.h>
#include <build_in_progress/HL/dgraph/dgraph_two_hop_label.h>
#include <build_in_progress/HL/dgraph/dgraph_shortest_paths.h>
#include <unordered_map>

void print_vector_pair_int(vector<pair<int, int>>& input_vector)
{

    std::cout << "print_vector_pair_int:" << std::endl;
    for (int i = 0; i < input_vector.size(); i++)
    {
        std::cout << "item: |" << input_vector[i].first << "," << input_vector[i].second << "|" << std::endl;
    }
}

void dgraph_v1_check_correctness(dgraph_case_info_v1& case_info, dgraph_v_of_v<two_hop_weight_type>& instance_graph,
    int iteration_source_times, int iteration_terminal_times)
{
    /*below is for checking whether the above labels are right (by randomly computing shortest paths)

    this function can only be used when 0 to n-1 is in the graph, i.e., the graph is an ideal graph

    */

    boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(instance_graph.INs.size() - 1) };

    // graph_hash_of_mixed_weighted_print(instance_graph);

    for (int yy = 0; yy < iteration_source_times; yy++)
    {
        int source = dist(boost_random_time_seed);
        unordered_map<int, two_hop_weight_type> distances;
        source = 3; //cout << "source = " << source << endl;

        dgraph_shortest_paths_source_to_all(instance_graph, source, distances);

        for (int xx = 0; xx < iteration_terminal_times; xx++)
        {

            int terminal = dist(boost_random_time_seed);

            terminal = 4; //cout << "terminal = " << terminal << endl;

            two_hop_weight_type dis =
                dgraph_v1_extract_shortest_distance(
                    case_info.L_in, case_info.L_out, instance_graph, source, terminal);

            if (abs(dis - distances[terminal]) > 1e-4 &&
                (dis < std::numeric_limits<double>::max() || distances[terminal] < std::numeric_limits<double>::max()))
            {
                cout << "source = " << source << endl;
                cout << "terminal = " << terminal << endl;
                cout << "source vector:" << endl;
                for (auto it = case_info.L_out[source].begin(); it != case_info.L_out[source].end(); it++)
                {
                    cout << "<" << it->vertex << "," << it->distance << ">";
                }
                cout << endl;
                cout << "terminal vector:" << endl;
                for (auto it = case_info.L_in[terminal].begin(); it != case_info.L_in[terminal].end(); it++)
                {
                    cout << "<" << it->vertex << "," << it->distance << ">";
                }
                cout << endl;

                cout << "dis = " << dis << endl;
                cout << "distances[terminal] = " << distances[terminal] << endl;
                cout << "abs(dis - distances[terminal]) > 1e-5!" << endl;
                getchar();
            }

            // cout << 0 << endl;
            // cout << source << " " << terminal << endl;
            // getchar();

        }
    }

}

void test_dgraph()
{
    /*parameters*/
    int iteration_graph_times = 100, iteration_source_times = 100, iteration_terminal_times = 100;

    int V = 1000, E = 5000, precision = 1, thread_num = 5;
    //int V = 5, E = 7, precision = 1, thread_num = 3;

    two_hop_weight_type ec_min = 0.1, ec_max = 1; // set ec_min=ec_max=1 for testing unweighted PLL_with_non_adj_reduction

    double avg_index_time = 0, avg_index_size_per_v = 0;

    bool use_PLL = 1; // 1: PLL 0: PSL

    /*reduction method selection*/
    dgraph_case_info_v1 mm;

    mm.use_canonical_repair = 1;

    /*iteration*/
    for (int i = 0; i < iteration_graph_times; i++)
    {
        cout << i << endl;

        /*input and output; below is for generating random new graph, or read saved graph*/
        // int generate_new_graph = 1;
        int generate_new_graph = 1;

        dgraph_v_of_v<two_hop_weight_type> instance_graph;

        if (generate_new_graph == 1)
        {
            //id应该是排序好的
            instance_graph = dgraph_generate_random_dgraph(V, E, ec_min, ec_max, precision, boost_random_time_seed);
            sort_graph(instance_graph);
            dgraph_save_dgraph_with_weight("random_dgraph_test.txt", instance_graph);
        }
        else
        {
            dgraph_read_dgraph_with_weight("PLL_test.txt", instance_graph);
            //read_dgraph_with_weight("random_dgraph_test.txt", instance_graph);
            //read_dgraph_with_weight("canonical_repair_text.txt", instance_graph);
            //read_dgraph_with_weight("canonical_repair_text_2.txt", instance_graph);
        }


        //cout << "begin sort" << endl;
        //sort_graph(instance_graph);
        //cout << "after sort" << endl;

        auto begin = std::chrono::high_resolution_clock::now();
        try
        {
            if (use_PLL)
            {
                dgraph_PLL(instance_graph, V, thread_num, mm);
                //dgraph_PLL(instance_graph, V, thread_num, mm, 1);
            }
            else
                //dgraph_PSL(instance_graph, V, 1, mm);
                //dgraph_PSL_v2(instance_graph, V, 3, mm);
                dgraph_PSL_v3(instance_graph, V, thread_num, mm);

            if (0)
            {
                //cout << "mm.time_initialization: " << mm.time_initialization << "s" << endl;
                //cout << "mm.time_2019R1: " << mm.time_2019R1 << "s" << endl;
                //cout << "mm.time_2019R2_or_enhanced_pre: " << mm.time_2019R2_or_enhanced_pre << "s" << endl;
                //cout << "mm.time_2019R2_or_enhanced_fixlabels: " << mm.time_2019R2_or_enhanced_fixlabels << "s" << endl;
                //cout << "mm.time_generate_labels: " << mm.time_generate_labels << "s" << endl;
                //cout << "mm.time_canonical_repair1: " << mm.time_canonical_repair1 << "s" << endl;
                //cout << "mm.time_canonical_repair2: " << mm.time_canonical_repair2 << "s" << endl;
                //cout << "mm.time_update_old_IDs_in_labels: " << mm.time_update_old_IDs_in_labels << "s" << endl;
            }
        }
        catch (string s)
        {
            cout << s << endl;
            if (use_PLL)
                dgraph_clear_global_values_PLL();
            else
                dgraph_clear_global_values_PSL();

            continue;
        }
        auto end = std::chrono::high_resolution_clock::now();
        double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
        avg_index_time = avg_index_time + runningtime / iteration_graph_times;

        /*debug*/
        if (0)
        {

        }

        dgraph_v1_check_correctness(mm, instance_graph, iteration_source_times,
            iteration_terminal_times);

        long long int index_size = 0;
        for (auto it = mm.L_in.begin(); it != mm.L_in.end(); it++)
        {
            index_size = index_size + (*it).size();
        }

        for (auto it = mm.L_out.begin(); it != mm.L_out.end(); it++)
        {
            index_size = index_size + (*it).size();
        }

        avg_index_size_per_v = avg_index_size_per_v + (double)index_size / V / iteration_graph_times;

        mm.clear_labels();
    }




    cout << "avg_index_time = " << avg_index_time << "s" << endl;



}