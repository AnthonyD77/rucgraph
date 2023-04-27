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

#include <build_in_progress/HL/dgraph/test_dgraph.h>

int main()
{
    test_dgraph_PLL_PSL();
}

------------------------------------------------------------------------------------------
Commends for running the above cpp file on Linux:

g++ -std=c++17 -I/home/boost_1_75_0 -I/home/dengs/dgraph try.cpp -lpthread -Ofast -o A
./A
rm A

(optional to put the above commends in run.sh, and then use the comment: sh run.sh)
*/

#include <unordered_map>
#include <dgraph_v_of_v/dgraph_v_of_v.h>
#include <dgraph_v_of_v/dgraph_generate_random_dgraph.h>
#include <dgraph_v_of_v/dgraph_save_dgraph.h>
#include <dgraph_v_of_v/dgraph_read_dgraph.h>
#include <build_in_progress/HL/dgraph/dgraph_change_IDs.h>
#include <dgraph_v_of_v/dgraph_shortest_paths.h>
#include <build_in_progress/HL/dgraph/dgraph_two_hop_label.h>
#include <build_in_progress/HL/dgraph/dgraph_PLL.h>
#include <build_in_progress/HL/dgraph/dgraph_PSL.h>
#include <build_in_progress/HL/dgraph/dgraph_CT.h>


/*check_correctness*/

template <typename weight_type>
void dgraph_v1_check_correctness(dgraph_case_info_v1& case_info, dgraph_case_info_v2& case_info2, dgraph_v_of_v<weight_type>& instance_graph, int iteration_source_times, int iteration_terminal_times, bool use_CT) {

    /*below is for checking whether the above labels are right (by randomly computing shortest distances)*/

    boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(instance_graph.INs.size() - 1) };

    for (int yy = 0; yy < iteration_source_times; yy++) {
        int source = dist(boost_random_time_seed);
        // source = 27; //cout << "source = " << source << endl;

        auto distances = dgraph_shortest_distances_source_to_all(instance_graph, source);

        for (int xx = 0; xx < iteration_terminal_times; xx++) {
            int terminal = dist(boost_random_time_seed);
            // terminal = 4; //cout << "terminal = " << terminal << endl;

            weight_type dis;

            if (use_CT) {
                dis = CT_extract_distance(case_info2, source, terminal);
            }
            else {
                dis = dgraph_v1_extract_shortest_distance(case_info.L_in, case_info.L_out, source, terminal);
            }

            if (abs(dis - distances[terminal]) > 1e-4) {
                cout << "dis = " << dis << endl;
                cout << "distances[terminal] = " << distances[terminal] << endl;
                cout << "abs(dis - distances[terminal]) > 1e-5!" << endl;
                cout << "source = " << source << endl;
                cout << "terminal = " << terminal << endl;
                cout << "source vector:" << endl;
                for (auto it = case_info.L_out[source].begin(); it != case_info.L_out[source].end(); it++) {
                    cout << "<" << it->vertex << "," << it->distance << ">";
                }
                cout << endl;
                cout << "terminal vector:" << endl;
                for (auto it = case_info.L_in[terminal].begin(); it != case_info.L_in[terminal].end(); it++) {
                    cout << "<" << it->vertex << "," << it->distance << ">";
                }
                cout << endl;
                getchar();
            }
        }
    }
}

void test_dgraph_PLL_PSL() {
    /*parameters*/
    int iteration_graph_times = 100, iteration_source_times = 100, iteration_terminal_times = 100;
    int V = 100, E = 500, precision = 1, thread_num = 5;
    two_hop_weight_type ec_min = 1, ec_max = 1;

    double avg_index_time = 0, avg_index_bit_size = 0;

    bool use_PLL = 1; // 1: PLL 0: PSL

    dgraph_case_info_v1 mm;
    dgraph_case_info_v2 mm2;
    mm.use_canonical_repair = 1;
    mm.max_run_time_seconds = 1;
    mm.max_labal_bit_size = 4e5;

    /*iteration*/
    for (int i = 0; i < iteration_graph_times; i++) {
        cout << i << endl;

        /*input and output; below is for generating random new graph, or read saved graph*/
        int generate_new_graph = 1;

        dgraph_v_of_v<two_hop_weight_type> instance_graph;

        if (generate_new_graph == 1) {          
            instance_graph = dgraph_generate_random_dgraph(V, E, ec_min, ec_max, precision, boost_random_time_seed);
            vector<int> new2old;
            ThreadPool pool(thread_num);
            std::vector<std::future<int>> results;
            dgraph_change_IDs_sum_IN_OUT_degrees(instance_graph, new2old, pool, results); //id������õ�
            dgraph_save_dgraph("random_dgraph_test.txt", instance_graph);
        }
        else {
            dgraph_read_dgraph("random_dgraph_test.txt", instance_graph);
        }

        auto begin = std::chrono::high_resolution_clock::now();
        try {
            if (use_PLL) {
                dgraph_PLL(instance_graph, thread_num, mm);
            }
            else {
                dgraph_PSL(instance_graph, thread_num, mm);
            }
        }
        catch (string s) {
            cout << s << endl;
            dgraph_clear_global_values_PLL_PSL();
            mm.clear_labels();
            continue;
        }
        auto end = std::chrono::high_resolution_clock::now();
        double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
        avg_index_time = avg_index_time + runningtime / iteration_graph_times;

        dgraph_v1_check_correctness(mm, mm2, instance_graph, iteration_source_times, iteration_terminal_times, 0);

        avg_index_bit_size += (double)compute_label_bit_size(mm.L_out, mm.L_in) / iteration_graph_times;

        mm.clear_labels();
    }

    if (use_PLL) {
        cout << "V = " << V << " E = " << E << " thread_num = " << thread_num << " PLL avg_index_time = " << avg_index_time << "s avg_index_bit_size = " << avg_index_bit_size << "bit" << endl;
    }
    else {
        cout << "V = " << V << " E = " << E << " thread_num = " << thread_num << " PSL avg_index_time = " << avg_index_time << "s avg_index_bit_size = " << avg_index_bit_size << "bit" << endl;
    }
}

void test_dgraph_label_of_PLL_PSL_is_same_or_not()
{
    vector<int> new2old;
    for (int xx = 0; xx < 1e2; xx++) {
        int V = 1000, E = 5000, precision = 1, thread_num = 10;
        two_hop_weight_type ec_min = 0.1, ec_max = 1;
        dgraph_v_of_v<two_hop_weight_type> instance_graph;
        instance_graph = dgraph_generate_random_dgraph(V, E, ec_min, ec_max, precision, boost_random_time_seed);
        ThreadPool pool(thread_num);
        std::vector<std::future<int>> results;
        dgraph_change_IDs_sum_IN_OUT_degrees(instance_graph, new2old, pool, results); //id是排序好的

        dgraph_case_info_v1 mm;
        mm.use_canonical_repair = 1;
        dgraph_PLL(instance_graph, thread_num, mm);

        dgraph_case_info_v1 mm1;
        mm1.use_canonical_repair = 1;
        dgraph_PSL(instance_graph, thread_num, mm1);

        bool is_same = true;

        for (int i = 0; i < V; i++)
        {
            if (mm.L_in[i].size() == mm1.L_in[i].size() && mm.L_out[i].size() == mm1.L_out[i].size())
            {
                for (int j = 0; j < mm.L_in[i].size(); j++)
                {
                    if (mm.L_in[i][j].vertex != mm1.L_in[i][j].vertex || (mm.L_in[i][j].distance - mm1.L_in[i][j].distance > 1e-5))
                    {
                        is_same = false;
                        break;
                    }

                }

                for (int j = 0; j < mm.L_out[i].size(); j++)
                {

                    if (mm.L_out[i][j].vertex != mm1.L_out[i][j].vertex || (mm.L_out[i][j].distance - mm1.L_out[i][j].distance > 1e-5))
                    {
                        is_same = false;
                        break;
                    }
                }
            }
            else {
                is_same = false;
                break;
            }
        }

        cout << "is same? " << is_same << endl;
        if (!is_same) {
            getchar();
        }

        mm.clear_labels();
        mm1.clear_labels();
    }
}

void test_dgraph_CT() {
    /*parameters*/
    int iteration_graph_times = 30, iteration_source_times = 100, iteration_terminal_times = 100;

    int generate_new_graph = 1;

    int V = 1000, E = 10000, precision = 1, thread_num = 10;
    two_hop_weight_type ec_min = 0.1, ec_max = 1;
    double avg_index_time = 0, avg_index_size_per_v = 0;

    /*reduction method selection*/
    dgraph_case_info_v1 mm;
    dgraph_case_info_v2 ct_info;
    ct_info.thread_num = thread_num;
    ct_info.d = 5;
    ct_info.use_PLL = 1;
    ct_info.two_hop_order_method = 1;
    ct_info.max_bit_size = 1e7;
    ct_info.max_run_time_seconds = 0.5;

    /*iteration*/
    for (int i = 0; i < iteration_graph_times; i++)
    {
        cout << "iteration_graph_times:" << i << endl;

        dgraph_v_of_v<two_hop_weight_type> instance_graph;

        if (generate_new_graph == 1)
        {
            instance_graph = dgraph_generate_random_dgraph(V, E, ec_min, ec_max, precision, boost_random_time_seed);
            dgraph_save_dgraph("random_dgraph_test_CT.txt", instance_graph);
        }
        else
        {
            dgraph_read_dgraph("random_dgraph_test_CT.txt", instance_graph);
        }

        auto begin = std::chrono::high_resolution_clock::now();
        try {
            CT_dgraph(instance_graph, ct_info);
        }
        catch (string s)
        {
            cout << s << endl;
            clear_gloval_values_CT();
            ct_info.clear_labels();
            continue;
        }
        auto end = std::chrono::high_resolution_clock::now();
        double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
        avg_index_time += (double)runningtime / iteration_graph_times;

        if (0)
        {
            ct_info.two_hop_case_info.print_L();
            ct_info.print_Bags();
            ct_info.print_isIntree();
            ct_info.print_root();
        }
        
        dgraph_v1_check_correctness(mm, ct_info, instance_graph, iteration_source_times, iteration_terminal_times, 1);

        //ct_info.print_time();
        ct_info.clear_labels();
    }

    cout << "V = " << V << " E = " << E << " ct_info.d = " << ct_info.d << " thread_num = " << thread_num << " ct_info.use_PLL = "
        << ct_info.use_PLL << " avg_index_time = " << avg_index_time << "s" << endl;
}





/*test real data*/

#include <dgraph_v_of_v/dgraph_compare_two_dgraphs.h>
#include <dgraph_v_of_v/dgraph_binary_save_read_dgraph.h>

void dgraph_read_dgraph_from_txt(std::string instance_name, dgraph_v_of_v<two_hop_weight_type>& input_graph, int Jacard_or_Random)
{
    int vertex_num;
    input_graph.clear();
    std::string line_content;
    std::ifstream myfile(instance_name);
    if (myfile.is_open())
    {
        getline(myfile, line_content);
        int index1 = line_content.find(' ', 63);
        vertex_num = stoi(line_content.substr(63, index1 - 63));
        input_graph = dgraph_v_of_v<two_hop_weight_type>(vertex_num);
        getline(myfile, line_content);
        while (getline(myfile, line_content)) // read file line by line
        {
            std::vector<std::string> Parsed_content = parse_string(line_content, " ");
            int v1 = std::stoi(Parsed_content[0]);
            int v2 = std::stoi(Parsed_content[1]);
            two_hop_weight_type w;
            if (Jacard_or_Random) {
                w = std::stod(Parsed_content[3]);
            }
            else
                w = std::stod(Parsed_content[2]);

            input_graph.add_edge(v1, v2, w);
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

bool test_dgraph_weight(dgraph_v_of_v<two_hop_weight_type>& input_graph, int Jacard_or_random, int iteration_source_times) {

    boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(input_graph.INs.size() - 1) };

    for (int yy = 0; yy < iteration_source_times; yy++) {
        int source = dist(boost_random_time_seed);
        int source_out_size = input_graph.OUTs[source].size();
        for (int xx = 0; xx < source_out_size; xx++) {
            int terminal = input_graph.OUTs[source][xx].first;
            two_hop_weight_type dis = input_graph.OUTs[source][xx].second;

            if (Jacard_or_random) {
                int s_cap_t = 0;
                int terminal_in_size = input_graph.INs[terminal].size();
                auto s_pointer = input_graph.OUTs[source].begin();
                auto s_end = input_graph.OUTs[source].end();
                auto t_pointer = input_graph.INs[terminal].begin();
                auto t_end = input_graph.INs[terminal].end();
                while (s_pointer != s_end && t_pointer != t_end) {
                    if (s_pointer->first == t_pointer->first)
                    {
                        s_cap_t++;
                        s_pointer++;
                    }
                    else if (s_pointer->first > t_pointer->first)
                    {
                        t_pointer++;
                    }
                    else {
                        s_pointer++;
                    }
                }

                int s_cup_t = terminal_in_size + source_out_size - s_cap_t;
                two_hop_weight_type ec = 1.0 - (two_hop_weight_type)s_cap_t / (s_cup_t);

                if (ec - dis < 1e-5)
                    continue;

                return false;
            }
            else {
                if (dis >= 0 && dis <= 100)
                    continue;

                return false;
            }
        }
    }

    return true;
}

void test_real_data() {

    vector<string> datas = { "soc-Epinions1" };

    for (auto data_name : datas) {

        string path = "/home/malu/DHL_exp/" + data_name + "/"; // server
        path = ""; // local

        int iteration_source_times = 100;

        dgraph_v_of_v<two_hop_weight_type> input_txt_graph;
        dgraph_read_dgraph_from_txt(path + data_name + ".txt", input_txt_graph, 1);
        dgraph_v_of_v<two_hop_weight_type> g;
        dgraph_binary_read_dgraph(path + data_name + "_Jacard.bin", g);
        if (!dgraph_compare_two_dgraphs_not_eaxct_same_weight(g, input_txt_graph)) {
            cout << "1 input_txt_graph != g" << endl;
            getchar();
        }
        if (!test_dgraph_weight(g, 1, iteration_source_times)) {
            cout << "1 test_dgraph_weight is not right" << endl;
            getchar();
        }

        dgraph_read_dgraph_from_txt(path + data_name + ".txt", input_txt_graph, 0);
        dgraph_binary_read_dgraph(path + data_name + "_random.bin", g);
        if (!dgraph_compare_two_dgraphs_not_eaxct_same_weight(g, input_txt_graph)) {
            cout << "0 input_txt_graph != g" << endl;
            getchar();
        }
        if (!test_dgraph_weight(g, 0, iteration_source_times)) {
            cout << "0 test_dgraph_weight is not right" << endl;
            getchar();
        }
    }
}







/*compare*/

double avg_query_time(int query_times, int N, vector<vector<two_hop_label>>& L_in, vector<vector<two_hop_label>>& L_out) {

    boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(N-1) };

    double avg_query_time_ms = 0;

    for (int i = 0; i < query_times; i++) {
        int source = dist(boost_random_time_seed), terminal = dist(boost_random_time_seed);;
        auto begin = std::chrono::high_resolution_clock::now();
        dgraph_v1_extract_shortest_distance(L_in, L_out, source, terminal);
        auto end = std::chrono::high_resolution_clock::now();
        double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e6; // ms
        avg_query_time_ms+= runningtime / query_times;
    }

    return avg_query_time_ms;
}

void compare_PLL_PSL() {

    /*
    PLL is faster in sparse weighted graphs;
    PSL is faster in dense unweighted graphs;
    */

    /*parameters*/
    int iteration_graph_times = 20;
    int V = 1000, E = 10000, precision = 1, thread_num = 5;
    two_hop_weight_type ec_min = 1, ec_max = 10;

    double PLL_avg_index_time = 0, PSL_avg_index_time = 0;

    ThreadPool pool(thread_num);
    std::vector<std::future<int>> results;

    dgraph_case_info_v1 mm;
    mm.use_canonical_repair = 1;

    /*iteration*/
    for (int i = 0; i < iteration_graph_times; i++) {
        cout << i << endl;

        /*input and output; below is for generating random new graph, or read saved graph*/
        int generate_new_graph = 1;

        dgraph_v_of_v<two_hop_weight_type> instance_graph;

        if (generate_new_graph == 1) {
            instance_graph = dgraph_generate_random_dgraph(V, E, ec_min, ec_max, precision, boost_random_time_seed);
            vector<int> new2old;           
            dgraph_change_IDs_sum_IN_OUT_degrees(instance_graph, new2old, pool, results); //id������õ�
            dgraph_save_dgraph("random_dgraph_test.txt", instance_graph);
        }
        else {
            dgraph_read_dgraph("random_dgraph_test.txt", instance_graph);
        }

        double runningtime1 = 0, runningtime2 = 0;
        if (1) {
            auto begin = std::chrono::high_resolution_clock::now();
            dgraph_PLL(instance_graph, thread_num, mm);
            auto end = std::chrono::high_resolution_clock::now();
            runningtime1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
            PLL_avg_index_time += runningtime1 / iteration_graph_times;
            mm.clear_labels();
        }

        if (1) {
            auto begin = std::chrono::high_resolution_clock::now();
            dgraph_PSL(instance_graph, thread_num, mm);
            auto end = std::chrono::high_resolution_clock::now();
            runningtime2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
            PSL_avg_index_time += runningtime2 / iteration_graph_times;
            mm.clear_labels();
        }

        dgraph_case_info_v2 mm2;
        mm2.thread_num = thread_num;
        mm2.d = 0;
        CT_dgraph(instance_graph, mm2);
        cout << "PLL is faster? " << (runningtime1 < runningtime2) << " " << mm2.V_log_V_uk << endl;
    }

    cout << "V = " << V << " E = " << E << " thread_num = " << thread_num << " PLL_avg_index_time = " << PLL_avg_index_time << "s" 
        << " PSL_avg_index_time = " << PSL_avg_index_time << "s" << endl << "PLL_avg_index_time/PSL_avg_index_time = " << PLL_avg_index_time/PSL_avg_index_time << endl;
}

void compare_different_sorting_method() {
    /*parameters*/
    int iteration_graph_times = 100, query_times_per_g = 1e4;
    int V = 1000, E = 5000, precision = 1, thread_num = 10;
    two_hop_weight_type ec_min = 0.1, ec_max = 1;

    double degree_order_avg_index_time = 0, degree_order_avg_label_bit_size = 0, degree_order_avg_query_time = 0,
        weighted_degree_order_avg_index_time = 0, weighted_degree_order_avg_label_bit_size = 0, weighted_degree_order_avg_query_time = 0;

    bool use_PLL = 1; // 1: PLL 0: PSL

    dgraph_case_info_v1 mm;
    mm.use_canonical_repair = 1;

    vector<int> new2old;
    ThreadPool pool(thread_num);
    std::vector<std::future<int>> results;

    /*iteration*/
    for (int i = 0; i < iteration_graph_times; i++) {
        cout << i << endl;

        dgraph_v_of_v<two_hop_weight_type> instance_graph = dgraph_generate_random_dgraph(V, E, ec_min, ec_max, precision, boost_random_time_seed);

        /*degree order*/
        if (1) {
            dgraph_change_IDs_sum_IN_OUT_degrees(instance_graph, new2old, pool, results); 
            auto begin = std::chrono::high_resolution_clock::now();
            try {
                if (use_PLL) {
                    dgraph_PLL(instance_graph, thread_num, mm);
                }
                else {
                    dgraph_PSL(instance_graph, thread_num, mm);
                }
            }
            catch (string s) {
                cout << s << endl;
                dgraph_clear_global_values_PLL_PSL();
                continue;
            }
            auto end = std::chrono::high_resolution_clock::now();
            double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
            degree_order_avg_index_time += runningtime / iteration_graph_times;
            degree_order_avg_label_bit_size += (double)compute_label_bit_size(mm.L_in,mm.L_out) / iteration_graph_times;
            degree_order_avg_query_time += (double)avg_query_time(query_times_per_g, V, mm.L_in, mm.L_out) / iteration_graph_times;
            mm.clear_labels();
        }

        /*weighted degree order*/
        if (1) {
            dgraph_change_IDs_weighted_degrees(instance_graph, new2old, pool, results);
            auto begin = std::chrono::high_resolution_clock::now();
            try {
                if (use_PLL) {
                    dgraph_PLL(instance_graph, thread_num, mm);
                }
                else {
                    dgraph_PSL(instance_graph, thread_num, mm);
                }
            }
            catch (string s) {
                cout << s << endl;
                dgraph_clear_global_values_PLL_PSL();
                continue;
            }
            auto end = std::chrono::high_resolution_clock::now();
            double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
            weighted_degree_order_avg_index_time += runningtime / iteration_graph_times;
            weighted_degree_order_avg_label_bit_size += (double)compute_label_bit_size(mm.L_in, mm.L_out) / iteration_graph_times;
            weighted_degree_order_avg_query_time += (double)avg_query_time(query_times_per_g, V, mm.L_in, mm.L_out) / iteration_graph_times;
            mm.clear_labels();
        }

    }

    cout << "V = " << V << " E = " << E << " thread_num = " << thread_num << " use_PLL = " << use_PLL << endl;
    cout << "degree_order:  avg_index_time = " << degree_order_avg_index_time << "s" << " avg_label_bit_size = " << degree_order_avg_label_bit_size << "bit" 
        << " avg_query_time = " << degree_order_avg_query_time << "ms" << endl;
    cout << "weighted_degree_order: avg_index_time = " << weighted_degree_order_avg_index_time << "s" << " avg_label_bit_size = " << weighted_degree_order_avg_label_bit_size << "bit" 
        << " avg_query_time = " << weighted_degree_order_avg_query_time << "ms" << endl;
}





