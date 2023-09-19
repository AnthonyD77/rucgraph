#pragma once
#include <map>
#include <shared_mutex>
#include <tool_functions/ThreadPool.h>
#include <graph_v_of_v_idealID/graph_v_of_v_idealID.h>

using namespace std;

/*PLL label format*/
class two_hop_label_v1
{
  public:
    int vertex, parent_vertex;
    int hop;
    double distance;
};

/*
    global values that should be cleared after usig PLL or PSL
    unique code for this file: 599
*/
string reach_limit_error_string_MB = "reach limit error MB";
string reach_limit_error_string_time = "reach limit error time";
long long int max_labal_size_599;
long long int labal_size_599;
int max_N_ID_for_mtx_599 = 1e7;
double max_run_time_nanoseconds_599;
auto begin_time_599 = std::chrono::high_resolution_clock::now();
graph_v_of_v_idealID ideal_graph_599;
map<pair<int, int>, int> new_edges_with_middle_v_599;
map<pair<int, int>, double> new_edges_with_origin_ec_599;
vector<vector<two_hop_label_v1>> L_temp_599;
vector<std::shared_mutex> mtx_599(max_N_ID_for_mtx_599);
queue<int> Qid_599;
vector<vector<pair<double, int>>> P_dij_599;
vector<vector<pair<double, int>>> T_dij_599;
vector<vector<two_hop_label_v1>> incremental_label_vectors_599;
vector<int> reduction_measures_2019R2;

void graph_v_of_v_idealID_two_hop_clear_global_values()
{
    vector<vector<two_hop_label_v1>>().swap(L_temp_599);
    ideal_graph_599.clear();
    queue<int>().swap(Qid_599);
    vector<vector<pair<double, int>>>().swap(P_dij_599);
    vector<vector<pair<double, int>>>().swap(T_dij_599);
    vector<int>().swap(reduction_measures_2019R2);
    vector<vector<two_hop_label_v1>>().swap(incremental_label_vectors_599);
    map<pair<int, int>, int>().swap(new_edges_with_middle_v_599);
    map<pair<int, int>, double>().swap(new_edges_with_origin_ec_599);
}

/* global querying values, used in the query func */
map<int, vector<pair<int, double>>>  R2_reduced_vertices;

void graph_v_of_v_idealID_two_hop_clear_global_values2()
{
    map<int, vector<pair<int, double>>>().swap(R2_reduced_vertices);
}

class graph_v_of_v_idealID_two_hop_case_info_v1
{
  public:
    /*hop bounded*/
    int upper_k = 0;
    double value_M = 0;
    bool use_M = 0;
    bool print_label_before_canonical_fix = 0;

    /*use reduction info*/
    bool use_2019R2 = false;
    bool use_enhanced2019R2 = false;
    bool use_non_adj_reduc_degree = false;
    int max_degree_MG_enhanced2019R2 = 100;
    int MG_num = 0;

    /*running time records*/
    double time_initialization = 0;
    double time_reduction = 0;
    double time_generate_labels = 0;
    double time_update_predecessors = 0;
    double time_canonical_repair = 0;
    double time_update_labels = 0;

    /*running limits*/
    long long int max_labal_size = 1e12;
    double max_run_time_seconds = 1e12;

    /*labels*/
    vector<int> reduction_measures_2019R2;
    vector<vector<two_hop_label_v1>> L;

    /*canonical_repair info*/
    bool use_canonical_repair = false;
    long long int label_size_before_canonical_repair = 0;
    long long int label_size_after_canonical_repair = 0;
    double canonical_repair_remove_label_ratio = 0;

    /*compute label size; this should equal label_size_after_canonical_repair when use_canonical_repair==true*/
    long long int compute_label_bit_size()
    {
        long long int size = 0;
        size = size + reduction_measures_2019R2.size() * 4;
        for (auto it = L.begin(); it != L.end(); it++)
        {
            size = size + (*it).size() * sizeof(two_hop_label_v1);
        }
        return size;
    }

    /*clear labels*/
    void clear_labels()
    {
        vector<int>().swap(reduction_measures_2019R2);
        vector<vector<two_hop_label_v1>>().swap(L);
    }

    /*printing*/
    void print_L()
    {
        cout << "print_L:" << endl;
        for (int i = 0; i < L.size(); i++)
        {
            cout << "L[" << i << "]=";
            for (int j = 0; j < L[i].size(); j++)
            {
                cout << "{" << L[i][j].vertex << "," << L[i][j].distance << "," << L[i][j].parent_vertex << "," << L[i][j].hop << "}";
            }
            cout << endl;
        }
    }
    void print_reduction_measures_2019R2()
    {
        cout << "print_reduction_measures_2019R2:" << endl;
        for (int i = 0; i < reduction_measures_2019R2.size(); i++)
        {
            cout << "reduction_measures_2019R2[" << i << "]=" << reduction_measures_2019R2[i] << endl;
        }
    }
};

bool compare_pair_second_large_to_small(const pair<int, int> &i, pair<int, int> &j)
{
    return i.second > j.second;
}

bool compare_two_hop_label_small_to_large(two_hop_label_v1 &i, two_hop_label_v1 &j)
{
    if (i.vertex != j.vertex) {
        return i.vertex < j.vertex;
    } else {
        return i.hop < j.hop;
    }
}