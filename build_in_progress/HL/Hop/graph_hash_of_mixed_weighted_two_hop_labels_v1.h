#pragma once
#include <iostream>
#include <mutex>
#include <unordered_map>
#include <vector>
#include <map>
#include <tool_functions/ThreadPool.h>
#include <shared_mutex>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
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
auto begin_time_599 = std::chrono::high_resolution_clock::now();
double max_run_time_nanoseconds_599;
bool this_parallel_PLL_PSL_is_running_599 = false;
bool if_continue_599;
graph_v_of_v_idealID ideal_graph_599;
map<pair<int, int>, int> new_edges_with_middle_v;
map<pair<int, int>, double> new_edges_with_origin_ec;
vector<vector<two_hop_label_v1>> L_599;
vector<vector<two_hop_label_v1>> L_temp_599;
std::shared_mutex mtx_599_1, mtx_599_2;
int max_N_ID_for_mtx_599 = 1e7;               
vector<std::shared_mutex> mtx_599(max_N_ID_for_mtx_599);
queue<int> Qid_599;
vector<vector<pair<double, int>>> P_dij_599;
vector<vector<pair<double, int>>> T_dij_599;
vector<vector<pair<double, int>>> T_bfs_599;
vector<vector<int>> pre_599;
vector<int> vertexID_new_to_old_599;
vector<int> pos_599;
vector<int> pos_2_599;
vector<int> increment_599;
vector<vector<bool>> dirty_tag_599;
vector<int> reduction_measures_2019R2;
vector<int> f_2019R1_new_ID;
vector<vector<two_hop_label_v1>> incremental_label_vectors;
vector<pair<int, double>> min_adjs_new_IDs;

void graph_hash_of_mixed_weighted_two_hop_clear_global_values()
{
    this_parallel_PLL_PSL_is_running_599 = false;
    vector<vector<two_hop_label_v1>>().swap(L_599);
    vector<vector<two_hop_label_v1>>().swap(L_temp_599);
    ideal_graph_599.clear();
    queue<int>().swap(Qid_599);
    vector<vector<pair<double, int>>>().swap(P_dij_599);
    vector<vector<pair<double, int>>>().swap(T_dij_599);
    vector<vector<pair<double, int>>>().swap(T_bfs_599);
    vector<vector<int>>().swap(pre_599);
    vector<int>().swap(vertexID_new_to_old_599);
    vector<int>().swap(pos_599);
    vector<int>().swap(pos_2_599);
    vector<int>().swap(increment_599);
    vector<vector<bool>>().swap(dirty_tag_599);
    vector<int>().swap(reduction_measures_2019R2);
    vector<int>().swap(f_2019R1_new_ID);
    vector<vector<two_hop_label_v1>>().swap(incremental_label_vectors);
    vector<pair<int, double>>().swap(min_adjs_new_IDs);
    map<pair<int, int>, int>().swap(new_edges_with_middle_v);
    map<pair<int, int>, double>().swap(new_edges_with_origin_ec);
}

/* global querying values, used in the query func */
map<int, vector<pair<int, double>>>  R2_reduced_vertices;

void graph_hash_of_mixed_weighted_two_hop_clear_global_values2()
{
    map<int, vector<pair<int, double>>>().swap(R2_reduced_vertices);
}


class graph_hash_of_mixed_weighted_two_hop_case_info_v1
{
  public:
    /*hop bounded*/
    int upper_k = 0;
    double value_M = 0;
    bool use_M = 0;
    bool use_hbdij = 1;
    bool print_label_before_canonical_fix = 0;

    /*use reduction info*/
    bool use_2019R1 = false;
    bool use_2019R2 = false;
    bool use_enhanced2019R2 = false;
    bool use_non_adj_reduc_degree = false;
    bool use_dummy_dij_search_in_PLL = false;
    int max_degree_MG_enhanced2019R2 = 100;
    int reduce_V_num_2019R1 = 0, MG_num = 0;

    /*running time records*/
    double time_2019R1 = 0;
    double time_2019R2_or_enhanced_pre = 0;
    double time_2019R2_or_enhanced_fixlabels = 0; // s
    double time_initialization = 0;
    double time_generate_labels = 0;
    double time_canonical_repair1 = 0;
    double time_canonical_repair2 = 0;
    double time_update_old_IDs_in_labels = 0;

    /*running limits*/
    long long int max_labal_size = 1e12; // 2-hop-label num
    double max_run_time_seconds = 1e12;  // s

    /*labels*/
    vector<int> reduction_measures_2019R2; // for 2019 R2
    vector<int> reduction_measures_2019R1; // for 2019 R1;  11 means equivalent_1 relation (no edge between), 12 means equivalent_2 relation (edge between)
    vector<int> f_2019R1;                  // for 2019 R1
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
        size = size + reduction_measures_2019R1.size() * 4;
        size = size + f_2019R1.size() * 4;
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
        vector<int>().swap(reduction_measures_2019R1);
        vector<int>().swap(f_2019R1);
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
    void print_reduction_measures_2019R1()
    {
        cout << "print_reduction_measures_2019R1:" << endl;
        for (int i = 0; i < reduction_measures_2019R1.size(); i++)
        {
            cout << "reduction_measures_2019R1[" << i << "]=" << reduction_measures_2019R1[i] << endl;
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
    void print_f_2019R1()
    {
        cout << "print_f_2019R1:" << endl;
        for (int i = 0; i < f_2019R1.size(); i++)
        {
            cout << "f_2019R1[" << i << "]=" << f_2019R1[i] << endl;
        }
    }
};

/*common functions shared by PLL and PSL*/
bool compare_pair_second_large_to_small(const pair<int, int> &i, pair<int, int> &j)
{
    /*< is nearly 10 times slower than >*/
    return i.second > j.second; // < is from small to big; > is from big to small.  sort by the second item of pair<int, int>
}

bool compare_two_hop_label_small_to_large(two_hop_label_v1 &i, two_hop_label_v1 &j)
{
    return i.vertex < j.vertex; // < is from small to big; > is from big to small
}

/*the following locks are used in PLL search process and canonical_repair*/

/*
    canonical repair
*/

double graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc_for_canonical_repair_with_value_M(int u, int v, int hop_cst, double value_M)
{
    /* return query(u,v,h) using L_{>r(v)}[u] and L[v] */
    /*return std::numeric_limits<double>::max() is not connected*/

    if (u == v)
    {
        return 0;
    }

    double distance = std::numeric_limits<double>::max(); // if disconnected, return this large value

    vector<two_hop_label_v1>::iterator vector1_check_pointer, vector2_check_pointer, vector1_check_pointer_end, vector2_check_pointer_end;

    vector1_check_pointer = L_temp_599[u].begin();
    vector1_check_pointer_end = L_temp_599[u].end();
    vector2_check_pointer_end = L_temp_599[v].end();

    while (vector1_check_pointer != vector1_check_pointer_end)
    {
        if (vector1_check_pointer->vertex > v)
            break;
        vector2_check_pointer = L_temp_599[v].begin();
        while (vector2_check_pointer != vector2_check_pointer_end)
        {
            if (vector2_check_pointer->vertex > vector1_check_pointer->vertex)
                break;
            if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
            {
                int tmp_hop = int(vector1_check_pointer->distance / value_M) + int(vector2_check_pointer->distance / value_M);
                if ( tmp_hop <= hop_cst)
                {
                    double dis = vector1_check_pointer->distance + vector2_check_pointer->distance - value_M * tmp_hop;
                    if (distance > dis)
                    {
                        distance = dis;
                    }
                }
            }
            vector2_check_pointer++;
        }
        vector1_check_pointer++;
    }

    return distance;
}

double graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc_for_canonical_repair(int u, int v, int hop_cst)
{
    /* return query(u,v,h) using L_{>r(v)}[u] and L[v] */
    /*return std::numeric_limits<double>::max() is not connected*/

    if (u == v)
    {
        return 0;
    }

    double distance = std::numeric_limits<double>::max(); // if disconnected, return this large value

    vector<two_hop_label_v1>::iterator vector1_check_pointer, vector2_check_pointer, vector1_check_pointer_end, vector2_check_pointer_end;

    vector1_check_pointer = L_temp_599[u].begin();
    vector1_check_pointer_end = L_temp_599[u].end();
    vector2_check_pointer_end = L_temp_599[v].end();

    while (vector1_check_pointer != vector1_check_pointer_end)
    {
        if (vector1_check_pointer->vertex > v)
            break;
        vector2_check_pointer = L_temp_599[v].begin();
        while (vector2_check_pointer != vector2_check_pointer_end)
        {
            if (vector2_check_pointer->vertex > vector1_check_pointer->vertex)
                break;
            if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
            {
                if (vector1_check_pointer->hop + vector2_check_pointer->hop <= hop_cst)
                {
                    double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
                    if (distance > dis)
                    {
                        distance = dis;
                    }
                }
            }
            vector2_check_pointer++;
        }
        vector1_check_pointer++;
    }

    return distance;
}

void canonical_repair_element1_with_value_M(int u, double value_M)
{
    auto it = L_temp_599[u].begin(), end = L_temp_599[u].end();
    for (; it != end; it++)
    {
        int v = it->vertex;
        if (v == u)
        {
            incremental_label_vectors[u].push_back(*it);
            continue;
        }
        
        double v_distance = it->distance - int(it->distance / value_M) * value_M;
        /* query in the canonical repair has nothing to do with reduction */
        double query_dis = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc_for_canonical_repair_with_value_M(u, v, it->hop, value_M);
        // cout << "check " << v << " in L[" << u << "],  dist = " << it->distance << ",  query = " << query_dis << endl;
        if (query_dis + 1e-5 >= v_distance)
        {                                                  
            incremental_label_vectors[u].push_back(*it);
        }
    }
}

void canonical_repair_element1(int u)
{

    auto it = L_temp_599[u].begin(), end = L_temp_599[u].end();
    for (; it != end; it++)
    {
        int v = it->vertex;
        if (v == u)
        {
            incremental_label_vectors[u].push_back(*it);
            continue;
        }
        
        /* query in the canonical repair has nothing to do with reduction */
        double query_dis = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc_for_canonical_repair(u, v, it->hop);
        // cout << "check " << v << " in L[" << u << "],  dist = " << it->distance << ",  query = " << query_dis << endl;
        if (query_dis + 1e-5 >= it->distance)
        {                                                  
            incremental_label_vectors[u].push_back(*it);
        }
    }
}

void canonical_repair_element2(int target_v)
{
    L_temp_599[target_v] = incremental_label_vectors[target_v];
    vector<two_hop_label_v1>(L_temp_599[target_v]).swap(L_temp_599[target_v]);
}

void canonical_repair_multi_threads_with_value_M(long long int &label_size_before_canonical_repair, long long int &label_size_after_canonical_repair, double &canonical_repair_remove_label_ratio, int num_of_threads, double value_M)
{

    int max_N_ID = L_temp_599.size();
    incremental_label_vectors.resize(max_N_ID);

    ThreadPool pool(num_of_threads);
    std::vector<std::future<int>> results; // return typename: xxx

    /*find labels_to_be_removed*/
    for (int target_v = 0; target_v < max_N_ID; target_v++)
    {
        int size = L_temp_599[target_v].size();
        if (size > 0)
        {
            results.emplace_back(
                pool.enqueue([target_v, value_M] { // pass const type value j to thread; [] can be empty
                    canonical_repair_element1_with_value_M(target_v, value_M);
                    return 1;
                }));
        }
    }
    for (auto &&result : results)
        result.get(); // all threads finish here
    results.clear();

    /*remove labels_to_be_removed*/
    label_size_before_canonical_repair = 0;
    label_size_after_canonical_repair = 0;
    for (int target_v = 0; target_v < max_N_ID; target_v++)
    {
        int old_size = L_temp_599[target_v].size();
        int new_size = incremental_label_vectors[target_v].size();
        label_size_before_canonical_repair = label_size_before_canonical_repair + old_size;
        label_size_after_canonical_repair = label_size_after_canonical_repair + new_size;
        if (new_size < old_size)
        {
            // canonical_repair_element2(target_v);
            results.emplace_back(
                pool.enqueue([target_v] {
                    canonical_repair_element2(target_v);
                    return 1;
                }));
        }
    }
    for (auto &&result : results)
        result.get(); // all threads finish here

    canonical_repair_remove_label_ratio = (double)(label_size_before_canonical_repair - label_size_after_canonical_repair) / label_size_before_canonical_repair;
}

void canonical_repair_multi_threads(long long int &label_size_before_canonical_repair, long long int &label_size_after_canonical_repair, double &canonical_repair_remove_label_ratio, int num_of_threads)
{

    int max_N_ID = L_temp_599.size();
    incremental_label_vectors.resize(max_N_ID);

    ThreadPool pool(num_of_threads);
    std::vector<std::future<int>> results; // return typename: xxx

    /*find labels_to_be_removed*/
    for (int target_v = 0; target_v < max_N_ID; target_v++)
    {
        int size = L_temp_599[target_v].size();
        if (size > 0)
        {
            results.emplace_back(
                pool.enqueue([target_v] { // pass const type value j to thread; [] can be empty
                    canonical_repair_element1(target_v);
                    return 1; // return to results; the return type must be the same with results
                }));
        }
    }
    for (auto &&result : results)
        result.get(); // all threads finish here
    results.clear();

    /*remove labels_to_be_removed*/
    label_size_before_canonical_repair = 0;
    label_size_after_canonical_repair = 0;
    for (int target_v = 0; target_v < max_N_ID; target_v++)
    {
        int old_size = L_temp_599[target_v].size();
        int new_size = incremental_label_vectors[target_v].size();
        label_size_before_canonical_repair = label_size_before_canonical_repair + old_size;
        label_size_after_canonical_repair = label_size_after_canonical_repair + new_size;
        if (new_size < old_size)
        {
            // canonical_repair_element2(target_v);
            results.emplace_back(
                pool.enqueue([target_v] {
                    canonical_repair_element2(target_v);
                    return 1;
                }));
        }
    }
    for (auto &&result : results)
        result.get(); // all threads finish here

    canonical_repair_remove_label_ratio = (double)(label_size_before_canonical_repair - label_size_after_canonical_repair) / label_size_before_canonical_repair;
}

/*
    codes for querying distances for hop-bounded
*/

double graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(vector<vector<two_hop_label_v1>> &L, int source, int terminal, int hop_cst, double value_M)
{

    /*return std::numeric_limits<double>::max() is not connected*/

    if (source == terminal)
    {
        return 0;
    }

    double distance = std::numeric_limits<double>::max(); // if disconnected, return this large value
    auto vector1_check_pointer = L[source].begin();
    auto vector2_check_pointer = L[terminal].begin();
    auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
    while (vector1_check_pointer != pointer_L_s_end)
    {
        vector2_check_pointer = L[terminal].begin();
        while (vector2_check_pointer != pointer_L_t_end)
        {
            if (vector2_check_pointer->vertex > vector1_check_pointer->vertex)
                break;
            if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
            {
                if (value_M != 0)
                {
                    int tmp_hop = int(vector1_check_pointer->distance / value_M) + int(vector2_check_pointer->distance / value_M);
                    if ( tmp_hop <= hop_cst)
                    {
                        double dis = vector1_check_pointer->distance + vector2_check_pointer->distance - value_M * tmp_hop;
                        if (distance > dis)
                        {
                            distance = dis;
                        }
                    }
                } else {
                    // cout << "hop_cst: " << hop_cst << endl;
                    // cout << vector1_check_pointer->hop << "--" << vector2_check_pointer->hop << endl;
                    if (vector1_check_pointer->hop + vector2_check_pointer->hop <= hop_cst)
                    {
                        double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
                        // cout << "dis: " << dis << endl;
                        if (distance > dis)
                        {
                            distance = dis;
                        }
                    }
                }
            }
            vector2_check_pointer++;
        }
        vector1_check_pointer++;
    }

    return distance;
}

double graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, int source, int terminal, int hop_cst, double value_M = 0)
{
	if (source == terminal) {
		return 0;
	}

	double min_selected_distance = std::numeric_limits<double>::max();
   

	if (reduction_measures_2019R2[source] == 2)
	{
		if (reduction_measures_2019R2[terminal] == 2)
		{
			/*"Both reduced"*/
            // cout << "case 1" << endl;
            auto s_adj_begin = R2_reduced_vertices[source].begin();
            auto s_adj_end = R2_reduced_vertices[source].end();
            auto t_adj_begin = R2_reduced_vertices[terminal].begin();
            auto t_adj_end = R2_reduced_vertices[terminal].end();
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
				{
                    double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, it2->first, hop_cst - 2, value_M);
                    if (x == std::numeric_limits<double>::max()) {
                        continue;
                    }
                    else {
                        min_selected_distance = min(min_selected_distance, x + double(it1->second) + double(it2->second));
                    }
				}
			}
		}
		else
		{
			/*"Only source reduced"*/
            // cout << "case 2" << endl;
            auto s_adj_begin = R2_reduced_vertices[source].begin();
            auto s_adj_end = R2_reduced_vertices[source].end();
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
                double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, terminal, hop_cst - 1, value_M);
                // cout << "dist from " << terminal << " to " << it1->first << ": " << x << " with hop_cst " << hop_cst - 1 << endl;
                if (x == std::numeric_limits<double>::max()) {
                    continue;
                }
                else {
                    min_selected_distance = min(min_selected_distance, x + double(it1->second));
                }
			}
		}
	}
	else
	{
		if (reduction_measures_2019R2[terminal] == 2)
		{
			/*"Only terminal reduced"*/
            // cout << "case 3" << endl;
            auto t_adj_begin = R2_reduced_vertices[terminal].begin();
            auto t_adj_end = R2_reduced_vertices[terminal].end();
			for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
			{
                double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, source, it2->first, hop_cst - 1, value_M);
                // cout << "neighber " << it2->first << " with distance " << x << endl;
                if (x == std::numeric_limits<double>::max()) {
                    continue;
                }
                else {
                    min_selected_distance = min(min_selected_distance, x + double(it2->second));
                }
			}
		}
		else
		{
			/*"Nothing happened"*/
            // cout << "case 4" << endl;
			min_selected_distance = min( min_selected_distance, graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, source, terminal, hop_cst, value_M) );
		}
	}

	return min_selected_distance;
}

/*
    codes for querying paths for hop-bounded
*/
vector<pair<int, int>> graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, int source, int terminal, int hop_cst, double value_M = 0)
{

	vector<pair<int, int>> paths;
	if (source == terminal)
	{
		return paths;
	}

	double min_dis = std::numeric_limits<double>::max();
	vector<pair<int, int>> partial_edges(2);

	if (reduction_measures_2019R2[source] == 2)
	{
        auto s_adj_begin = R2_reduced_vertices[source].begin();
        auto s_adj_end = R2_reduced_vertices[source].end();
        /*"Both reduced"*/
		if (reduction_measures_2019R2[terminal] == 2)
		{
            auto t_adj_begin = R2_reduced_vertices[terminal].begin();
            auto t_adj_end = R2_reduced_vertices[terminal].end();
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
				{
                    double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, it2->first, hop_cst - 2, value_M) + double(it1->second) + double(it2->second);
                    /*After removing the two edges, it becomes the fourth case: nothing happened*/
                    if (min_dis > x)
                    {
                        min_dis = x;
                        partial_edges[0] = { it1->first, source };
                        partial_edges[1] = { it2->first, terminal };
                    }
				}
			}

			if (min_dis == std::numeric_limits<double>::max())
			{   /* disconnected */
				return paths;
			}
			paths.push_back(partial_edges[0]);
			paths.push_back(partial_edges[1]);

			/* turn to case 4 */
			vector<pair<int, int>> new_edges;
            if (value_M != 0)
                new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, partial_edges[0].first, partial_edges[1].first, hop_cst - 2, value_M);
            else
                new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, partial_edges[0].first, partial_edges[1].first, hop_cst - 2);
			if (new_edges.size() > 0)
			{
				for (int i = new_edges.size() - 1; i >= 0; i--)
				{
					paths.push_back(new_edges[i]);
				}
			}
		}
		/*"Only source reduced"*/
		else
		{
            // cout << "????????????????" << endl;
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
                double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, terminal, hop_cst - 1, value_M) + double(it1->second) - value_M;
                // cout << "neighber " << it1->first << ", dist = " << x << endl;
                // cout << "x: " << x << endl;
                if (min_dis > x)
                {
                    min_dis = x;
                    partial_edges[0] = { it1->first, source };
                }
			}
			if (min_dis == std::numeric_limits<double>::max())
			{
				return paths;
			}
            // cout << "push edge " << partial_edges[0].first << "," << partial_edges[0].second << endl;
			paths.push_back(partial_edges[0]);

            // cout << "hop_cst: " << hop_cst << endl;

            /* turn to case 4 */
			vector<pair<int, int>> new_edges;
            if (value_M != 0)
                new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, partial_edges[0].first, terminal, hop_cst - 1, value_M);
            else
                new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, partial_edges[0].first, terminal, hop_cst - 1);
			if (new_edges.size() > 0)
			{
				for (int i = new_edges.size() - 1; i >= 0; i--)
				{
					paths.push_back(new_edges[i]);
				}
			}
		}
	}
	else
	{
        /*"Only terminal reduced"*/
		if (reduction_measures_2019R2[terminal] == 2)
		{
            auto t_adj_begin = R2_reduced_vertices[terminal].begin();
	        auto t_adj_end = R2_reduced_vertices[terminal].end();
			for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
			{
                double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, source, it2->first, hop_cst - 1, value_M) + double(it2->second);
                if (min_dis > x)
                {
                    min_dis = x;
                    partial_edges[0] = { it2->first, terminal };
                }
			}
			if (min_dis == std::numeric_limits<double>::max())
			{
				return paths;
			}
			paths.push_back(partial_edges[0]);

			/* turn to case 4 */
			vector<pair<int, int>> new_edges;
            if (value_M != 0)
                new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, source, partial_edges[0].first, hop_cst - 1, value_M);
            else
                new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, source, partial_edges[0].first, hop_cst - 1);
            
			if (new_edges.size() > 0)
			{
				for (int i = new_edges.size() - 1; i >= 0; i--)
				{
					paths.push_back(new_edges[i]);
				}
			}
		}
		/* Nothing happened */
		/* In this case, the problem that the removed vertices appear in the path needs to be solved */
		else
		{
			int vector1_capped_v_parent = 0, vector2_capped_v_parent = 0;
			double distance = std::numeric_limits<double>::max();
			bool connected = false;
			auto vector1_check_pointer = L[source].begin();
			auto vector2_check_pointer = L[terminal].begin();
			auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
			// while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
			
			// 	if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
			// 	{
			// 		connected = true;
			// 		double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
			// 		if (distance > dis)
			// 		{
			// 			distance = dis;
			// 			vector1_capped_v_parent = vector1_check_pointer->parent_vertex;
			// 			vector2_capped_v_parent = vector2_check_pointer->parent_vertex;
			// 		}
			// 		vector1_check_pointer++;
			// 	}
			// 	else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex)
			// 	{
			// 		vector2_check_pointer++;
			// 	}
			// 	else
			// 	{
			// 		vector1_check_pointer++;
			// 	}
			// }
            while (vector1_check_pointer != pointer_L_s_end)
            {
                vector2_check_pointer = L[terminal].begin();
                while (vector2_check_pointer != pointer_L_t_end)
                {
                    if (vector2_check_pointer->vertex > vector1_check_pointer->vertex)
                        break;
                    if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
                    {
                        if (value_M != 0)
                        {
                            int tmp_hop = int(vector1_check_pointer->distance / value_M) + int(vector2_check_pointer->distance / value_M);
                            if (tmp_hop <= hop_cst)
                            {
                                connected = true;
                                double dis = vector1_check_pointer->distance + vector2_check_pointer->distance - value_M * tmp_hop;
                                if (distance > dis)
                                {
                                    distance = dis;
                                    vector1_capped_v_parent = vector1_check_pointer->parent_vertex;
                                    vector2_capped_v_parent = vector2_check_pointer->parent_vertex;
                                }
                            }
                        } else {
                            // cout << hop_cst << endl;
                            if (vector1_check_pointer->hop + vector2_check_pointer->hop <= hop_cst)
                            {
                                connected = true;
                                double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
                                // cout << "dis: " << dis << endl;
                                if (distance > dis)
                                {
                                    distance = dis;
                                    vector1_capped_v_parent = vector1_check_pointer->parent_vertex;
                                    vector2_capped_v_parent = vector2_check_pointer->parent_vertex;
                                }
                            }
                        }
                    }
                    vector2_check_pointer++;
                }
                vector1_check_pointer++;
            }

            // cout << "dis = " << distance << endl;
            // cout << "vector1_capped_v_parent: " << vector1_capped_v_parent << endl;
            // cout << "vector2_capped_v_parent: " << vector2_capped_v_parent << endl;

			if (connected)
			{
				if (source != vector1_capped_v_parent)
				{
					paths.push_back({ source, vector1_capped_v_parent });
					source = vector1_capped_v_parent;
                    hop_cst--;
				}
				if (terminal != vector2_capped_v_parent)
				{
					paths.push_back({ terminal, vector2_capped_v_parent });
					terminal = vector2_capped_v_parent;
                    hop_cst--;
				}
			}
			else
			{
				return paths;
			}

            /* find new */
            // cout << "turn to " << source << "-" << terminal << endl;
            vector<pair<int, int>> new_edges;
            if (value_M != 0)
			    new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, source, terminal, hop_cst, value_M);
            else
                new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, source, terminal, hop_cst);

			if (new_edges.size() > 0)
			{
				for (int i = new_edges.size() - 1; i >= 0; i--)
				{
					paths.push_back(new_edges[i]);
				}
			}
		}
	}

	return paths;
}


// /*query paths (final)*/
// void graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(vector<pair<int, int>>& paths, int old_, int new_)
// {
// 	/*
// 		this function can change specific vector on paths, used for shortesd paths query
// 		take care, we assume that the 'old_' vector only appear once
// 	*/
// 	for (auto it = paths.begin(); it != paths.end(); it++)
// 	{
// 		if (it->first == old_)
// 		{
// 			it->first = new_;
// 			break;
// 		}
// 		if (it->second == old_)
// 		{
// 			it->second = new_;
// 			break;
// 		}
// 	}
// }

// vector<pair<int, int>> graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path
// (vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& reduction_measures_2019R1,
// 	vector<int>& f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
// {

// 	vector<pair<int, int>> paths;
// 	if (source == terminal)
// 	{
// 		return paths;
// 	}

// 	if (reduction_measures_2019R1[source] == 11)
// 	{
// 		/* case 2 */
// 		if (reduction_measures_2019R1[terminal] == 11)
// 		{
// 			if (f_2019R1[source] == f_2019R1[terminal])
// 			{
// 				pair<int, double> s_min_adj = min_adjs[source];
// 				paths.push_back({ source, s_min_adj.first });
// 				paths.push_back({ terminal, s_min_adj.first });
// 				return paths;
// 			}
// 			else
// 			{
// 				paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
// 				graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
// 				graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
// 				return paths;
// 			}
// 		}
// 		/* case 3 */
// 		else if (reduction_measures_2019R1[terminal] == 12)
// 		{
// 			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
// 			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
// 			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
// 			return paths;
// 		}
// 		/* case 1 */
// 		else {
// 			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], terminal);
// 			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
// 			return paths;
// 		}
// 	}
// 	else if (reduction_measures_2019R1[source] == 12)
// 	{
// 		/* case 5 -- same with case 3 */
// 		if (reduction_measures_2019R1[terminal] == 11)
// 		{
// 			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
// 			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
// 			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
// 			return paths;
// 		}
// 		/* case 6 */
// 		else if (reduction_measures_2019R1[terminal] == 12)
// 		{
// 			if (f_2019R1[source] == f_2019R1[terminal])
// 			{
// 				pair<int, double> s_min_adj = min_adjs[source];
// 				double s_t_weight = double(graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, terminal));
// 				if (double(s_min_adj.second * 2) < s_t_weight)
// 				{
// 					paths.push_back({ source, s_min_adj.first });
// 					paths.push_back({ terminal, s_min_adj.first });
// 					return paths;
// 				}
// 				else
// 				{
// 					paths.push_back({ source, terminal });
// 					return paths;
// 				}
// 			}
// 			else // f_2019R1[source] != f_2019R1[terminal]
// 			{
// 				paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
// 				graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
// 				graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
// 				return paths;
// 			}
// 		}
// 		/* case 4 */
// 		else {
// 			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], terminal);
// 			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
// 			return paths;
// 		}

// 	}
// 	else {
// 		/* case 8 -- same with case 1 */
// 		if (reduction_measures_2019R1[terminal] == 11)
// 		{
// 			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, f_2019R1[terminal]);
// 			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
// 			return paths;
// 		}
// 		/* case 9 -- same with case 4 */
// 		else if (reduction_measures_2019R1[terminal] == 12)
// 		{
// 			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, f_2019R1[terminal]);
// 			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
// 			return paths;
// 		}
// 		/* case 7 */
// 		else {
// 			return graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, terminal);
// 		}

// 	}

// }

// /*query predecessors in instance_graph*/

// pair<int, int> graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1 /* we assume that source and terminal are not reduced by 2019R1*/
// (vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
// {

// 	/* we assume that source and terminal are not reduced by 2019R1*/

// 	if (source == terminal)
// 	{
// 		return { source, terminal }; // source_predecessor and terminal_predecessor;
// 	}

// 	double min_dis = std::numeric_limits<double>::max();
// 	vector<pair<int, int>> partial_edges(2);

// 	auto s_adj_begin = adjs[source].begin();
// 	auto s_adj_end = adjs[source].end();
// 	auto t_adj_begin = adjs[terminal].begin();
// 	auto t_adj_end = adjs[terminal].end();

// 	if (reduction_measures_2019R2[source] == 2)
// 	{
// 		if (reduction_measures_2019R2[terminal] == 2)
// 		{
// 			/*"Both reduced"*/
// 			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
// 			{
// 				for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
// 				{
// 					if (f_2019R1[it1->first] == it1->first && f_2019R1[it2->first] == it2->first) {
// 						double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, it2->first) + double(it1->second) + double(it2->second);
// 						/*After removing the two edges, it becomes the fourth case: nothing happened*/
// 						if (min_dis > x)
// 						{
// 							min_dis = x;
// 							partial_edges[0] = { it1->first, source };
// 							partial_edges[1] = { it2->first, terminal };
// 						}
// 					}
// 				}
// 			}

// 			if (min_dis == std::numeric_limits<double>::max())
// 			{ // disconnected
// 				return { source, terminal };
// 			}

// 			return { partial_edges[0].first, partial_edges[1].first };
// 		}
// 		/*"Only source reduced"*/
// 		else
// 		{
// 			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
// 			{
// 				if (f_2019R1[it1->first] == it1->first) {
// 					double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, terminal) + double(it1->second);
// 					if (min_dis > x)
// 					{
// 						min_dis = x;
// 						partial_edges[0] = { it1->first, source };
// 					}
// 				}
// 			}
// 			if (min_dis == std::numeric_limits<double>::max())
// 			{ // disconnected
// 				return { source, terminal };
// 			}

// 			pair<int, int> recursive_two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, partial_edges[0].first, terminal);

// 			return { partial_edges[0].first , recursive_two_predecessors.second };
// 		}
// 	}
// 	else
// 	{
// 		/*"Only terminal reduced"*/
// 		if (reduction_measures_2019R2[terminal] == 2)
// 		{
// 			for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
// 			{
// 				if (f_2019R1[it2->first] == it2->first) {
// 					double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, source, it2->first) + double(it2->second);
// 					if (min_dis > x)
// 					{
// 						min_dis = x;
// 						partial_edges[0] = { it2->first, terminal };
// 					}
// 				}
// 			}

// 			if (min_dis == std::numeric_limits<double>::max())
// 			{ // disconnected
// 				return { source, terminal };
// 			}

// 			pair<int, int> recursive_two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, partial_edges[0].first);

// 			return { recursive_two_predecessors.first, partial_edges[0].first };
// 		}
// 		/*
// 		"Nothing happened"
// 		In this case, the problem that the removed vertices appear in the path needs to be solved
// 		*/
// 		else
// 		{
// 			int vector1_capped_v_parent = 0, vector2_capped_v_parent = 0;
// 			double distance = std::numeric_limits<double>::max(); // if disconnected, return this large value
// 			bool connected = false;
// 			auto vector1_check_pointer = L[source].begin();
// 			auto vector2_check_pointer = L[terminal].begin();
// 			auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
// 			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
// 			{
// 				if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
// 				{
// 					connected = true;
// 					double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
// 					if (distance > dis)
// 					{
// 						distance = dis;
// 						vector1_capped_v_parent = vector1_check_pointer->parent_vertex;
// 						vector2_capped_v_parent = vector2_check_pointer->parent_vertex;
// 					}
// 					vector1_check_pointer++;
// 				}
// 				else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex)
// 				{
// 					vector2_check_pointer++;
// 				}
// 				else
// 				{
// 					vector1_check_pointer++;
// 				}
// 			}

// 			if (connected)
// 			{
// 				/*the following code will not induce redundent edge, since for each vertex v_k, there is a label (v_k,0,v_k) in L[v_k];
// 				you must ascending from both source and terminal, otherwise you may not be able to extract SP */
// 				pair<int, int> recursive_two_predecessors;
// 				if (source != vector1_capped_v_parent) // ascending from source
// 				{
// 					recursive_two_predecessors.first = vector1_capped_v_parent;
// 				}
// 				else {
// 					recursive_two_predecessors.first = source;
// 				}
// 				if (terminal != vector2_capped_v_parent) // ascending from terminal
// 				{
// 					recursive_two_predecessors.second = vector2_capped_v_parent;
// 				}
// 				else {
// 					recursive_two_predecessors.second = terminal;
// 				}
// 				return recursive_two_predecessors;
// 			}
// 			else
// 			{
// 				return { source, terminal };
// 			}
// 		}
// 	}
// }

// pair<int, int> graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors
// (vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& reduction_measures_2019R1,
// 	vector<int>& f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
// {

// 	pair<int, int> two_predecessors; // source_predecessor and terminal_predecessor;
// 	/*
// 	if (source_predecessor == source && terminal_predecessor == terminal){
// 	no edge is in path
// 	}

// 	if (source_predecessor == terminal && terminal_predecessor == source){
// 	one edge is in path
// 	}

// 	*/

// 	if (source == terminal)
// 	{
// 		return { source, terminal };
// 	}

// 	if (reduction_measures_2019R1[source] == 11)
// 	{
// 		/* case 2 */
// 		if (reduction_measures_2019R1[terminal] == 11)
// 		{
// 			if (f_2019R1[source] == f_2019R1[terminal])
// 			{
// 				pair<int, double> s_min_adj = min_adjs[source];
// 				return { s_min_adj.first , s_min_adj.first };
// 			}
// 			else
// 			{
// 				two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
// 				if (two_predecessors.first == f_2019R1[source]) {
// 					two_predecessors.first = source;
// 				}
// 				if (two_predecessors.second == f_2019R1[source]) {
// 					two_predecessors.second = source;
// 				}
// 				if (two_predecessors.first == f_2019R1[terminal]) {
// 					two_predecessors.first = terminal;
// 				}
// 				if (two_predecessors.second == f_2019R1[terminal]) {
// 					two_predecessors.second = terminal;
// 				}
// 				return two_predecessors;
// 			}
// 		}
// 		/* case 3 */
// 		else if (reduction_measures_2019R1[terminal] == 12)
// 		{
// 			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
// 			if (two_predecessors.first == f_2019R1[source]) {
// 				two_predecessors.first = source;
// 			}
// 			if (two_predecessors.second == f_2019R1[source]) {
// 				two_predecessors.second = source;
// 			}
// 			if (two_predecessors.first == f_2019R1[terminal]) {
// 				two_predecessors.first = terminal;
// 			}
// 			if (two_predecessors.second == f_2019R1[terminal]) {
// 				two_predecessors.second = terminal;
// 			}
// 			return two_predecessors;
// 		}
// 		/* case 1 */
// 		else {
// 			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], terminal);
// 			if (two_predecessors.first == f_2019R1[source]) {
// 				two_predecessors.first = source;
// 			}
// 			if (two_predecessors.second == f_2019R1[source]) {
// 				two_predecessors.second = source;
// 			}
// 			return two_predecessors;
// 		}
// 	}
// 	else if (reduction_measures_2019R1[source] == 12)
// 	{
// 		/* case 5 -- same with case 3 */
// 		if (reduction_measures_2019R1[terminal] == 11)
// 		{
// 			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
// 			if (two_predecessors.first == f_2019R1[source]) {
// 				two_predecessors.first = source;
// 			}
// 			if (two_predecessors.second == f_2019R1[source]) {
// 				two_predecessors.second = source;
// 			}
// 			if (two_predecessors.first == f_2019R1[terminal]) {
// 				two_predecessors.first = terminal;
// 			}
// 			if (two_predecessors.second == f_2019R1[terminal]) {
// 				two_predecessors.second = terminal;
// 			}
// 			return two_predecessors;
// 		}
// 		/* case 6 */
// 		else if (reduction_measures_2019R1[terminal] == 12)
// 		{
// 			if (f_2019R1[source] == f_2019R1[terminal])
// 			{
// 				pair<int, double> s_min_adj = min_adjs[source];
// 				double s_t_weight = double(graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, terminal));
// 				if (double(s_min_adj.second * 2) < s_t_weight)
// 				{
// 					return { s_min_adj.first , s_min_adj.first };
// 				}
// 				else
// 				{
// 					return { terminal, source };
// 				}
// 			}
// 			else // f_2019R1[source] != f_2019R1[terminal]
// 			{
// 				two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
// 				if (two_predecessors.first == f_2019R1[source]) {
// 					two_predecessors.first = source;
// 				}
// 				if (two_predecessors.second == f_2019R1[source]) {
// 					two_predecessors.second = source;
// 				}
// 				if (two_predecessors.first == f_2019R1[terminal]) {
// 					two_predecessors.first = terminal;
// 				}
// 				if (two_predecessors.second == f_2019R1[terminal]) {
// 					two_predecessors.second = terminal;
// 				}
// 				return two_predecessors;
// 			}
// 		}
// 		/* case 4 */
// 		else {
// 			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], terminal);
// 			if (two_predecessors.first == f_2019R1[source]) {
// 				two_predecessors.first = source;
// 			}
// 			if (two_predecessors.second == f_2019R1[source]) {
// 				two_predecessors.second = source;
// 			}
// 			return two_predecessors;
// 		}
// 	}
// 	else {
// 		/* case 8 -- same with case 1 */
// 		if (reduction_measures_2019R1[terminal] == 11)
// 		{
// 			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, f_2019R1[terminal]);
// 			if (two_predecessors.first == f_2019R1[terminal]) {
// 				two_predecessors.first = terminal;
// 			}
// 			if (two_predecessors.second == f_2019R1[terminal]) {
// 				two_predecessors.second = terminal;
// 			}
// 			return two_predecessors;
// 		}
// 		/* case 9 -- same with case 4 */
// 		else if (reduction_measures_2019R1[terminal] == 12)
// 		{
// 			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, f_2019R1[terminal]);
// 			if (two_predecessors.first == f_2019R1[terminal]) {
// 				two_predecessors.first = terminal;
// 			}
// 			if (two_predecessors.second == f_2019R1[terminal]) {
// 				two_predecessors.second = terminal;
// 			}
// 			return two_predecessors;
// 		}
// 		/* case 7 */
// 		else {
// 			return graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, terminal);
// 		}

// 	}

// }