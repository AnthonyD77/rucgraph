#pragma once
#include <build_in_progress/HL/Hop/graph_v_of_v_idealID_HB_two_hop_labels_v1.h>

double graph_v_of_v_idealID_two_hop_v1_extract_distance_no_reduc_for_canonical_repair_with_value_M(int u, int v, int hop_cst, double value_M)
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

double graph_v_of_v_idealID_two_hop_v1_extract_distance_no_reduc_for_canonical_repair(int u, int v, int hop_cst)
{
    /* return query(u,v,h) using L_{>=r(v)}[u] and L[v] */
    /*return std::numeric_limits<double>::max() is not connected*/

    if (u == v)
    {
        return 0;
    }

    double distance = std::numeric_limits<double>::max();
    vector<two_hop_label_v1>::iterator vector1_check_pointer, vector2_check_pointer, vector1_check_pointer_end, vector2_check_pointer_end;

    vector1_check_pointer = L_temp_599[u].begin();
    vector1_check_pointer_end = L_temp_599[u].end();
    vector2_check_pointer_end = L_temp_599[v].end();

    while (vector1_check_pointer != vector1_check_pointer_end)
    {
        if (vector1_check_pointer->vertex >= v && vector1_check_pointer->hop > hop_cst)
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
            incremental_label_vectors_599[u].push_back(*it);
            continue;
        }
        
        double v_distance = it->distance - int(it->distance / value_M) * value_M;
        /* query in the canonical repair has nothing to do with reduction */
        double query_dis = graph_v_of_v_idealID_two_hop_v1_extract_distance_no_reduc_for_canonical_repair_with_value_M(u, v, it->hop, value_M);
        // cout << "check " << v << " in L[" << u << "],  dist = " << it->distance << ",  query = " << query_dis << endl;
        if (query_dis + 1e-5 >= v_distance)
        {                                                  
            incremental_label_vectors_599[u].push_back(*it);
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
            incremental_label_vectors_599[u].push_back(*it);
            continue;
        }
        
        /* query in the canonical repair has nothing to do with reduction */
        double query_dis = graph_v_of_v_idealID_two_hop_v1_extract_distance_no_reduc_for_canonical_repair(u, v, it->hop);
        if (query_dis >= it->distance)
        {                                                  
            incremental_label_vectors_599[u].push_back(*it);
        }
    }
}

void canonical_repair_element2(int target_v)
{
    L_temp_599[target_v] = incremental_label_vectors_599[target_v];
    vector<two_hop_label_v1>(L_temp_599[target_v]).swap(L_temp_599[target_v]);
}

void canonical_repair_multi_threads_with_value_M(long long int &label_size_before_canonical_repair, long long int &label_size_after_canonical_repair, double &canonical_repair_remove_label_ratio, int num_of_threads, double value_M)
{

    int max_N_ID = L_temp_599.size();
    incremental_label_vectors_599.resize(max_N_ID);

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
        int new_size = incremental_label_vectors_599[target_v].size();
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
    incremental_label_vectors_599.resize(max_N_ID);

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
        int new_size = incremental_label_vectors_599[target_v].size();
        label_size_before_canonical_repair = label_size_before_canonical_repair + old_size;
        label_size_after_canonical_repair = label_size_after_canonical_repair + new_size;
        if (new_size < old_size)
        {
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
