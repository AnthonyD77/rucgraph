#pragma once
#include <iostream>
#include <tool_functions/ThreadPool.h>
#include <shared_mutex>
#include <chrono>
#include <boost/heap/fibonacci_heap.hpp>
#include <dgraph_v_of_v/dgraph_v_of_v.h>
#include <build_in_progress/HL/dgraph/dgraph_dijkstra.h>

/*
    to include this file, please copy
    #include <dgraph_v_of_v/dgraph_PLL.h>
*/

void dgraph_PLL(dgraph_v_of_v<two_hop_weight_type> &input_graph, int N, int num_of_threads, dgraph_case_info_v1 &case_info)
{
    ideal_dgraph = input_graph;
    revse_dgraph = input_graph.reverse_graph();

    L_temp_in.resize(N);
    L_temp_out.resize(N);

    ThreadPool pool(num_of_threads);
    std::vector<std::future<int>> results;

    dij_dist.resize(num_of_threads);
    L_vk_tmp.resize(num_of_threads);

    for (int i = 0; i < num_of_threads; i++)
    {
        dij_dist[i].resize(N);
        L_vk_tmp[i].resize(N);
        for (int j = 0; j < N; j++)
        {
            dij_dist[i][j] = std::numeric_limits<two_hop_weight_type>::max();
            L_vk_tmp[i][j] = std::numeric_limits<two_hop_weight_type>::max();
        }
        Qid_595.push(i);
    }
    
    for (int v_k = 0; v_k < N; v_k++)
    {   
        results.emplace_back(pool.enqueue([v_k, N] {

                dgraph_pruned_dijkstra_v1(v_k, N, 0);
                dgraph_pruned_dijkstra_v1(v_k, N, 1);

                //dgraph_pruned_dijkstra(v_k, N, 0);
                //dgraph_pruned_dijkstra(v_k, N, 1);
                return 1;
        }));
    }

    for (auto &&result : results)
        result.get();

    results.clear();

    for (int i = 0; i < N; i++)
    {
        sort(L_temp_in[i].begin(), L_temp_in[i].end(),
             compare_two_hop_label_vertex_small_to_large); // sort is necessary
        sort(L_temp_out[i].begin(), L_temp_out[i].end(),
             compare_two_hop_label_vertex_small_to_large); // sort is necessary
    }

    case_info.L_in = L_temp_in;
    case_info.L_out = L_temp_out;
    
    //label_output_to_file("test_dgraph_label.txt", L_temp_in, L_temp_out);

    if (case_info.use_canonical_repair)
    {
        case_info.label_size_before_canonical_repair = case_info.compute_label_bit_size();
        canonical_repair_multi_threads(num_of_threads);
        case_info.L_in = L_temp_in;
        case_info.L_out = L_temp_out;
        case_info.label_size_after_canonical_repair = case_info.compute_label_bit_size();
    }

    label_output_to_file("test_dgraph_label.txt", L_temp_in, L_temp_out);
    dgraph_clear_global_values_PLL();
   
    if (case_info.use_canonical_repair)
    {
        dgraph_clear_global_values_in_canonical_repair();
    }
}
