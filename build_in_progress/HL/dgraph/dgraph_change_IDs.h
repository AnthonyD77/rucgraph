#pragma once
#include <tool_functions/ThreadPool.h>
#include <dgraph_v_of_v/dgraph_v_of_v.h>
#include <cmath>

template <typename weight_type>
void dgraph_change_IDs_element(dgraph_v_of_v<weight_type>& graph, int N, vector<int>& old_to_new, ThreadPool& pool, std::vector<std::future<int>>& results) {

    std::vector<std::vector<std::pair<int, weight_type>>> New_INs(N);
    std::vector<std::vector<std::pair<int, weight_type>>> New_OUTs(N);
    auto* g = &graph;
    auto* old_to_new_p = &old_to_new;
    auto* New_INs_p = &New_INs;
    auto* New_OUTs_p = &New_OUTs;

    for (int i = 0; i < N; i++) {
        results.emplace_back(pool.enqueue([i, g, old_to_new_p, New_INs_p, New_OUTs_p] {
            int in_size = g->INs[i].size();
            std::vector<std::pair<int, weight_type>> v1 = g->INs[i];
            int out_size = g->OUTs[i].size();
            std::vector<std::pair<int, weight_type>> v2 = g->OUTs[i];
            for (int j = 0; j < in_size; j++) {
                v1[j].first = (*old_to_new_p)[v1[j].first];
            }
            (*New_INs_p)[(*old_to_new_p)[i]] = v1;

            for (int j = 0; j < out_size; j++) {
                v2[j].first = (*old_to_new_p)[v2[j].first];
            }
            (*New_OUTs_p)[(*old_to_new_p)[i]] = v2;
            return 1; }));
    }
    for (auto&& result : results)
        result.get();
    results.clear();

    graph.INs = New_INs; // O(E)
    graph.OUTs = New_OUTs; // O(E)
}


template <typename weight_type>
void dgraph_change_IDs_sum_IN_OUT_degrees(dgraph_v_of_v<weight_type>& graph, vector<int> &new_to_old, ThreadPool& pool, std::vector<std::future<int>>& results) {

    /*time complexity: O(E+V*logV)*/

    int N = graph.INs.size();
    new_to_old.resize(N);
    vector<pair<int, int>> value(N);

    for (int i = 0; i < N; i++) {
        value[i] = { graph.OUTs[i].size() + graph.INs[i].size() , i }; // increasing order of sums of in and out degrees
    }  
    sort(value.begin(), value.end()); // 升序排列 O(V*logV)

    vector<int> old_to_new(N);
    for (int j = 0; j < N; j++) {
        old_to_new[value[j].second] = N - 1 - j; // new ID 0 has the largest sum of in and out degrees
        new_to_old[N - 1 - j] = value[j].second;
    }

    dgraph_change_IDs_element(graph, N, old_to_new, pool, results);
}


template <typename weight_type>
void dgraph_change_IDs_weighted_degrees(dgraph_v_of_v<weight_type>& graph, vector<int> &new_to_old, ThreadPool& pool, std::vector<std::future<int>>& results) {

    /*time complexity: O(E+V*logV)*/

    int N = graph.INs.size();
    new_to_old.resize(N);
    vector<pair<int, int>> value(N);
    auto* g = &graph;
    auto* value_p = &value;

    for (int i = 0; i < N; i++) {
        results.emplace_back(pool.enqueue([i, g, value_p] {
            weight_type w = 0;
            for (auto x : g->OUTs[i]) {
                w += log(1 / x.second);
            }
            for (auto x : g->INs[i]) {
                w += log(1 / x.second);
            }
            (*value_p)[i] = { w , i }; // increasing order of values
            return 1; }));
    }
    for (auto&& result : results)
        result.get();
    results.clear();
    sort(value.begin(), value.end()); // 升序排列 O(V*logV)

    vector<int> old_to_new(N);
    for (int j = 0; j < N; j++) {
        old_to_new[value[j].second] = N - 1 - j; // new ID 0 has the largest value
        new_to_old[N - 1 - j] = value[j].second;
    }

    dgraph_change_IDs_element(graph, N, old_to_new, pool, results);
}
