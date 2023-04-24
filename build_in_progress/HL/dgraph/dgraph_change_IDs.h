#pragma once

#include <dgraph_v_of_v/dgraph_v_of_v.h>
#include <cmath>

template <typename weight_type>
void dgraph_change_IDs_element(dgraph_v_of_v<weight_type>& graph, int N, vector<int>& old_to_new) {

    std::vector<std::vector<std::pair<int, weight_type>>> New_INs(N);
    std::vector<std::vector<std::pair<int, weight_type>>> New_OUTs(N);
    for (int i = 0; i < N; i++) {
        int in_size = graph.INs[i].size();
        std::vector<std::pair<int, weight_type>> v1 = graph.INs[i];
        int out_size = graph.OUTs[i].size();
        std::vector<std::pair<int, weight_type>> v2 = graph.OUTs[i];
        for (int j = 0; j < in_size; j++) {
            v1[j].first = old_to_new[v1[j].first];
        }
        New_INs[old_to_new[i]] = v1;

        for (int j = 0; j < out_size; j++) {
            v2[j].first = old_to_new[v2[j].first];
        }
        New_OUTs[old_to_new[i]] = v2;
    }

    graph.INs = New_INs; // O(E)
    graph.OUTs = New_OUTs; // O(E)
}


template <typename weight_type>
void dgraph_change_IDs_sum_IN_OUT_degrees(dgraph_v_of_v<weight_type>& graph, vector<int> &new_to_old) {

    /*time complexity: O(E+V*logV)*/

    int N = graph.INs.size();
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

    dgraph_change_IDs_element(graph, N, old_to_new);
}


template <typename weight_type>
void dgraph_change_IDs_weighted_degrees(dgraph_v_of_v<weight_type>& graph, vector<int> &new_to_old) {

    /*time complexity: O(E+V*logV)*/

    int N = graph.INs.size();
    vector<pair<int, int>> value(N);

    for (int i = 0; i < N; i++) {
        weight_type w = 0;
        for (auto x : graph.OUTs[i]) {         
            w += log(1 / x.second);
        }
        for (auto x : graph.INs[i]) {
            w += log(1 / x.second);
        }
        value[i] = { w , i }; // increasing order of values
    }
    sort(value.begin(), value.end()); // 升序排列 O(V*logV)

    vector<int> old_to_new(N);
    for (int j = 0; j < N; j++) {
        old_to_new[value[j].second] = N - 1 - j; // new ID 0 has the largest value
        new_to_old[N - 1 - j] = value[j].second;
    }

    dgraph_change_IDs_element(graph, N, old_to_new);
}
