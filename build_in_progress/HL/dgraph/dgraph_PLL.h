#pragma once
#include <iostream>
#include <tool_functions/ThreadPool.h>
#include <shared_mutex>
#include <chrono>
#include <unordered_map>
#include <boost/heap/fibonacci_heap.hpp>
#include <dgraph_v_of_v/dgraph_v_of_v.h>
#include <build_in_progress/HL/dgraph/dgraph_two_hop_label.h>

/* struct used for dijkstra extra_min */
struct node_for_dij
{
public:
    int vertex;
    two_hop_weight_type priority_value;
};

bool operator<(node_for_dij const& x, node_for_dij const& y)
{
    return x.priority_value > y.priority_value; // < is the max-heap; > is the min heap
}
typedef typename boost::heap::fibonacci_heap<node_for_dij>::handle_type dgraph_heap_pointer;

vector<vector<dgraph_heap_pointer>> Q_pointers;



void dgraph_pruned_dijkstra(int v_k, int N, int in_out, dgraph_v_of_v<two_hop_weight_type>* input_graph) {

    mtx_595[max_N_ID_for_mtx_595 - 1].lock();
    int thread_id = Qid_595.front();
    Qid_595.pop();
    mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

    vector<int> dij_dist_changed_vertices, L_vk_changed_vertices; 

    /*hash label distances in L_temp_in[v_k] or L_temp_out[v_k]*/
    mtx_595[v_k].lock();
    if (in_out == 0) {
        int L_out_vk_size = L_temp_out[v_k].size();
        for (int i = 0; i < L_out_vk_size; i++) {
            int x = L_temp_out[v_k][i].vertex;
            L_vk_tmp[thread_id][x] = L_temp_out[v_k][i].distance;
            L_vk_changed_vertices.push_back(x);
        }
    }
    else {
        int L_in_vk_size = L_temp_in[v_k].size();
        for (int i = 0; i < L_in_vk_size; i++) {
            int x = L_temp_in[v_k][i].vertex;
            L_vk_tmp[thread_id][x] = L_temp_in[v_k][i].distance;
            L_vk_changed_vertices.push_back(x);
        }
    }
    mtx_595[v_k].unlock();
   
    node_for_dij node;
    boost::heap::fibonacci_heap<node_for_dij> Q;

    two_hop_label xx;
    node.vertex = v_k;
    node.priority_value = 0;
    Q_pointers[thread_id][v_k] = Q.push(node);

    dij_dist[thread_id][v_k] = 0;
    dij_dist_changed_vertices.push_back(v_k);

    long long int new_label_num = 0;

    while (Q.size()) {
        node = Q.top();
        Q.pop();
        int u = node.vertex;
        two_hop_weight_type dist_u = node.priority_value;

        if (v_k <= u) {
            two_hop_weight_type query_vk_u = std::numeric_limits<two_hop_weight_type>::max();
            mtx_595[u].lock();
            if (in_out == 0) {
                int L_in_u_size = L_temp_in[u].size();
                for (int i = 0; i < L_in_u_size; i++) {
                    two_hop_weight_type dis = L_temp_in[u][i].distance + L_vk_tmp[thread_id][L_temp_in[u][i].vertex];
                    if (query_vk_u > dis) {
                        query_vk_u = dis;
                    }
                }
            }
            else { // reverse case
                int L_out_u_size = L_temp_out[u].size();
                for (int i = 0; i < L_out_u_size; i++) {
                    two_hop_weight_type dis = L_temp_out[u][i].distance + L_vk_tmp[thread_id][L_temp_out[u][i].vertex];
                    if (query_vk_u > dis) {
                        query_vk_u = dis;
                    }
                }
            }
            mtx_595[u].unlock();
            if (query_vk_u <= dist_u)
                continue;

            
            xx.vertex = v_k;
            xx.distance = dist_u;
            new_label_num++;
            if (in_out == 0) {             
                mtx_595[u].lock();
                L_temp_in[u].push_back(xx); /* add new label - in part */
                mtx_595[u].unlock();

                int u_adj_size = input_graph->OUTs[u].size();
                for (int i = 0; i < u_adj_size; i++) {
                    int adj_v = input_graph->OUTs[u][i].first;
                    two_hop_weight_type new_dist = input_graph->OUTs[u][i].second + dist_u;
                    if (dij_dist[thread_id][adj_v] == std::numeric_limits<two_hop_weight_type>::max()) {
                        node.vertex = adj_v;
                        node.priority_value = new_dist;
                        Q_pointers[thread_id][adj_v] = Q.push(node);
                        dij_dist[thread_id][adj_v] = new_dist;
                        dij_dist_changed_vertices.push_back(adj_v);
                    }
                    else { // v is already in the dij_dist, only need to check if update the dij_dist
                        if (dij_dist[thread_id][adj_v] > new_dist) {
                            node.vertex = adj_v;
                            node.priority_value = new_dist;
                            Q.update(Q_pointers[thread_id][adj_v], node);
                            dij_dist[thread_id][adj_v] = new_dist;
                        }
                    }
                }
            }
            else {
                mtx_595[u].lock();
                L_temp_out[u].push_back(xx); /* add new label - out part */
                mtx_595[u].unlock();

                int u_adj_size = input_graph->INs[u].size();
                for (int i = 0; i < u_adj_size; i++) {
                    int adj_v = input_graph->INs[u][i].first;
                    two_hop_weight_type new_dist = input_graph->INs[u][i].second + dist_u;
                    if (dij_dist[thread_id][adj_v] == std::numeric_limits<two_hop_weight_type>::max()) {
                        node.vertex = adj_v;
                        node.priority_value = new_dist;
                        Q_pointers[thread_id][adj_v] = Q.push(node);
                        dij_dist[thread_id][adj_v] = new_dist;
                        dij_dist_changed_vertices.push_back(adj_v);
                    }
                    else { // v is already in the dij_dist, only need to check if update the dij_dist
                        if (dij_dist[thread_id][adj_v] > new_dist) {
                            node.vertex = adj_v;
                            node.priority_value = new_dist;
                            Q.update(Q_pointers[thread_id][adj_v], node);
                            dij_dist[thread_id][adj_v] = new_dist;
                        }
                    }
                }
            }
        }
    }

    /*recover global hash vectors*/
    for (int i : dij_dist_changed_vertices) {
        dij_dist[thread_id][i] = std::numeric_limits<two_hop_weight_type>::max();
    }
    for (int i : L_vk_changed_vertices) {
        L_vk_tmp[thread_id][i] = std::numeric_limits<two_hop_weight_type>::max();
    }

    /* recycle Qid_595 and update labal_size_595*/
    mtx_595[max_N_ID_for_mtx_595 - 1].lock();
    Qid_595.push(thread_id);
    labal_size_595 = labal_size_595 + new_label_num;
    mtx_595[max_N_ID_for_mtx_595 - 1].unlock();
}

void dgraph_PLL(dgraph_v_of_v<two_hop_weight_type>& input_graph, int num_of_threads, dgraph_case_info_v1 &case_info) {

    int N = input_graph.INs.size();
    L_temp_in.resize(N);
    L_temp_out.resize(N);

    ThreadPool pool(num_of_threads);
    std::vector<std::future<int>> results;

    /*dgraph_pruned_dijkstra*/
    queue<int>().swap(Qid_595);
    dij_dist.resize(num_of_threads);
    L_vk_tmp.resize(num_of_threads);
    Q_pointers.resize(num_of_threads);
    for (int i = 0; i < num_of_threads; i++) {
        dij_dist[i].resize(N, std::numeric_limits<two_hop_weight_type>::max());
        L_vk_tmp[i].resize(N, std::numeric_limits<two_hop_weight_type>::max());
        Q_pointers[i].resize(N);
        Qid_595.push(i);
    }
    auto it_g = &input_graph;
    for (int v_k = 0; v_k < N; v_k++) {   
        results.emplace_back(pool.enqueue([v_k, N, it_g] {
                dgraph_pruned_dijkstra(v_k, N, 0, it_g);
                dgraph_pruned_dijkstra(v_k, N, 1, it_g);
                return 1; }));
    }
    for (auto &&result : results)
        result.get();
    results.clear();

    /*sort labels*/
    for (int i = 0; i < N; i++) {
        results.emplace_back(pool.enqueue([i] {
            sort(L_temp_in[i].begin(), L_temp_in[i].end(),
                compare_two_hop_label_vertex_small_to_large); // sort is necessary
            sort(L_temp_out[i].begin(), L_temp_out[i].end(),
                compare_two_hop_label_vertex_small_to_large); // sort is necessary
            return 1; }));
    }
    for (auto&& result : results)
        result.get();
    results.clear();

    /*canonical_repair*/
    if (case_info.use_canonical_repair) {
        case_info.label_size_before_canonical_repair = compute_label_bit_size(L_temp_in, L_temp_out);
        canonical_repair_multi_threads(num_of_threads, &case_info.L_in, &case_info.L_out);
        case_info.label_size_after_canonical_repair = compute_label_bit_size(case_info.L_in, case_info.L_out);
    }
    else {
        case_info.L_in = L_temp_in;
        case_info.L_out = L_temp_out;
    }

    dgraph_clear_global_values_PLL();
    vector<vector<dgraph_heap_pointer>>().swap(Q_pointers);
}
