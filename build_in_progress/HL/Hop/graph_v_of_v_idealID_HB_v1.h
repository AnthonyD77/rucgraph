#pragma once
#include <boost/heap/fibonacci_heap.hpp>
#include <build_in_progress/HL/Hop/graph_v_of_v_idealID_HB_canonical_repair.h>

using namespace std;

struct HBPLL_v1_node {
public:
    int vertex, parent_vertex, hop;
    double priority_value;
};

bool operator<(HBPLL_v1_node const &x, HBPLL_v1_node const &y) {
    return x.priority_value > y.priority_value;
}
typedef typename boost::heap::fibonacci_heap<HBPLL_v1_node>::handle_type graph_v_of_v_idealID_HL_PLL_v1_handle_t_for_sp;

void graph_v_of_v_idealID_HL_HB_v1_thread_function_HBDIJ(int v_k, int N, int upper_k) {
    bool print = false;

    mtx_599[max_N_ID_for_mtx_599 - 1].lock();
    int used_id = Qid_599.front();
    Qid_599.pop();
    mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

    queue<int> P_changed_vertices, T_changed_vertices;

    HBPLL_v1_node node;
    node.vertex = v_k;
    node.parent_vertex = v_k;
    node.hop = 0;
    node.priority_value = 0;

    /* use unordered_map to avoid the influence of upper_k */
    boost::heap::fibonacci_heap<HBPLL_v1_node> Q;
    two_hop_label_v1 xx;
    Q.push(node);

    /* 
        P_dij_599 stores the shortest distance from vk to any other vertices with its hop_cst,
        note that the hop_cst is determined by the shortest distance
     */
    P_dij_599[used_id][v_k] = {0, 0};
    P_changed_vertices.push(v_k);

    /* T_dij_599 stores the label (dist and hop) of vertex v_k */
    mtx_599[v_k].lock();
    int L_v_k_size = L_temp_599[v_k].size();
    for (int i = 0; i < L_v_k_size; i++) {
        int L_v_k_i_vertex = L_temp_599[v_k][i].vertex;
        // TODO: 这个地方有问题, 没有考虑同 vertex 不同 hop 的情况
        T_dij_599[used_id][L_v_k_i_vertex].first = L_temp_599[v_k][i].distance;
        T_dij_599[used_id][L_v_k_i_vertex].second = L_temp_599[v_k][i].hop;
        T_changed_vertices.push(L_v_k_i_vertex);
    }
    mtx_599[v_k].unlock();

    long long int new_label_num = 0;

    if (v_k == -1)
        print = 1;
    if (print)
        cout << "-----begin from " << v_k << "-----" << endl;

    while (Q.size() > 0) {
        node = Q.top();
        Q.pop();
        int u = node.vertex;

        if (print)
            cout << "\t visit " << u << " with hop " << node.hop << "  and dist  " << node.priority_value << endl;

        if (v_k <= u) {  // rank pruning, r(v_k) > r(u)
            int u_parent = node.parent_vertex;
            int u_hop = node.hop;
            double P_u = node.priority_value;
            double P_u_with_error = P_u + 1e-5;
            double query_v_k_u = std::numeric_limits<double>::max();

            /* there are two upper_k judge in total */
            if (u_hop > upper_k)
                break;

#ifdef _WIN32
            mtx_599[u].lock();
            auto L_u_size = L_temp_599[u].size();  // a vector<PLL_with_non_adj_reduction_sorted_label>
            mtx_599[u].unlock();
            for (int i = 0; i < L_u_size; i++) {
                mtx_599[u].lock();  // put lock in for loop is very slow, but it may be the only way under Windows
                double dis = L_temp_599[u][i].distance + T_bfs_599[used_id][L_temp_599[u][i].vertex];
                mtx_599[u].unlock();
                if (query_v_k_u > dis) {
                    query_v_k_u = dis;
                }
            }  //求query的值
#else
            mtx_599[u].lock();
            auto L_u_size1 = L_temp_599[u].size();
            for (int i = 0; i < L_u_size1; i++) {
                if (L_temp_599[u][i].hop + T_dij_599[used_id][L_temp_599[u][i].vertex].second <= u_hop) {
                    double dis = L_temp_599[u][i].distance + T_dij_599[used_id][L_temp_599[u][i].vertex].first;
                    if (query_v_k_u > dis) {
                        query_v_k_u = dis;
                    }
                }
            }
            if (print)
                cout << "\t query " << v_k << " and " << u << " = " << query_v_k_u << endl;
            mtx_599[u].unlock();
#endif

            if (print) {
                cout << "\t P_u_with_error: " << P_u_with_error << " || query_v_k_u: " << query_v_k_u << endl;
            }

            if (P_u_with_error < query_v_k_u) {  //pruning
                xx.vertex = v_k;
                xx.distance = P_u;
                xx.parent_vertex = u_parent;
                xx.hop = u_hop;

                mtx_599[u].lock();
                L_temp_599[u].push_back(xx);
                if (print)
                    cout << "\t\t add " << v_k << " into " << u << endl;
                mtx_599[u].unlock();
                new_label_num++;

                /* update adj */
                if (print)
                    cout << "\t update adj" << endl;
                int u_adj_size = ideal_graph_599[u].size();
                for (int i = 0; i < u_adj_size; i++) {
                    int adj_v = ideal_graph_599[u][i].first;
                    double ec = ideal_graph_599[u][i].second;

                    /* update node info */
                    node.vertex = adj_v;
                    node.parent_vertex = u;
                    node.priority_value = P_u + ec;
                    node.hop = u_hop + 1;

                    if (print)
                        cout << "\t neighber: " << adj_v << ",, dist: " << node.priority_value << ",, hop: " << node.hop << endl;
                    // << "P_dij_599[used_id][adj_v].first: " << P_dij_599[used_id][adj_v].first << endl
                    // << "P_dij_599[used_id][adj_v].second: " << P_dij_599[used_id][adj_v].second << endl;

                    /* check if the edge is generated by R2: hop+1 */
                    double new_dist = node.priority_value;
                    if (new_edges_with_origin_ec_599.find({u, adj_v}) != new_edges_with_origin_ec_599.end() && adj_v != u_parent && adj_v != v_k) /* imp check to prevent infinite loop  */
                    {
                        if (new_edges_with_origin_ec_599[{u, adj_v}] != std::numeric_limits<double>::max()) {
                            if (ec < new_edges_with_origin_ec_599[{u, adj_v}]) {
                                /* add 1-edge (origin) info, here do not need to modify parent */
                                node.priority_value = new_edges_with_origin_ec_599[{u, adj_v}] + P_u;
                                node.parent_vertex = -u;
                                if (u == 0)
                                    node.parent_vertex = std::numeric_limits<int>::max();
                                Q.push(node);

                                // if (u==0 && adj_v==6)
                                //     cout << "\t push " << node.vertex << "," << node.priority_value << "," << node.parent_vertex << "," << node.hop << " into Q" << endl;

                                node.parent_vertex = u;
                            }
                        }
                        /* update node info of 2-edge */
                        node.hop++;
                        node.priority_value = new_dist;
                    }
                    /* beyond upper_k, then stop expansion */
                    if (node.hop > upper_k) {
                        if (print)
                            cout << "\t beyond upper_k" << endl;
                        break;
                    }

                    /* 
                        vertices not reached yet:
                        just add the distance and hop info
                    */
                    if (P_dij_599[used_id][adj_v].first == std::numeric_limits<double>::max()) {
                        Q.push(node);
                        P_dij_599[used_id][adj_v].first = node.priority_value;
                        P_dij_599[used_id][adj_v].second = node.hop;
                        P_changed_vertices.push(adj_v);
                        if (print)
                            cout << "\t 1 add new " << adj_v << " with " << node.priority_value << endl;
                    }
                    /*
                        vertices already reached:
                        1. smaller distance, then update info
                        2. greater distance but smaller hop, add new info, do not update P_dij_599
                    */
                    else {
                        if (node.priority_value < P_dij_599[used_id][adj_v].first) {
                            Q.push(node);
                            // Q.update(Q_handles[node.hop][adj_v], node);
                            P_dij_599[used_id][adj_v].first = node.priority_value;
                            P_dij_599[used_id][adj_v].second = node.hop;
                            if (print)
                                cout << "\t 2 add neighbor " << adj_v << " with " << node.priority_value << endl;
                        } else if (node.hop < P_dij_599[used_id][adj_v].second) {
                            Q.push(node);
                            if (print)
                                cout << "\t 3 add neighbor " << adj_v << " with " << node.priority_value << endl;
                        }
                    }
                }
                /* stop update adj */
            }
        }
    }

    while (P_changed_vertices.size() > 0) {
        // P_dij_599[used_id][P_changed_vertices.front()].clear();
        P_dij_599[used_id][P_changed_vertices.front()].first = std::numeric_limits<double>::max();
        P_changed_vertices.pop();
    }
    while (T_changed_vertices.size() > 0) {
        T_dij_599[used_id][T_changed_vertices.front()].first = std::numeric_limits<double>::max();  // reverse-allocate T values
        T_changed_vertices.pop();
    }

    mtx_599[v_k].lock();
    vector<two_hop_label_v1>(L_temp_599[v_k]).swap(L_temp_599[v_k]);  // swap释放vector中多余空间： https://blog.csdn.net/qq_41929943/article/details/103190891
    mtx_599[v_k].unlock();

    mtx_599[max_N_ID_for_mtx_599 - 1].lock();
    Qid_599.push(used_id);
    labal_size_599 = labal_size_599 + new_label_num;
    mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

    if (labal_size_599 > max_labal_size_599) {
        throw reach_limit_error_string_MB;  // after catching error, must call terminate_procedures_599(), otherwise this PLL cannot be reused
    }

    if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time_599).count() > max_run_time_nanoseconds_599) {
        throw reach_limit_error_string_time;  // after catching error, must call terminate_procedures_599(), otherwise this PLL cannot be reused
    }
}

void graph_v_of_v_idealID_HB_v1_sort_labels_thread(vector<vector<two_hop_label_v1>> *output_L, int v_k, double value_M) {
    sort(L_temp_599[v_k].begin(), L_temp_599[v_k].end(), compare_two_hop_label_small_to_large);
    if (value_M != 0) {
        int size_vk = L_temp_599[v_k].size();
        for (int i = 0; i < size_vk; i++) {
            L_temp_599[v_k][i].distance += L_temp_599[v_k][i].hop * value_M;
        }
    }
    (*output_L)[v_k] = L_temp_599[v_k];
    vector<two_hop_label_v1>().swap(L_temp_599[v_k]);  // clear new labels for RAM efficiency
}

vector<vector<two_hop_label_v1>> graph_v_of_v_idealID_HB_v1_sort_labels(int N, int max_N_ID, int num_of_threads, double value_M = 0) {
    vector<vector<two_hop_label_v1>> output_L(max_N_ID);
    vector<vector<two_hop_label_v1>> *p = &output_L;

    ThreadPool pool(num_of_threads);
    std::vector<std::future<int>> results;
    for (int v_k = 0; v_k < N; v_k++) {
        results.emplace_back(pool.enqueue([p, v_k, value_M] {
            graph_v_of_v_idealID_HB_v1_sort_labels_thread(p, v_k, value_M);
            return 1;
        }));
    }
    for (auto &&result : results)
        result.get();

    return output_L;
}

void graph_v_of_v_idealID_HB_v1(graph_v_of_v_idealID &input_graph, int max_N_ID, bool weighted, int num_of_threads, graph_v_of_v_idealID_two_hop_case_info_v1 &case_info) {
    //----------------------------------- step 1: initialization -----------------------------------
    cout << "step 1: initialization" << endl;

    auto begin = std::chrono::high_resolution_clock::now();
    /*information prepare*/
    labal_size_599 = 0;
    begin_time_599 = std::chrono::high_resolution_clock::now();
    max_run_time_nanoseconds_599 = case_info.max_run_time_seconds * 1e9;
    max_labal_size_599 = case_info.max_labal_size;

    if (max_N_ID > max_N_ID_for_mtx_599) {
        cout << "max_N_ID > max_N_ID_for_mtx_599; max_N_ID_for_mtx_599 is too small!" << endl;
        exit(1);
    }

    L_temp_599.resize(max_N_ID);
    int N = input_graph.size();

    /* thread info */
    ThreadPool pool(num_of_threads);
    std::vector<std::future<int>> results;
    int num_of_threads_per_push = num_of_threads * 100;

    ideal_graph_599 = input_graph;

    auto end = std::chrono::high_resolution_clock::now();
    case_info.time_initialization = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;

    //----------------------------------------------- step 2: reduction ---------------------------------------------------------------
    cout << "step 2: reduction" << endl;

    /*redcution: add and remove certain edges*/
    case_info.reduction_measures_2019R2.clear();  // for using this function multiple times
    case_info.reduction_measures_2019R2.resize(max_N_ID, 0);

    begin = std::chrono::high_resolution_clock::now();
    if (case_info.use_2019R2) {
        case_info.MG_num = 0;
        for (int x = 0; x < N; x++) {
            if (ideal_graph_599[x].size() > 0) {
                if (x > ideal_graph_599[x][ideal_graph_599[x].size() - 1].first) {
                    case_info.reduction_measures_2019R2[x] = 2;
                    R2_reduced_vertices[x] = *new vector<pair<int, double>>(ideal_graph_599[x]);
                }
            }
        }
    } else if (case_info.use_enhanced2019R2) {
        case_info.MG_num = 0;
        for (int x = 0; x < N; x++) {
            if (ideal_graph_599[x].size() > 0) {
                if (x > ideal_graph_599[x][ideal_graph_599[x].size() - 1].first) {
                    case_info.reduction_measures_2019R2[x] = 2;
                    R2_reduced_vertices[x] = *new vector<pair<int, double>>(ideal_graph_599[x]);
                }
            }
        }
        int bound = case_info.max_degree_MG_enhanced2019R2;
        for (int x = N - 1; x >= 0; x--) {
            if (case_info.reduction_measures_2019R2[x] == 0 && ideal_graph_599[x].size() <= bound) {
                bool no_adj_MG_vertices = true;
                for (auto it = ideal_graph_599[x].begin(); it != ideal_graph_599[x].end(); it++) {
                    if (case_info.reduction_measures_2019R2[it->first] == 2) {
                        no_adj_MG_vertices = false;
                        break;
                    }
                }
                if (no_adj_MG_vertices) {
                    case_info.reduction_measures_2019R2[x] = 2;
                    R2_reduced_vertices[x] = *new vector<pair<int, double>>(ideal_graph_599[x]);
                }
            }
        }
    } else if (case_info.use_non_adj_reduc_degree) {
        case_info.MG_num = 0;
        int bound = case_info.max_degree_MG_enhanced2019R2;
        for (int x = N - 1; x >= 0; x--) {
            if (case_info.reduction_measures_2019R2[x] == 0 && ideal_graph_599[x].size() <= bound) {
                bool no_adj_MG_vertices = true;
                for (auto it = ideal_graph_599[x].begin(); it != ideal_graph_599[x].end(); it++) {
                    if (case_info.reduction_measures_2019R2[it->first] == 2) {
                        no_adj_MG_vertices = false;
                        break;
                    }
                }
                if (no_adj_MG_vertices) {
                    case_info.reduction_measures_2019R2[x] = 2;
                    R2_reduced_vertices[x] = *new vector<pair<int, double>>(ideal_graph_599[x]);
                }
            }
        }
    }
    for (int x = N - 1; x >= 0; x--) {
        if (case_info.reduction_measures_2019R2[x] == 2) {
            /*add edge*/
            auto it1 = ideal_graph_599[x].begin();
            for (int m = ideal_graph_599[x].size() - 1; m > 0; m--) {
                for (int n = m - 1; n >= 0; n--) {
                    double s_vPLUSv_t = (it1 + m)->second + (it1 + n)->second;
                    int e1 = (it1 + m)->first;
                    int e2 = (it1 + n)->first;
                    double origin_ec = graph_v_of_v_idealID_edge_weight(ideal_graph_599, e1, e2);
                    if (s_vPLUSv_t < origin_ec) {
                        graph_v_of_v_idealID_add_edge(ideal_graph_599, e1, e2, s_vPLUSv_t);
                        new_edges_with_middle_v_599[{e1, e2}] = x;
                        new_edges_with_middle_v_599[{e2, e1}] = x;
                        /* only record the origin edge, in case of the edge is updated many times */
                        if (new_edges_with_origin_ec_599.find({e1, e2}) == new_edges_with_origin_ec_599.end()) {
                            new_edges_with_origin_ec_599[{e1, e2}] = origin_ec;
                            new_edges_with_origin_ec_599[{e2, e1}] = origin_ec;
                        }
                    }
                }
            }
            /*remove edge*/
            case_info.MG_num++;
            graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_599, x);
        }
    }
    end = std::chrono::high_resolution_clock::now();
    case_info.time_reduction = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;  // s

    //----------------------------------------------- step 3: generate labels ---------------------------------------------------------------
    cout << "step 3: generate labels" << endl;
    begin = std::chrono::high_resolution_clock::now();

    /*seaching shortest paths*/
    int upper_k = case_info.upper_k == 0 ? std::numeric_limits<int>::max() : case_info.upper_k;

    P_dij_599.resize(num_of_threads);
    T_dij_599.resize(num_of_threads);
    for (int i = 0; i < num_of_threads; i++) {
        P_dij_599[i].resize(N);
        T_dij_599[i].resize(N);
        for (int j = 0; j < N; j++) {
            P_dij_599[i][j] = {std::numeric_limits<double>::max(), 0};
            T_dij_599[i][j] = {std::numeric_limits<double>::max(), 0};
        }
        Qid_599.push(i);
    }
    int push_num = 0;
    for (int v_k = 0; v_k < N; v_k++) {
        if (ideal_graph_599[v_k].size() > 0) {
            results.emplace_back(
                pool.enqueue([v_k, N, upper_k] {
                    graph_v_of_v_idealID_HL_HB_v1_thread_function_HBDIJ(v_k, N, upper_k);
                    return 1;
                }));
            push_num++;
        }
        if (push_num % num_of_threads_per_push == 0) {
            for (auto &&result : results)
                result.get();
            results.clear();
        }
    }

    for (auto &&result : results)
        result.get();

    end = std::chrono::high_resolution_clock::now();
    case_info.time_generate_labels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;  // s

    /*
        update predecessors for this non_adj_reduction,
        this update is for correct recursive direction.
    */
    begin = std::chrono::high_resolution_clock::now();
    for (auto it = new_edges_with_middle_v_599.begin(); it != new_edges_with_middle_v_599.end(); it++) {
        int e1 = it->first.first;
        int e2 = it->first.second;
        int middle_k = it->second;
        for (int j = L_temp_599[e1].size() - 1; j >= 0; j--) {
            if (L_temp_599[e1][j].parent_vertex == e2) {
                L_temp_599[e1][j].parent_vertex = middle_k;
            }
        }
    }
    for (auto it = new_edges_with_middle_v_599.begin(); it != new_edges_with_middle_v_599.end(); it++) {
        int e1 = it->first.first;
        int e2 = it->first.second;
        int middle_k = it->second;
        for (int j = L_temp_599[e1].size() - 1; j >= 0; j--) {
            if (L_temp_599[e1][j].parent_vertex < 0) {
                L_temp_599[e1][j].parent_vertex *= (-1);
                continue;
            }
            if (L_temp_599[e1][j].parent_vertex == std::numeric_limits<int>::max()) {
                L_temp_599[e1][j].parent_vertex = 0;
                continue;
            }
        }
    }
    end = std::chrono::high_resolution_clock::now();
    case_info.time_update_predecessors = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;  // s

    //----------------------------------------------- step 4: canonical_repair---------------------------------------------------------------
    cout << "step 4: canonical_repair" << endl;

    if (case_info.use_canonical_repair) {
        begin = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < max_N_ID; i++) {
            sort(L_temp_599[i].begin(), L_temp_599[i].end(), compare_two_hop_label_small_to_large);  // sort is necessary
        }

        /* print label before canonical repair if necessary */
        if (case_info.print_label_before_canonical_fix) {
            cout << "label_before_canonical_fix:" << endl;
            for (int i = 0; i < L_temp_599.size(); i++) {
                cout << "L[" << i << "]:\t";
                int ii = i;
                for (int j = 0; j < L_temp_599[ii].size(); j++) {
                    cout << "{" << L_temp_599[ii][j].vertex << "," << L_temp_599[ii][j].distance << "," << L_temp_599[ii][j].parent_vertex << "," << L_temp_599[ii][j].hop << "}";
                }
                cout << endl;
            }
        }

        canonical_repair_multi_threads(case_info.label_size_before_canonical_repair,
                                       case_info.label_size_after_canonical_repair,
                                       case_info.canonical_repair_remove_label_ratio, num_of_threads);
        end = std::chrono::high_resolution_clock::now();
        case_info.time_canonical_repair = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
    }

    //----------------------------------------------- step 5: update_labels ---------------------------------------------------------------
    cout << "step 5: update_old_IDs_in_labels" << endl;
    begin = std::chrono::high_resolution_clock::now();

    /* sort lables */
    if (case_info.use_M)
        case_info.L = graph_v_of_v_idealID_HB_v1_sort_labels(N, max_N_ID, num_of_threads, case_info.value_M);
    else
        case_info.L = graph_v_of_v_idealID_HB_v1_sort_labels(N, max_N_ID, num_of_threads);

    end = std::chrono::high_resolution_clock::now();
    case_info.time_update_labels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;  // s

    graph_v_of_v_idealID_two_hop_clear_global_values();
}