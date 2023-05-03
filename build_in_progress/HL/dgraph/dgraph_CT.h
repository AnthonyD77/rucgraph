#pragma once
#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <unordered_set>
#include <dgraph_v_of_v/dgraph_v_of_v.h>
#include <build_in_progress/HL/dgraph/dgraph_two_hop_label.h>
#include <build_in_progress/HL/dgraph/dgraph_PLL.h>
#include <build_in_progress/HL/dgraph/dgraph_PSL.h>

struct node_degree
{
    int vertex, degree;
    bool operator<(const node_degree &nd) const
    {
        return degree > nd.degree;
    }
};

class dgraph_case_info_v2 {
  public:
    /* parameters */
    int thread_num = 1;
    bool use_PLL = 1;
    int d = 0;
    int two_hop_order_method = 0;

    /* labels */
    dgraph_case_info_v1 two_hop_case_info;
    std::vector<std::vector<std::pair<int, two_hop_weight_type>>> Bags_in, Bags_out;
    std::vector<bool> isIntree;
    std::vector<int> root;
    std::vector<std::vector<int>> tree_st;   // for lca
    std::vector<std::vector<int>> tree_st_r; // for lca
    std::vector<int> first_pos;              // for lca
    std::vector<int> lg;                     // for lca
    std::vector<int> dep;                    // for lca

    /*clear labels*/
    void clear_labels() {
        two_hop_case_info.clear_labels();
        std::vector<std::vector<pair<int, two_hop_weight_type>>>().swap(Bags_in);
        std::vector<std::vector<pair<int, two_hop_weight_type>>>().swap(Bags_out);
        vector<bool>().swap(isIntree);
        vector<int>().swap(root);
        vector<vector<int>>().swap(tree_st);
        vector<vector<int>>().swap(tree_st_r);
        vector<int>().swap(first_pos);
        vector<int>().swap(lg);
        vector<int>().swap(dep);
    }

    /*running limits*/
    long long int max_bit_size = 1e12;
    double max_run_time_seconds = 1e12; // s

    /*indexing times*/
    double time1_initialization = 0;
    double time2_tree_decomposition = 0;
    double time3_tree_indexs = 0;
    double time4_lca = 0;
    double time5_core_indexs_prepare1 = 0;
    double time5_core_indexs_prepare2 = 0;
    double time5_core_indexs_prepare3 = 0;
    double time5_core_indexs_prepare4 = 0;
    double time5_core_indexs = 0;
    double time5_core_indexs_post = 0;
    double time6_post = 0;
    double time_total = 0;

    /*core graph info*/
    long long int core_graph_V = 0;
    long long int core_graph_E = 0;
    long long int uk = 0;
    double V_log_V_uk = 0;
    long long int core_graph_label_bit_size = 0;  // bit size 1
    long long int pre_core_graph_label_bit_size = 0; // bit size 0

    long long int total_label_bit_size() {
        return pre_core_graph_label_bit_size + core_graph_label_bit_size;
    }

    /*printing*/
    void print_Bags()
    {
        cout << "print_Bags:" << endl;
        for (int i = 0; i < Bags_in.size(); i++)
        {
            cout << i << ": ";
            for (int j = 0; j < Bags_in[i].size(); j++)
            {
                cout << Bags_in[i][j].first << ", ";
            }
            cout << "||";
            for (int j = 0; j < Bags_out[i].size(); j++)
            {
                cout << Bags_out[i][j].first << ", ";
            }
            cout << endl;
        }
    }

    void print_root()
    {
        cout << "print_root:" << endl;
        for (int i = 0; i < root.size(); i++)
        {
            cout << "root[" << i << "]=" << root[i] << endl;
        }
    }

    void print_isIntree()
    {
        cout << "print_isIntree:" << endl;
        for (int i = 0; i < isIntree.size(); i++)
        {
            cout << "isIntree[" << i << "]=" << isIntree[i] << endl;
        }
    }

    void print_time() {
        cout << "time1_initialization: " << time1_initialization << endl;
        cout << "time2_tree_decomposition: " << time2_tree_decomposition << endl;
        cout << "time3_tree_indexs: " << time3_tree_indexs << endl;
        cout << "time4_lca: " << time4_lca << endl;
        cout << "time5_core_indexs_prepare1: " << time5_core_indexs_prepare1 << endl;
        cout << "time5_core_indexs_prepare2: " << time5_core_indexs_prepare2 << endl;
        cout << "time5_core_indexs_prepare3: " << time5_core_indexs_prepare3 << endl;
        cout << "time5_core_indexs_prepare4: " << time5_core_indexs_prepare4 << endl;
        cout << "time5_core_indexs: " << time5_core_indexs << endl;
        cout << "time5_core_indexs_post: " << time5_core_indexs_post << endl;
        cout << "time6_post: " << time6_post << endl;
        cout << "time_total: " << time_total << endl;
        cout << "time1_PLL_PSL_initialization: " << two_hop_case_info.time1_PLL_PSL_initialization << endl;
        cout << "time2_PLL_PSL_label_generation: " << two_hop_case_info.time2_PLL_PSL_label_generation << endl;
        cout << "time3_PLL_PSL_label_postprocess: " << two_hop_case_info.time3_PLL_PSL_label_postprocess << endl;
        cout << "time4_PLL_PSL_label_canonical: " << two_hop_case_info.time4_PLL_PSL_label_canonical << endl;
        cout << "time5_PLL_PSL_total: " << two_hop_case_info.time5_PLL_PSL_total << endl;
    }

    void print_label_bit_size() {
        cout << "pre_core_graph_label_bit_size: " << pre_core_graph_label_bit_size << endl;
        cout << "core_graph_label_bit_size: " << core_graph_label_bit_size << endl;
        cout << "total_label_bit_size: " << pre_core_graph_label_bit_size + core_graph_label_bit_size << endl;
    }
};

/*global values*/
dgraph_v_of_v<two_hop_weight_type> global_dgraph_CT;
vector<pair<int, int>> infty_edge;
dgraph_case_info_v1 two_hop_case_info_sorted;
vector<int> new_to_old;
int global_N;
long long int sum_uk = 0;
bool this_parallel_CT_is_running = false;

void clear_gloval_values_CT() {
    this_parallel_CT_is_running = false;
    sum_uk = 0;
    global_dgraph_CT.clear();
    vector<pair<int, int>>().swap(infty_edge);
    vector<int>().swap(new_to_old);
    two_hop_case_info_sorted.clear_labels();
    dgraph_clear_global_values_PLL_PSL();
}




long long int compute_CT_label_bit_size(dgraph_case_info_v2& info, ThreadPool& pool, std::vector<std::future<int>>& results) {

    long long int ssize = 0;
    int N = info.Bags_in.size();
    auto& ii = info;

    results.emplace_back(pool.enqueue([&ssize, &ii, N] {
        long long int x = compute_label_bit_size(ii.two_hop_case_info.L_in, ii.two_hop_case_info.L_out);
        x += N * (sizeof(bool) + 3 * sizeof(int)) + ii.lg.size() * sizeof(int); // isIntree, root, first_pos, dep, lg
        mtx.lock();
        ssize += x;
        mtx.unlock();
        return 1; }));
    results.emplace_back(pool.enqueue([&ssize, &ii, N] {
        long long int x = 0;
        for (int i = N - 1; i >= 0; i--) {
            x += ii.Bags_in[i].size() * sizeof(std::pair<int, two_hop_weight_type>);
        }
        mtx.lock();
        ssize += x;
        mtx.unlock();
        return 1; }));
    results.emplace_back(pool.enqueue([&ssize, &ii, N] {
        long long int x = 0;
        for (int i = N - 1; i >= 0; i--) {
            x += ii.Bags_out[i].size() * sizeof(std::pair<int, two_hop_weight_type>);
        }
        mtx.lock();
        ssize += x;
        mtx.unlock();
        return 1; }));
    results.emplace_back(pool.enqueue([&ssize, &ii, N] {
        long long int x = 0;
        for (int i = N - 1; i >= 0; i--) {
            x += ii.Bags_out[i].size() * sizeof(std::pair<int, two_hop_weight_type>);
        }
        mtx.lock();
        ssize += x;
        mtx.unlock();
        return 1; }));
    results.emplace_back(pool.enqueue([&ssize, &ii] {
        long long int x = 0;
        for (int i = ii.tree_st.size() - 1; i >= 0; i--) {
            x += ii.tree_st[i].size() * sizeof(int);
        }
        mtx.lock();
        ssize += x;
        mtx.unlock();
        return 1; }));
    results.emplace_back(pool.enqueue([&ssize, &ii] {
        long long int x = 0;
        for (int i = ii.tree_st_r.size() - 1; i >= 0; i--) {
            x += ii.tree_st_r[i].size() * sizeof(int);
        }
        mtx.lock();
        ssize += x;
        mtx.unlock();
        return 1; }));
    for (auto&& result : results)
        result.get();
    results.clear();

    return ssize;
}

void Label_new_to_old_parallel(int newID, vector<vector<two_hop_label>> &L_in, vector<vector<two_hop_label>>& L_out) {

    int oldID = new_to_old[newID];

    int size = two_hop_case_info_sorted.L_in[newID].size();
    for (int i = 0; i < size; i++) {
        two_hop_case_info_sorted.L_in[newID][i].vertex = new_to_old[two_hop_case_info_sorted.L_in[newID][i].vertex];
    }
    sort(two_hop_case_info_sorted.L_in[newID].begin(), two_hop_case_info_sorted.L_in[newID].end(), compare_two_hop_label_vertex_small_to_large);
    L_in[oldID].swap(two_hop_case_info_sorted.L_in[newID]); // two_hop_case_info_sorted.L_in[newID] becomes empty

    size = two_hop_case_info_sorted.L_out[newID].size();
    for (int i = 0; i < size; i++) {
        two_hop_case_info_sorted.L_out[newID][i].vertex = new_to_old[two_hop_case_info_sorted.L_out[newID][i].vertex];
    }
    sort(two_hop_case_info_sorted.L_out[newID].begin(), two_hop_case_info_sorted.L_out[newID].end(), compare_two_hop_label_vertex_small_to_large);
    L_out[oldID].swap(two_hop_case_info_sorted.L_out[newID]);
}

void substitute_parallel(int u, int w, two_hop_weight_type ec) {
    mtx_595[u].lock();                  // u_out  
    auto old_ec = global_dgraph_CT.edge_weight(u, w);
    if (old_ec == std::numeric_limits<two_hop_weight_type>::max() || old_ec > ec) { // ec may be std::numeric_limits<two_hop_weight_type>::max()
        mtx_595[w + global_N].lock();       // w_in
        global_dgraph_CT.add_edge(u, w, ec);
        mtx_595[w + global_N].unlock();
    }
    mtx_595[u].unlock();   
}

void substitute_parallel_2(int u, int w) {
    mtx_595[w].lock();                  // w_out
    bool contain_wu = global_dgraph_CT.contain_edge(w, u);
    mtx_595[w].unlock();

    mtx_595[u].lock();                  // u_out  
    if (global_dgraph_CT.contain_edge(u, w) == 0 && contain_wu == 0) {
        mtx_595[w + global_N].lock();       // w_in
        edge_lock.lock();
        global_dgraph_CT.add_edge(u, w, std::numeric_limits<two_hop_weight_type>::max());
        infty_edge.push_back({u, w});
        edge_lock.unlock();
        mtx_595[w + global_N].unlock();
    }  
    mtx_595[u].unlock();   
}

void graph_hash_of_mixed_weighted_binary_operations_erase_infty(std::vector<std::pair<int, two_hop_weight_type>>& input_vector, int key) {

    /*erase key from vector; time complexity O(log n + size()-position ), which is O(n) in the worst case, as
    the time complexity of erasing an element from a vector is the number of elements behind this element*/

    if (input_vector.size() > 0) {
        int left = 0, right = input_vector.size() - 1;

        while (left <= right) {
            int mid = left + ((right - left) / 2);
            if (input_vector[mid].first == key && input_vector[mid].second == std::numeric_limits<two_hop_weight_type>::max()) {
                input_vector.erase(input_vector.begin() + mid);
                break;
            }
            else if (input_vector[mid].first > key) {
                right = mid - 1;
            }
            else {
                left = mid + 1;
            }
        }
    }

}

void dgraph_v_of_v_remove_edge_infty(int v1, int v2) {

    /*edge direction: v1 to v2*/
    mtx_595[v1].lock(); // out lock
    graph_hash_of_mixed_weighted_binary_operations_erase_infty(global_dgraph_CT.OUTs[v1], v2); 
    mtx_595[v1].unlock();
    mtx_595[v2 + global_N].lock(); // in lock
    graph_hash_of_mixed_weighted_binary_operations_erase_infty(global_dgraph_CT.INs[v2], v1); 
    mtx_595[v2 + global_N].unlock();
}

void dfs(int &total, vector<int> &first_pos, int x, vector<vector<int>> &son, vector<int> &dfn)
{
    total++;
    dfn[total] = x;
    first_pos[x] = total;
    int s_size = son[x].size();
    for (int i = 0; i < s_size; i++)
    {
        dfs(total, first_pos, son[x][i], son, dfn); // this is a recursive function
        total++;
        dfn[total] = x;
    }
}


/*choose_PLL_PSL*/

void dgraph_dijstar_calculate_uk(dgraph_v_of_v<two_hop_weight_type>* input_graph_pointer, int v_k) {

    int uk = 0;
    auto& INs = input_graph_pointer->INs;
    auto& OUTs = input_graph_pointer->OUTs;
    int N = INs.size();
    vector<two_hop_weight_type> distances(N, std::numeric_limits<two_hop_weight_type>::max());

    node_for_shortest_distances node;
    boost::heap::fibonacci_heap<node_for_shortest_distances> Q;
    vector<heap_pointer_dgraph_shortest_distances_source_to_all> Q_pointer(N);

    /*out search*/
    distances[v_k] = 0;
    node.vertex = v_k;
    node.priority_value = 0;
    Q_pointer[v_k] = Q.push(node);
    while (Q.size() > 0) {
        node = Q.top();
        Q.pop();
        int u = node.vertex;
        auto dist_u = node.priority_value;

        int u_adj_size = OUTs[u].size();
        for (int i = 0;i < u_adj_size;i++) {
            int adj_v = OUTs[u][i].first;
            auto new_dist = OUTs[u][i].second + dist_u;
            if (distances[adj_v] == std::numeric_limits<two_hop_weight_type>::max()) {
                node.vertex = adj_v;
                node.priority_value = new_dist;
                Q_pointer[adj_v] = Q.push(node);
                distances[adj_v] = node.priority_value;
            }
            else {
                if (distances[adj_v] > new_dist) {
                    uk++;
                    node.vertex = adj_v;
                    node.priority_value = new_dist;
                    Q.update(Q_pointer[adj_v], node);
                    distances[adj_v] = node.priority_value;
                }
            }
        }
    }

    /*in search*/
    distances.assign(N, std::numeric_limits<two_hop_weight_type>::max());
    distances[v_k] = 0;
    node.vertex = v_k;
    node.priority_value = 0;
    Q_pointer[v_k] = Q.push(node);
    while (Q.size() > 0) {
        node = Q.top();
        Q.pop();
        int u = node.vertex;
        auto dist_u = node.priority_value;

        int u_adj_size = INs[u].size();
        for (int i = 0;i < u_adj_size;i++) {
            int adj_v = INs[u][i].first;
            auto new_dist = INs[u][i].second + dist_u;
            if (distances[adj_v] == std::numeric_limits<two_hop_weight_type>::max()) {
                node.vertex = adj_v;
                node.priority_value = new_dist;
                Q_pointer[adj_v] = Q.push(node);
                distances[adj_v] = node.priority_value;
            }
            else {
                if (distances[adj_v] > new_dist) {
                    uk++;
                    node.vertex = adj_v;
                    node.priority_value = new_dist;
                    Q.update(Q_pointer[adj_v], node);
                    distances[adj_v] = node.priority_value;
                }
            }
        }
    }


    mtx_595[0].lock();
    sum_uk += uk;
    mtx_595[0].unlock();
}

void choose_PLL_PSL(dgraph_case_info_v2& case_info, dgraph_v_of_v<two_hop_weight_type>& dgraph, ThreadPool& pool, std::vector<std::future<int>>& results) {

    int N = dgraph.INs.size();
    auto& V = case_info.core_graph_V; // number of not isolated vertices
    auto& E = case_info.core_graph_E;
    for (int i = 0; i < N; i++) {
        if (dgraph.INs[i].size()) {
            V++;
            E += dgraph.INs[i].size();
        }
    }
    int log2V = ceil(log2(V)); 
    auto* dgraph_pointer = &dgraph;
    for (int i = 0; i < log2V; i++) {
        results.emplace_back(pool.enqueue([dgraph_pointer, i] {
            dgraph_dijstar_calculate_uk(dgraph_pointer, i);
            return 1; }));
    }
    for (auto&& result : results)
        result.get();
    results.clear();

    case_info.uk = sum_uk;
    case_info.V_log_V_uk = log2(V) * V / (double) sum_uk; // in unweighted or special weighted graphs, sum_uk=0, and V_log_V_uk=inf
}


/*indexing function*/
void CT_dgraph(dgraph_v_of_v<two_hop_weight_type> &input_graph, dgraph_case_info_v2 &case_info) {

    mtx_595[max_N_ID_for_mtx_595 - 1].lock();
    if (this_parallel_CT_is_running == true) {
        cout << "CT_dgraph cannot be run parallelly, due to the above (static) globel values" << endl;
        exit(1);
    }
    this_parallel_CT_is_running = true;
    mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

    //--------------------------- step 1: initialization ---------------------------
    auto begin1 = std::chrono::high_resolution_clock::now();

    auto &Bags_in = case_info.Bags_in;
    auto &Bags_out = case_info.Bags_out;
    auto &isIntree = case_info.isIntree;
    auto &root = case_info.root;
    auto &tree_st = case_info.tree_st;
    auto &tree_st_r = case_info.tree_st_r;
    auto &first_pos = case_info.first_pos;
    auto &lg = case_info.lg;
    auto &dep = case_info.dep;
    int N = input_graph.INs.size();
    global_N = N + 1;
    global_dgraph_CT = input_graph;
    isIntree.resize(N, 0);
    /* initialize queue */
    priority_queue<node_degree> q;
    node_degree nd;
    for (int i = 0; i < N; i++) {
        nd.degree = global_dgraph_CT.degree(i);
        nd.vertex = i;
        q.push(nd);
    }

    case_info.time1_initialization = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin1).count() / 1e9;
    //---------------------------------------------------------------------------------

    //--------------------------- step 2: MDE-based tree decomposition ---------------------------
    auto begin2 = std::chrono::high_resolution_clock::now();

    /* MDE */
    int bound_lambda = N;
    Bags_in.resize(N);
    Bags_out.resize(N);
    vector<int> node_order(N + 1); // merging ID to original ID
    ThreadPool pool(case_info.thread_num);
    std::vector<std::future<int>> results;
    for (int i = 1; i <= N; i++) {
        /* 1. find v_i */
        while (1) {
            nd = q.top();
            q.pop();
            if (!isIntree[nd.vertex] && global_dgraph_CT.degree(nd.vertex) == nd.degree)
                break;
        }
        int v_i = nd.vertex;
        /* 2. Determine whether to break */
        if (nd.degree >= case_info.d) {
            bound_lambda = i - 1;
            q.push(nd);
            break;
        }
        /* 3. update CT info & construct Bags_in/out */
        Bags_in[v_i] = global_dgraph_CT.INs[v_i], Bags_out[v_i] = global_dgraph_CT.OUTs[v_i];
        auto* vi_in = &Bags_in[v_i];
        auto* vi_out = &Bags_out[v_i];
        int in_size = vi_in->size(), out_size = vi_out->size();
        isIntree[v_i] = 1;
        node_order[i] = v_i;
        /* 4. add edge */
        for (int u = 0; u < in_size; u++) {
            results.emplace_back(
                pool.enqueue([u, out_size, vi_in, vi_out] { // pass const type value j to thread; [] can be empty
                    int _u = (*vi_in)[u].first;
                    /* u -> vi -> w */
                    for (int w = 0; w < out_size; w++) {
                        int _w = (*vi_out)[w].first;
                        if (_u != _w) {
                            two_hop_weight_type new_ec = (*vi_in)[u].second + (*vi_out)[w].second;
                            substitute_parallel(_u, _w, new_ec);
                        }
                    }
                    return 1;
                    }));
            results.emplace_back(
                pool.enqueue([u, in_size, vi_in] { // pass const type value j to thread; [] can be empty
                    int _u = (*vi_in)[u].first;
                    /* u -> vi, u2 -> vi */
                    for (int u2 = u + 1; u2 < in_size; u2++) {
                        int _u2 = (*vi_in)[u2].first;
                        substitute_parallel_2(_u, _u2);
                    }
                    return 1;
                    }));
        }
        for (int w = 0; w < out_size; w++) {
            results.emplace_back(pool.enqueue([w, vi_out, out_size] {
                int _w = (*vi_out)[w].first;
                /* vi -> w, vi -> w2 */
                for (int w2 = w + 1; w2 < out_size; w2++) {
                    int _w2 = (*vi_out)[w2].first;
                    substitute_parallel_2(_w, _w2);
                }
                return 1;
                }));  
        }
        for (auto &&result : results) {
            result.get();
        }
        results.clear();
        /* 5. delete edge*/
        for (int u = 0; u < in_size; u++) {
            int _u = (*vi_in)[u].first;
            results.emplace_back(pool.enqueue([_u, v_i] {
                graph_hash_of_mixed_weighted_binary_operations_erase(global_dgraph_CT.OUTs[_u], v_i); // no need to lock here
                return 1;
            }));
        }
        for (int w = 0; w < out_size; w++) {
            int _w = (*vi_out)[w].first;
            results.emplace_back(pool.enqueue([v_i, _w] {
                graph_hash_of_mixed_weighted_binary_operations_erase(global_dgraph_CT.INs[_w], v_i); // no need to lock here
                return 1;
            }));
        }
        std::vector<std::pair<int, two_hop_weight_type>>().swap(global_dgraph_CT.INs[v_i]);
        std::vector<std::pair<int, two_hop_weight_type>>().swap(global_dgraph_CT.OUTs[v_i]);
        for (auto &&result : results) {
            result.get();
        }
        results.clear();
        /* 6. update priority queue */
        for (int u = 0; u < in_size; u++) {
            nd.vertex = (*vi_in)[u].first;
            nd.degree = global_dgraph_CT.degree(nd.vertex);
            q.push(nd);
        }
        for (int w = 0; w < out_size; w++) {
            nd.vertex = (*vi_out)[w].first;
            nd.degree = global_dgraph_CT.degree(nd.vertex);
            q.push(nd);
        }
    }

    case_info.time2_tree_decomposition = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin2).count() / 1e9;
    // --------------------------------------------------------------------------------------------------------

    //--------------------------- step 3: generate CT-tree indexs ---------------------------
    auto begin3 = std::chrono::high_resolution_clock::now();

    /* some variables */
    vector<vector<two_hop_label>> L1_in(N), L1_out(N);
    vector<int> order_mapping(N + 1);
    vector<int> fa(N);
    vector<two_hop_weight_type> temp_dis_in(N);
    vector<two_hop_weight_type> temp_dis_out(N);
    vector<int> index_node_in(N);
    vector<int> index_node_out(N);
    vector<int> islabel_in(N, 0);
    vector<int> islabel_out(N, 0);
    vector<vector<int>> son(N);

    int labelnum = 0;
    root.resize(N);
    dep.resize(N);
    first_pos.resize(N);

    /* order_mapping: origin to merging */
    for (int i = 1; i <= bound_lambda; i++)
        order_mapping[node_order[i]] = i;
    vector<bool> popped_isIncore_node(N, 0);
    /* sort the core node */
    for (int i = bound_lambda + 1; i <= N; i++) {
        struct node_degree nd;
        while (1) {
            nd = q.top();
            q.pop();
            if (!isIntree[nd.vertex] && !popped_isIncore_node[nd.vertex] &&
                global_dgraph_CT.degree(nd.vertex) == nd.degree)
                break;
        }
        node_order[i] = nd.vertex;
        popped_isIncore_node[nd.vertex] = 1;
        order_mapping[nd.vertex] = i;
    }

    /* tree node index */
    for (int i = bound_lambda; i >= 1; i--) {
        /* determine father */
        int v_i = node_order[i]; // merging ID to original ID
        fa[v_i] = INT_MAX;       // merging ID
        int Bags_in_size = Bags_in[v_i].size();
        int Bags_out_size = Bags_out[v_i].size();
        int Bags_size = Bags_in_size + Bags_out_size;
        for (auto it = Bags_in[v_i].begin(); it != Bags_in[v_i].end(); it++)
            if (order_mapping[it->first] < fa[v_i])
                fa[v_i] = order_mapping[it->first];
        for (auto it = Bags_out[v_i].begin(); it != Bags_out[v_i].end(); it++)
            if (order_mapping[it->first] < fa[v_i])
                fa[v_i] = order_mapping[it->first];
        /* interface nodes */
        if (fa[v_i] > bound_lambda || Bags_size == 0) // a root in the forest, interface
        {
            root[v_i] = v_i; // original ID
            fa[v_i] = -1;
            dep[v_i] = 0;
            /* interface nodes only need to add neighborn as hub */
            two_hop_label xx;
            for (auto it = Bags_in[v_i].begin(); it != Bags_in[v_i].end(); it++)
            {
                xx.vertex = it->first;
                xx.distance = it->second;
                L1_in[v_i].push_back(xx);
            }
            for (auto it = Bags_out[v_i].begin(); it != Bags_out[v_i].end(); it++)
            {
                xx.vertex = it->first;
                xx.distance = it->second;
                L1_out[v_i].push_back(xx);
            }
        }
        /* non-root tree nodes */
        else {
            fa[v_i] = node_order[fa[v_i]];
            root[v_i] = root[fa[v_i]];
            dep[v_i] = dep[fa[v_i]] + 1;
            son[fa[v_i]].push_back(v_i);

            int index_node_num_in = 0;
            int index_node_num_out = 0;
            labelnum++;

            int r = root[v_i];

            /* add interface's adj */
            /* here, add into both Lin an Lout is very important */
            for (auto it = Bags_in[r].begin(); it != Bags_in[r].end(); it++) {
                islabel_in[it->first] = labelnum;
                index_node_num_in++;
                index_node_in[index_node_num_in] = it->first;
                temp_dis_in[it->first] = std::numeric_limits<two_hop_weight_type>::max();
                islabel_out[it->first] = labelnum;
                index_node_num_out++;
                index_node_out[index_node_num_out] = it->first;
                temp_dis_out[it->first] = std::numeric_limits<two_hop_weight_type>::max();
            }
            for (auto it = Bags_out[r].begin(); it != Bags_out[r].end(); it++) {
                islabel_in[it->first] = labelnum;
                index_node_num_in++;
                index_node_in[index_node_num_in] = it->first;
                temp_dis_in[it->first] = std::numeric_limits<two_hop_weight_type>::max();
                islabel_out[it->first] = labelnum;
                index_node_num_out++;
                index_node_out[index_node_num_out] = it->first;
                temp_dis_out[it->first] = std::numeric_limits<two_hop_weight_type>::max();
            }

            /* add ancestor */
            int tmp = fa[v_i]; // original id
            while (1)
            {
                if (islabel_in[tmp] != labelnum)
                {
                    islabel_in[tmp] = labelnum;
                    index_node_num_in++;
                    index_node_in[index_node_num_in] = tmp;
                    temp_dis_in[tmp] = std::numeric_limits<two_hop_weight_type>::max();
                }

                if (islabel_out[tmp] != labelnum)
                {
                    islabel_out[tmp] = labelnum;
                    index_node_num_out++;
                    index_node_out[index_node_num_out] = tmp;
                    temp_dis_out[tmp] = std::numeric_limits<two_hop_weight_type>::max();
                }

                if (tmp == r)
                    break;

                tmp = fa[tmp]; // fa[] is original
            }

            /* delta(u) */
            for (auto it = Bags_in[v_i].begin(); it != Bags_in[v_i].end(); it++)
            {
                if (islabel_in[it->first] != labelnum)
                { // line 33
                    islabel_in[it->first] = labelnum;
                    index_node_num_in++;
                    index_node_in[index_node_num_in] = it->first;
                    temp_dis_in[it->first] = it->second;
                }
                else
                { // line ?
                    temp_dis_in[it->first] = it->second;
                }
            }
            /* delta(w) */
            for (auto it = Bags_out[v_i].begin(); it != Bags_out[v_i].end(); it++)
            {
                if (islabel_out[it->first] != labelnum)
                { // line 34
                    islabel_out[it->first] = labelnum;
                    index_node_num_out++;
                    index_node_out[index_node_num_out] = it->first;
                    temp_dis_out[it->first] = it->second;
                }
                else
                { // line ?
                    temp_dis_out[it->first] = it->second;
                }
            }

            /* min vj in N^an-in */
            for (int j = 0; j < Bags_in_size; j++) { // line 33,38
                int v_j = Bags_in[v_i][j].first;
                if (isIntree[v_j])
                {
                    two_hop_weight_type delta_vj = Bags_in[v_i][j].second;
                    int Lj_in_size = L1_in[v_j].size();

                    for (int k = 0; k < Lj_in_size; k++)
                    {
                        int u = L1_in[v_j][k].vertex;

                        two_hop_weight_type delta_uvj = L1_in[v_j][k].distance;
                        if (islabel_in[u] == labelnum && delta_vj + delta_uvj < temp_dis_in[u])
                        {
                            temp_dis_in[u] = delta_vj + delta_uvj;
                        }
                    }
                }
            }
            /* min vj in N^an-out */
            for (int j = 0; j < Bags_out_size; j++) { // line 34,39
                int v_j = Bags_out[v_i][j].first;
                if (isIntree[v_j])
                {
                    two_hop_weight_type delta_vj = Bags_out[v_i][j].second;
                    int Lj_out_size = L1_out[v_j].size();

                    for (int k = 0; k < Lj_out_size; k++)
                    {
                        int w = L1_out[v_j][k].vertex;
                        two_hop_weight_type delta_vjw = L1_out[v_j][k].distance;
                        if (islabel_out[w] == labelnum && delta_vj + delta_vjw < temp_dis_out[w])
                        {
                            temp_dis_out[w] = delta_vj + delta_vjw;
                        }
                    }
                }
            }

            /*add correct indexes of Lines 29-30 of CT into L1, and possibly wrong distances for Lines 31-32 into L1*/
            int delete_num = 0;
            L1_in[v_i].resize(index_node_num_in);
            for (int j = 1; j <= index_node_num_in; j++)
            {
                two_hop_label xx;
                xx.vertex = index_node_in[j];
                xx.distance = temp_dis_in[index_node_in[j]];
                if (xx.distance == std::numeric_limits<two_hop_weight_type>::max())
                {
                    delete_num++;
                    continue;
                }
                L1_in[v_i][j - 1 - delete_num] = xx;
            }
            L1_in[v_i].resize(index_node_num_in - delete_num);

            delete_num = 0;
            L1_out[v_i].resize(index_node_num_out);
            for (int j = 1; j <= index_node_num_out; j++) {
                two_hop_label xx;
                xx.vertex = index_node_out[j];
                xx.distance = temp_dis_out[index_node_out[j]];
                if (xx.distance == std::numeric_limits<two_hop_weight_type>::max())
                {
                    delete_num++;
                    continue;
                }
                L1_out[v_i][j - 1 - delete_num] = xx;
            }
            L1_out[v_i].resize(index_node_num_out - delete_num);
        }
    }

    /*add distance-0 labels to tree nodes; this is needed in querying functions*/
    two_hop_label xx;
    for (int i = 0; i < N; i++) {
        if (isIntree[i]) {
            xx.vertex = i;
            xx.distance = 0;
            L1_in[i].push_back(xx);
            L1_out[i].push_back(xx);
        }
    }

    case_info.time3_tree_indexs = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin3).count() / 1e9;
    //-------------------------------------------------------------------------------------------------------

    //--------------------------- step 4: LCA ---------------------------
    auto begin4 = std::chrono::high_resolution_clock::now();

    /* LCA code; already get the root, the father and the depth, here is the preprocessing of querying LCA */
    int total = 0;

    vector<int> dfn(2 * N + 5);
    for (int i = 1; i <= bound_lambda; i++) {
        int v_i = node_order[i];
        if (root[v_i] == v_i)
            dfs(total, first_pos, v_i, son, dfn);
    }

    if (total > 0) {
        int multi_step = ceil(log(total) / log(2)) + 2;

        tree_st.resize(total + 5);
        tree_st_r.resize(total + 5);
        for (int i = 1; i <= total; i++) {
            tree_st[i].resize(multi_step + 2);
            tree_st_r[i].resize(multi_step + 2);
        }

        vector<int> pow_2(multi_step);

        pow_2[0] = 1;
        for (int i = 1; i < multi_step; i++)
            pow_2[i] = pow_2[i - 1] << 1;

        for (int i = 1; i <= total; i++) {
            tree_st[i][0] = dfn[i];
            tree_st_r[i][0] = dfn[i];
        }

        for (int j = 1; j < multi_step; j++)
            for (int i = 1; i <= total; i++) {
                int k = i + pow_2[j - 1];
                if (k > total || dep[tree_st[i][j - 1]] <= dep[tree_st[k][j - 1]])
                    tree_st[i][j] = tree_st[i][j - 1];
                else
                    tree_st[i][j] = tree_st[k][j - 1];
                k = i - pow_2[j - 1];
                if (k <= 0 || dep[tree_st_r[i][j - 1]] <= dep[tree_st_r[k][j - 1]])
                    tree_st_r[i][j] = tree_st_r[i][j - 1];
                else
                    tree_st_r[i][j] = tree_st_r[k][j - 1];
            }
    }

    lg.resize(total + 1);
    for (int i = 1; i <= total; i++)
        lg[i] = floor(log(i) / log(2));

    /*clear variables not used below*/
    priority_queue<node_degree>().swap(q);
    vector<int>().swap(node_order);
    vector<int>().swap(fa);
    vector<two_hop_weight_type>().swap(temp_dis_in);
    vector<two_hop_weight_type>().swap(temp_dis_out);
    vector<int>().swap(order_mapping);
    vector<int>().swap(islabel_in);
    vector<int>().swap(islabel_out);
    vector<vector<int>>().swap(son);
    vector<int>().swap(index_node_in);
    vector<int>().swap(index_node_out);
    vector<int>().swap(dfn);

    case_info.time4_lca = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin4).count() / 1e9;
    //---------------------------------------------------------------------------------

    //--------------------------- step 5: 2-hop labeling ---------------------------
    auto begin5_1 = std::chrono::high_resolution_clock::now();

    /*remove inf edges*/
    int infty_size = infty_edge.size();
    for (int i = 0; i < infty_size; i++) {
        int v1 = infty_edge[i].first, v2 = infty_edge[i].second;
        results.emplace_back(pool.enqueue([v1,v2] {
            dgraph_v_of_v_remove_edge_infty(v1, v2);
            return 1; }));
    }
    for (auto&& result : results)
        result.get();
    results.clear();

    case_info.time5_core_indexs_prepare1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin5_1).count() / 1e9;
    auto begin5_2 = std::chrono::high_resolution_clock::now();

    /*change IDs*/
    if (case_info.two_hop_order_method == 0) {
        dgraph_change_IDs_sum_IN_OUT_degrees(global_dgraph_CT, new_to_old);
    }
    else if (case_info.two_hop_order_method == 1) {
        dgraph_change_IDs_weighted_degrees(global_dgraph_CT, new_to_old);
    }

    case_info.time5_core_indexs_prepare2 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin5_2).count() / 1e9;
    auto begin5_3 = std::chrono::high_resolution_clock::now();

    /* construct 2-hop labels on core */
    two_hop_case_info_sorted = case_info.two_hop_case_info; // two_hop_case_info_sorted.use_canonical_repair is updated here 
    case_info.pre_core_graph_label_bit_size = compute_CT_label_bit_size(case_info, pool, results); // label bit size part 0
    two_hop_case_info_sorted.max_labal_bit_size = case_info.max_bit_size - case_info.pre_core_graph_label_bit_size;
    if (two_hop_case_info_sorted.max_labal_bit_size < 0) {
        throw reach_limit_error_string_MB;  // after catching error, must call clear_gloval_values_CT and clear CT labels
    }
    two_hop_case_info_sorted.max_run_time_seconds = case_info.max_run_time_seconds - (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin1).count() / 1e9);
    if (two_hop_case_info_sorted.max_run_time_seconds < 0) {
        throw reach_limit_error_string_time;  // after catching error, must call clear_gloval_values_CT and clear CT labels
    }

    case_info.time5_core_indexs_prepare3 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin5_3).count() / 1e9;
    auto begin5_4 = std::chrono::high_resolution_clock::now();

    choose_PLL_PSL(case_info, global_dgraph_CT, pool, results);

    case_info.time5_core_indexs_prepare4 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin5_4).count() / 1e9;
    auto begin5_5 = std::chrono::high_resolution_clock::now();

    if (case_info.use_PLL) {
        dgraph_PLL(global_dgraph_CT, case_info.thread_num, two_hop_case_info_sorted);
    }
    else {
    	dgraph_PSL(global_dgraph_CT, case_info.thread_num, two_hop_case_info_sorted);
    }

    case_info.time5_core_indexs = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin5_5).count() / 1e9;
    auto begin5_6 = std::chrono::high_resolution_clock::now();

    /* return to original ID */
    auto &L_in = case_info.two_hop_case_info.L_in;
    auto &L_out = case_info.two_hop_case_info.L_out; 
    L_in.resize(N);
    L_out.resize(N);
    for (int i = 0; i < N; i++) { /*on Linux, it is faster to make this part parallel*/
        results.emplace_back(pool.enqueue([i, &L_in, &L_out] {
            Label_new_to_old_parallel(i, L_in, L_out);
        return 1;
        }));  
    }
    case_info.two_hop_case_info.label_size_after_canonical_repair = two_hop_case_info_sorted.label_size_after_canonical_repair;
    case_info.two_hop_case_info.label_size_before_canonical_repair = two_hop_case_info_sorted.label_size_before_canonical_repair;
    case_info.two_hop_case_info.time1_PLL_PSL_initialization = two_hop_case_info_sorted.time1_PLL_PSL_initialization;
    case_info.two_hop_case_info.time2_PLL_PSL_label_generation = two_hop_case_info_sorted.time2_PLL_PSL_label_generation;
    case_info.two_hop_case_info.time3_PLL_PSL_label_postprocess = two_hop_case_info_sorted.time3_PLL_PSL_label_postprocess;
    case_info.two_hop_case_info.time4_PLL_PSL_label_canonical = two_hop_case_info_sorted.time4_PLL_PSL_label_canonical;
    case_info.two_hop_case_info.time5_PLL_PSL_total = two_hop_case_info_sorted.time5_PLL_PSL_total;
    for (auto &&result : results) {
        result.get();
    }
    results.clear();
    case_info.core_graph_label_bit_size = case_info.two_hop_case_info.label_size_after_canonical_repair; // label bit size part 1

    case_info.time5_core_indexs_post = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin5_6).count() / 1e9;
    //--------------------------------------------------------------------------------------------------------------------

    //--------------------------- step 6: ---------------------------
    auto begin6 = std::chrono::high_resolution_clock::now();

    /* merge tree_index: L1 into case_info.two_hop_case_info.L */
    for (int v_k = 0; v_k < N; v_k++) {
        if (L1_in[v_k].size() > 0) {
            vector<two_hop_label>(L1_in[v_k]).swap(L1_in[v_k]);
            case_info.two_hop_case_info.L_in[v_k] = L1_in[v_k];
            vector<two_hop_label>().swap(L1_in[v_k]);
        }
        if (L1_out[v_k].size() > 0) {
            vector<two_hop_label>(L1_out[v_k]).swap(L1_out[v_k]);
            case_info.two_hop_case_info.L_out[v_k] = L1_out[v_k];
            vector<two_hop_label>().swap(L1_out[v_k]);
        }
    }

    case_info.time6_post = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin6).count() / 1e9;
    //---------------------------------------------------------------------------------------------------------------------------------

    clear_gloval_values_CT();
    case_info.time_total = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin1).count() / 1e9;
}







/*query function*/

int lca(dgraph_case_info_v2 &case_info, int x, int y) {
    auto &first_pos = case_info.first_pos;

    if (first_pos[x] > first_pos[y]) {
        int t = x;
        x = y;
        y = t;
    }

    int len = first_pos[y] - first_pos[x] + 1;
    int j = case_info.lg[len];
    x = case_info.tree_st[first_pos[x]][j];
    y = case_info.tree_st_r[first_pos[y]][j];
    if (case_info.dep[x] < case_info.dep[y])
        return x;
    else
        return y; // return the vertex with minimum depth between x and y in the dfs sequence.
}

two_hop_weight_type CT_extract_distance(dgraph_case_info_v2 &case_info, int source, int terminal) {
    auto &L_in = case_info.two_hop_case_info.L_in;
    auto &L_out = case_info.two_hop_case_info.L_out;
    auto &Bags_in = case_info.Bags_in;
    auto &Bags_out = case_info.Bags_out;
    auto &isIntree = case_info.isIntree;
    auto &root = case_info.root;

    if (source == terminal)
        return 0;

    two_hop_weight_type distance = std::numeric_limits<two_hop_weight_type>::max();

    if (!isIntree[source] && !isIntree[terminal]) { /* both in core */
        return dgraph_v1_extract_shortest_distance(L_in, L_out, source, terminal);
    }
    else if (!isIntree[source] && isIntree[terminal]) { /* source in core, terminal in tree ; s -> root(t) -> t */
        int r_size = L_in[terminal].size();
        for (int i = 0; i < r_size; i++) {
            if (L_in[terminal][i].distance > distance)
                continue;
            int u = L_in[terminal][i].vertex;
            two_hop_weight_type dist_ut = L_in[terminal][i].distance;
            two_hop_weight_type dist_su = dgraph_v1_extract_shortest_distance(L_in, L_out, source, u);
            if (distance > dist_su + dist_ut)
                distance = dist_su + dist_ut;
        }
        return distance;
    }
    else if (isIntree[source] && !isIntree[terminal]) { /* source in tree, terminal in core ; s -> root(s) -> t */
        int r_size = L_out[source].size();
        for (int i = 0; i < r_size; i++) {
            if (L_out[source][i].distance > distance)
                continue;
            int u = L_out[source][i].vertex;
            two_hop_weight_type dist_su = L_out[source][i].distance;
            two_hop_weight_type dist_ut = dgraph_v1_extract_shortest_distance(L_in, L_out, u, terminal);
            if (distance > dist_su + dist_ut)
                distance = dist_su + dist_ut;
        }
        return distance;
    }
    else if (root[source] != root[terminal]) { /* in different tree */
        int r_s_size = L_out[source].size();
        int r_t_size = L_in[terminal].size();
        for (int i = 0; i < r_s_size; i++) {
            int u = L_out[source][i].vertex;
            double dist_su = L_out[source][i].distance;
            for (int j = 0; j < r_t_size; j++) {
                int w = L_in[terminal][j].vertex;
                double dist_wt = L_in[terminal][j].distance;
                double dist_uw = dgraph_v1_extract_shortest_distance(L_in, L_out, u, w);
                double dis = dist_su + dist_wt + dist_uw;
                if (dis < distance)
                    distance = dis;
            }
        }
        return distance;
    }
    else { /* in same tree */
        int grand = lca(case_info, source, terminal);

        /* d2 : s -> lca -> t */
        unordered_set<int> lca_in = {grand};  // do not use global vectors to accelerate it, for parallel use
        for (int i = Bags_in[grand].size() - 1; i >= 0; i--)
            lca_in.insert(Bags_in[grand][i].first);

        unordered_map<int, two_hop_weight_type> source_dis;
        for (int i = L_out[source].size() - 1; i >= 0; i--) {
            int v = L_out[source][i].vertex;
            if (lca_in.count(v) > 0) {
                source_dis[v] = L_out[source][i].distance;
            }
        }
        for (int i = L_in[terminal].size() - 1; i >= 0; i--) {
            int v = L_in[terminal][i].vertex;
            if (source_dis.count(v) > 0) {
                two_hop_weight_type d2 = source_dis[v] + L_in[terminal][i].distance;
                if (distance > d2) {
                    distance = d2;
                }
            }
        }

        /* d4 : s -> root -> root -> t */
        int r_in_size = L_out[source].size();
        int r_out_size = L_in[terminal].size();
        for (int i = 0; i < r_in_size; i++) {
            int u = L_out[source][i].vertex;
            double dist_su = L_out[source][i].distance;
            for (int j = 0; j < r_out_size; j++) {
                int w = L_in[terminal][j].vertex;
                double dist_wt = L_in[terminal][j].distance;
                double dist_uw = dgraph_v1_extract_shortest_distance(L_in, L_out, u, w);
                double d4 = dist_su + dist_uw + dist_wt;
                if (distance > d4) {
                    distance = d4;
                }
            }
        }
        return distance;
    }
}
