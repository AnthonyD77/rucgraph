#pragma once
#include <map>
#include <build_in_progress/HL/two_hop_v2_old/graph_hash_of_mixed_weighted_CT_v2.h>

/**
 * CT Extension to Path Queries
 *
 * from MLL 3.2.1
 */

int print = 0;

map<pair<int, int>, int> elimination_edge;
vector<int> elimination_vertex;

void clear_gloval_values_CT_extension() {
    graph_hash_of_mixed_weighted_two_hop_clear_global_values();
    graph_v_of_v_idealID().swap(global_ideal_graph_CT);
    map<pair<int, int>, int>().swap(elimination_edge);
    vector<int>().swap(elimination_vertex);
}

void graph_hash_of_mixed_weighted_to_ideal_graph_of_CT_extension(graph_hash_of_mixed_weighted& input_graph, int max_N_ID) {
    global_ideal_graph_CT.resize(max_N_ID);

    for (auto it1 = input_graph.hash_of_vectors.begin(); it1 != input_graph.hash_of_vectors.end(); it1++) {
        int i = it1->first;
        auto search = input_graph.hash_of_hashs.find(i);
        if (search != input_graph.hash_of_hashs.end()) {
            for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
                int j = it2->first;
                if (i < j) {
                    double ec = it2->second;
                    graph_v_of_v_idealID_add_edge(global_ideal_graph_CT, i, j, ec);
                }
            }
        }
        else {
            auto search2 = input_graph.hash_of_vectors.find(i);
            for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
                int j = it2->first;
                if (i < j) {
                    double ec = it2->second;
                    graph_v_of_v_idealID_add_edge(global_ideal_graph_CT, i, j, ec);
                }
            }
        }
    }
}

void substitute_parallel_extension(int v, int v_x, int v_y, double ec) {
    mtx_595[v_x].lock();
    int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(global_ideal_graph_CT[v_x], v_y);
    if (pos != -1) {
        if (ec < global_ideal_graph_CT[v_x][pos].second)  // need to substitute the length
        {
            global_ideal_graph_CT[v_x][pos].second = ec;
            elimination_edge[{v_x, v_y}] = v;
        }
    }
    else {
        graph_hash_of_mixed_weighted_binary_operations_insert(global_ideal_graph_CT[v_x], v_y, ec);
        elimination_edge[{v_x, v_y}] = v;
    }
    mtx_595[v_x].unlock();

    swap(v_x, v_y);

    mtx_595[v_x].lock();
    pos = graph_hash_of_mixed_weighted_binary_operations_search_position(global_ideal_graph_CT[v_x], v_y);
    if (pos != -1) {
        if (ec < global_ideal_graph_CT[v_x][pos].second)  // need to substitute the length
        {
            global_ideal_graph_CT[v_x][pos].second = ec;
            elimination_edge[{v_x, v_y}] = v;
        }
    }
    else {
        pos = graph_hash_of_mixed_weighted_binary_operations_insert(global_ideal_graph_CT[v_x], v_y, ec);
        elimination_edge[{v_x, v_y}] = v;
    }
    mtx_595[v_x].unlock();
}

void remove_parallel_extension(int v_x, int v_y) {
    int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(global_ideal_graph_CT[v_x], v_y);
    global_ideal_graph_CT[v_x].erase(global_ideal_graph_CT[v_x].begin() + pos);
}

void CT_extension(graph_hash_of_mixed_weighted& input_graph, int max_N_ID, graph_hash_of_mixed_weighted_CT_v2_case_info& case_info) {
    //--------------------------------- step 1: initialization ---------------------------
    cout << "step 1: initialization" << endl;
    auto begin1 = std::chrono::high_resolution_clock::now();

    auto& Bags = case_info.Bags;
    auto& isIntree = case_info.isIntree;
    auto& root = case_info.root;
    auto& tree_st = case_info.tree_st;
    auto& tree_st_r = case_info.tree_st_r;
    auto& first_pos = case_info.first_pos;
    auto& lg = case_info.lg;
    auto& dep = case_info.dep;

    int N = input_graph.hash_of_vectors.size();

    // used to extract the path
    graph_hash_of_mixed_weighted_to_ideal_graph_of_CT_extension(input_graph, max_N_ID);
    isIntree.resize(max_N_ID, 0);  // whether it is in the CT-tree

    /*
	priority_queue for maintaining the degrees of vertices  (we do not update degrees in q, so everytime you pop out a degree in q, you check whether it is the right one, ignore it if wrong)

	We need to maintain the degree of vertexes in graph dynamically, otherwise it would increase the cost in time
	*/
    priority_queue<node_degree> q;
    for (int i = 0; i < N; i++)  // change N to max_N_ID cause bugs, why?
    {
        node_degree nd;
        nd.degree = global_ideal_graph_CT[i].size();
        nd.vertex = i;
        q.push(nd);
    }

    auto end1 = std::chrono::high_resolution_clock::now();
    case_info.time_initialization = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9;
    //---------------------------------------------------------------------------------------------------

    //-------------------------------------------------- step 2: MDE-based tree decomposition ------------------------------------------------------------
    cout << "step 2: MDE-based tree decomposition" << endl;
    auto begin2 = std::chrono::high_resolution_clock::now();

    /*MDE-based tree decomposition; generating bags*/
    int bound_lambda = N;
    Bags.resize(N);
    vector<int> right_temp(N);
    vector<int> node_order(N + 1);  // merging ID to original ID

    ThreadPool pool(case_info.thread_num);
    std::vector<std::future<int>> results;  // return typename: xxx

    for (int i = 1; i <= N; i++)  // repeat
    {
        node_degree nd;
        while (1) {
            nd = q.top();
            q.pop();
            if (!isIntree[nd.vertex] && global_ideal_graph_CT[nd.vertex].size() == nd.degree) break;  // nd.vertex is the lowest degree vertex not in tree
        }
        int v_x = nd.vertex;  // the node with the minimum degree in G

        if (nd.degree >= case_info.d)  // reach the boudary
        {
            bound_lambda = i - 1;
            q.push(nd);
            case_info.tree_vertex_num = i - 1;
            break;  // until |Ni| >= d
        }

        elimination_vertex.push_back(v_x);
        isIntree[v_x] = 1;  // add to CT-tree
        node_order[i] = v_x;

        vector<pair<int, double>>& adj_temp = global_ideal_graph_CT[v_x];  // global_ideal_graph_CT is G_i-1
        int v_adj_size = global_ideal_graph_CT[v_x].size();
        for (int j = 0; j < v_adj_size; j++) {
            auto it = &global_ideal_graph_CT[v_x][j];
            Bags[v_x].push_back({it->first, it->second});  // Bags[v_x] stores adj vertices and weights of v_x
        }

        /*add new edge*/
        for (int j = 0; j < v_adj_size; j++) {
            int adj_j = adj_temp[j].first;
            int right_j = right_temp[j];
            for (int k = j + 1; k < v_adj_size; k++) {
                int adj_k = adj_temp[k].first;
                double new_ec = adj_temp[j].second + adj_temp[k].second;
                if (v_x == 9 && adj_j == 0 && adj_k == 7) {
                    cout << "new_ec = " << new_ec << endl;
                }
                int right_k = right_temp[k];
                results.emplace_back(
                    pool.enqueue([v_x, adj_j, adj_k, new_ec] {  // pass const type value j to thread; [] can be empty
                        substitute_parallel_extension(v_x, adj_j, adj_k, new_ec);
                        return 1;
                    }));
            }
        }
        for (auto&& result : results) {
            result.get();
        }
        results.clear();

        // delete edge related to v_x and update degree;
        for (int j = 0; j < v_adj_size; j++) {
            //update degree
            nd.vertex = adj_temp[j].first;
            nd.degree = global_ideal_graph_CT[nd.vertex].size() - 1;  // v_x will be removed from global_ideal_graph_CT[nd.vertex] below
            q.push(nd);
            // remove v_x
            int m = global_ideal_graph_CT[v_x][j].first;
            results.emplace_back(
                pool.enqueue([m, v_x] {  // pass const type value j to thread; [] can be empty
                    remove_parallel_extension(m, v_x);
                    return 1;
                }));
        }
        for (auto&& result : results) {
            result.get();
        }
        results.clear();
        // delete v_x from ideal graph directly
        vector<pair<int, double>>().swap(global_ideal_graph_CT[v_x]);
    }

    auto end2 = std::chrono::high_resolution_clock::now();
    case_info.time_tree_decomposition = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2).count() / 1e9;
    //--------------------------------------------------------------------------------------------------------

    //---------------------------------------------- step 3: generate CT-tree indexs ------------------------------------------------
    cout << "step 3: generate CT-tree indexs" << endl;
    auto begin3 = std::chrono::high_resolution_clock::now();

    /* generate CT-tree indexs */
    vector<vector<two_hop_label_v1>> L1(N);  // Labels for CT-tree index merge later, otherwise may increase query time on PLL
    vector<int> fa(N);
    root.resize(N);
    vector<double> temp_dis(N);
    vector<int> order_mapping(N + 1);  // original ID to merging ID
    for (int i = 1; i <= bound_lambda; i++) order_mapping[node_order[i]] = i;
    vector<bool> popped_isIncore_node(N, 0);

    for (int i = bound_lambda + 1; i <= N; i++)  // take advantage of the priority queue
    {
        struct node_degree nd;
        while (1) {
            nd = q.top();
            q.pop();
            if (!isIntree[nd.vertex] && !popped_isIncore_node[nd.vertex] && global_ideal_graph_CT[nd.vertex].size() == nd.degree) break;
        }
        node_order[i] = nd.vertex;  // merging ID to original ID
        popped_isIncore_node[nd.vertex] = 1;
        order_mapping[nd.vertex] = i;  // original ID to merging ID
    }

    dep.resize(N);
    vector<int> islabel(N, 0);
    first_pos.resize(N);
    vector<vector<int>> son(N);
    vector<int> index_node(N);

    vector<int> isneighour(N);
    int neighbournum = 0;
    vector<double> T_temp_n(N);

    int labelnum = 0;

    for (int i = bound_lambda; i >= 1; i--) {
        int v_x = node_order[i];  //  node_order(N + 1); // merging ID to original ID

        fa[v_x] = INT_MAX;  // merging ID
        int v_adj_size = Bags[v_x].size();

        for (int j = 0; j < v_adj_size; j++) {
            // renew fa[v_x] to be the smallest merging_ID (the lowest ancestor bag)
            if (order_mapping[Bags[v_x][j].first] < fa[v_x]) {
                fa[v_x] = order_mapping[Bags[v_x][j].first];
            }
        }

        if (fa[v_x] > bound_lambda || v_adj_size == 0) {
            // a root in the forest (a bag with interfaces)
            root[v_x] = v_x;
            fa[v_x] = -1;
            dep[v_x] = 0;

            two_hop_label_v1 xx;
            for (int j = 0; j < v_adj_size; j++)  // (bound_lambda - 1) - local distance to the interface
            {
                xx.vertex = Bags[v_x][j].first;
                xx.distance = Bags[v_x][j].second;
                xx.parent_vertex = -1;
                L1[v_x].push_back(xx);  // interpfact parts of tree indexes of root v_x
            }
        }
        else {
            //a non-root in the forest
            fa[v_x] = node_order[fa[v_x]];  //  node_order(N + 1); // merging ID to original ID
            root[v_x] = root[fa[v_x]];      // int i = bound_lambda; i >= 1; i--, already has the right root[fa[v_x]];  root[v_x] = root[fa[v_x]], from high to low to get roots
            dep[v_x] = dep[fa[v_x]] + 1;    // multi_fa[v_x][0] = fa[v_x]; for LCA
            son[fa[v_x]].push_back(v_x);    // for LCA

            int index_node_num = 0;
            labelnum++;

            int root_adj_size = Bags[root[v_x]].size();  // interfaces, already added above

            /*add interface node*/
            for (int j = 0; j < root_adj_size; j++)            // put the labels to interface in the beginnig position
            {                                                  // add interface
                islabel[Bags[root[v_x]][j].first] = labelnum;  // labelnum means that a hub vertex is added into bag v_x
                index_node_num++;
                index_node[index_node_num] = Bags[root[v_x]][j].first;
                temp_dis[Bags[root[v_x]][j].first] = std::numeric_limits<double>::max();  // initial dis
            }

            /*add ancestor node: representation nodes of all ancetor bags*/
            int v_y = v_x;
            while (fa[v_y] != -1) {
                // add ancestor
                if (islabel[fa[v_y]] != labelnum)  // fa[v_y] is not in bag v_x yet
                {
                    index_node_num++;
                    index_node[index_node_num] = fa[v_y];
                    islabel[fa[v_y]] = labelnum;                             // add fa[v_y] into bag v_x
                    temp_dis[fa[v_y]] = std::numeric_limits<double>::max();  // initial dis
                }
                v_y = fa[v_y];
            }

            /*corresponds to Line 30 of CT: the first value after min: delta_u = Bags[v_x][j].second*/
            for (int j = 0; j < v_adj_size; j++) {
                // add neighbours
                if (islabel[Bags[v_x][j].first] != labelnum) {
                    islabel[Bags[v_x][j].first] = labelnum;
                    index_node.push_back(Bags[v_x][j].first);
                    temp_dis[Bags[v_x][j].first] = Bags[v_x][j].second;
                }
                else {
                    temp_dis[Bags[v_x][j].first] = Bags[v_x][j].second;
                }
            }
            // query (bound_lambda - 1)_local_distance to ancestor or the interface through neighbours

            /*corresponds to Line 30 of CT: the second min value: dis_vj + L1[vj][k].distance*/
            for (int j = 0; j < v_adj_size; j++) {
                if (isIntree[Bags[v_x][j].first])  // isneighour and isintree --> can be used as an intermediate node to update labels
                {
                    double dis_vj = Bags[v_x][j].second;
                    int vj = Bags[v_x][j].first;
                    int Lj_size = L1[vj].size();
                    for (int k = 0; k < Lj_size; k++)  // update the (bound_lambda-1)_local_distance
                    {
                        if (islabel[L1[vj][k].vertex] == labelnum && dis_vj + L1[vj][k].distance < temp_dis[L1[vj][k].vertex]) {
                            temp_dis[L1[vj][k].vertex] = dis_vj + L1[vj][k].distance;
                        }
                    }
                }
            }

            /*add correct indexes of Lines 29-30 of CT into L1, and possibly wrong distances for Lines 31-32 into L1*/
            // add labels to L1
            L1[v_x].resize(index_node_num);
            for (int j = 1; j <= index_node_num; j++) {
                two_hop_label_v1 xx;
                xx.vertex = index_node[j];
                xx.distance = temp_dis[index_node[j]];
                xx.parent_vertex = -1;
                L1[v_x][j - 1] = xx;
            }

            /*Lines 31-32 of CT; update possibly wrong distances for Lines 31-32 in L1*/
            // update conversely
            neighbournum++;
            for (int j = 0; j < v_adj_size; j++) {
                isneighour[Bags[v_x][j].first] = neighbournum;
                T_temp_n[Bags[v_x][j].first] = Bags[v_x][j].second;
            }
            for (int j = 1; j <= index_node_num; j++) {
                int vj = index_node[j];
                int Lj_size = L1[vj].size();
                for (int k = 0; k < Lj_size; k++) {
                    int vk = L1[vj][k].vertex;
                    if ((isneighour[vk] == neighbournum) && (T_temp_n[vk] + L1[vj][k].distance < temp_dis[vj])) {
                        temp_dis[vj] = T_temp_n[vk] + L1[vj][k].distance;
                    }
                }
            }
            for (int j = 1; j <= index_node_num; j++) {
                if (temp_dis[index_node[j]] < L1[v_x][j - 1].distance) {
                    L1[v_x][j - 1].distance = temp_dis[index_node[j]];
                    L1[v_x][j - 1].parent_vertex = -1;
                }
            }
        }
    }

    /* add distance-0 labels to tree nodes; this is needed in querying functions */
    two_hop_label_v1 node;
    for (int i = 0; i < N; i++) {
        if (isIntree[i]) {
            node.vertex = i;
            node.distance = 0;
            node.parent_vertex = -1;
            L1[i].push_back(node);
        }
    }

    auto end3 = std::chrono::high_resolution_clock::now();
    case_info.time_tree_indexs = std::chrono::duration_cast<std::chrono::nanoseconds>(end3 - begin3).count() / 1e9;
    //-------------------------------------------------------------------------------------------------------

    //------------------------------------------------ step 4: LCA --------------------------------------------------------------
    cout << "step 4: LCA" << endl;
    auto begin4 = std::chrono::high_resolution_clock::now();

    /* LCA code; already get the root, the father and the depth, here is the preprocessing of querying LCA */
    int total = 0;
    vector<int> dfn(2 * N + 5);
    for (int i = 1; i <= bound_lambda; i++) {
        int v_x = node_order[i];
        if (root[v_x] == v_x) dfs(total, first_pos, v_x, son, dfn);
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
    vector<double>().swap(temp_dis);
    vector<int>().swap(order_mapping);
    vector<int>().swap(islabel);
    vector<vector<int>>().swap(son);
    vector<int>().swap(index_node);
    vector<int>().swap(isneighour);
    vector<double>().swap(T_temp_n);
    vector<int>().swap(dfn);
    vector<int>().swap(right_temp);
    vector<int>().swap(fa);

    auto end4 = std::chrono::high_resolution_clock::now();
    case_info.time_lca = std::chrono::duration_cast<std::chrono::nanoseconds>(end4 - begin4).count() / 1e9;
    //--------------------------------------------------------------------------------------------------------------------------

    //----------------------------------------------- step 5: 2-hop labeling -------------------------------------------
    cout << "step 5: 2-hop labeling" << endl;
    auto begin5 = std::chrono::high_resolution_clock::now();

    /*update limits*/
    double to_date_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end4 - begin1).count() / 1e9;
    case_info.two_hop_case_info.max_run_time_seconds = case_info.max_run_time_seconds - to_date_time;
    if (case_info.two_hop_case_info.max_run_time_seconds < 0) {
        throw reach_limit_error_string_time;
    }
    long long int to_date_bit_size = case_info.compute_label_bit_size();
    case_info.two_hop_case_info.max_labal_size = (case_info.max_bit_size - to_date_bit_size) / sizeof(two_hop_label_v1);  // this is slightly inaccurate, since reduction measures of R1 R2 are not counted
    if (case_info.two_hop_case_info.max_labal_size < 0) {
        throw reach_limit_error_string_MB;
    }

    /* construct 2-hop labels on core */
    /* graph_v_of_v_idealID to graph_hash_of_mixed_weighted */
    graph_hash_of_mixed_weighted hash_g;
    int size = global_ideal_graph_CT.size();
    for (int i = 0; i < size; i++) {
        int v_size = global_ideal_graph_CT[i].size();
        for (int x = 0; x < v_size; x++) {
            int j = global_ideal_graph_CT[i][x].first;
            if (i < j) {
                double ec = global_ideal_graph_CT[i][x].second;
                graph_hash_of_mixed_weighted_add_edge(hash_g, i, j, ec);
            }
        }
    }
    case_info.core_graph = hash_g;
    if (case_info.use_PLL == 1) {
        graph_hash_of_mixed_weighted_PLL_v1(hash_g, max_N_ID, 1, case_info.thread_num, case_info.two_hop_case_info);
    }
    else if (case_info.use_PLL == 0) {
        graph_hash_of_mixed_weighted_PSL_v1(hash_g, max_N_ID, case_info.thread_num, case_info.two_hop_case_info);
    }
    else if (case_info.use_PLL == -1) {
        VCPLL(hash_g, max_N_ID, 1, case_info.thread_num, case_info.two_hop_case_info);
    }
    auto end5 = std::chrono::high_resolution_clock::now();
    case_info.time_core_indexs = std::chrono::duration_cast<std::chrono::nanoseconds>(end5 - begin5).count() / 1e9;
    //--------------------------------------------------------------------------------------------------------------------

    //-------------------------------------------------- step 6: postprocessing -------------------------------------------------------------------
    cout << "step 6: postprocessing" << endl;
    auto begin6 = std::chrono::high_resolution_clock::now();

    /* merge tree_index: L1 into case_info.two_hop_case_info.L */
    for (int v_k = 0; v_k < N; v_k++) {
        if (L1[v_k].size() > 0) {
            vector<two_hop_label_v1>(L1[v_k]).swap(L1[v_k]);
            case_info.two_hop_case_info.L[v_k] = L1[v_k];
            vector<two_hop_label_v1>().swap(L1[v_k]);
        }
    }

    /* Extention: generate tree pred */
    auto begin_extension = std::chrono::high_resolution_clock::now();

    auto& L_ = case_info.two_hop_case_info.L;
    cout << "step 6.1: extention" << endl;

    cout << "elimination_vertex" << endl;
    for (auto it : elimination_vertex) {
        cout << it << ", ";
    }
    cout << endl;

    cout << "elimination_vertex" << endl;
    for (auto it : elimination_edge) {
        cout << it.first.first << ", " << it.first.second << " -> " << it.second << endl;
    }

    int cc = 0;
    for (auto it: Bags){
        cout << "Bags " << cc << ": ";
        for (auto itt : it) {
            cout << itt.first << " ";
        }
        cout << endl;
        cc++;
    }


    int u, v = 0, Lu_size = 0, v_in_Bagu = 0;
    int tree_size = elimination_vertex.size();
    for (int j = 0; j < tree_size; j++) {
        u = elimination_vertex[j];
        Lu_size = L_[u].size();
        for (int i = 0; i < Lu_size; i++) {
            // check if v is in Bag(u)
            v = L_[u][i].vertex;
            if (v == u) {
                continue;
            }
            v_in_Bagu = 0;
            for (auto it : Bags[u]) {
                // v is in Bag(u)
                if (it.first == v) {
                    double edge_uv = graph_hash_of_mixed_weighted_edge_weight(input_graph, u, v);
                    if ((edge_uv - L_[u][i].distance) < 1e-5){
                        // nothing
                    }
                    else {
                        for (auto it : elimination_vertex) {
                            if (it == u || it == v) {
                                continue;
                            }
                            double dist1 = CT_extract_distance(case_info, it, u);
                            double dist2 = CT_extract_distance(case_info, it, v);
                            if (dist1 != std::numeric_limits<double>::max() && dist2 != std::numeric_limits<double>::max() && abs(dist1 + dist2 - L_[u][i].distance) < 1e-5) {
                                L_[u][i].parent_vertex = it;
                                break;
                            }
                        }
                    }
                    v_in_Bagu = 1;
                    break;
                }
            }
            // v is not in Bag(u)
            if (!v_in_Bagu) {
                for (auto it : Bags[u]) {
                    double dist1 = CT_extract_distance(case_info, it.first, u);
                    double dist2 = CT_extract_distance(case_info, it.first, v);
                    if (u == 0) {
                        cout << it.first << ": " << "dist1=" << dist1 << ", dist2=" << dist2 << endl;
                    }
                    if (dist1 != std::numeric_limits<double>::max() && dist2 != std::numeric_limits<double>::max() && abs(dist1 + dist2 - L_[u][i].distance) < 1e-5) {
                        L_[u][i].parent_vertex = it.first;
                        break;
                    }
                }
            }
        }
    }

    auto end_extension = std::chrono::high_resolution_clock::now();
    case_info.time_extension = std::chrono::duration_cast<std::chrono::nanoseconds>(end_extension - begin_extension).count() / 1e9;
    /* extension end */

    for (int v_k = 0; v_k < N; v_k++) {
        if (isIntree[v_k]) {
            continue;
        }
        int Lvk_size = L_[v_k].size();
        for (int j = 0; j < Lvk_size; j++) {
            while (1) {
                if (elimination_edge.find({v_k, L_[v_k][j].parent_vertex}) == elimination_edge.end()) {
                    break;
                }
                L_[v_k][j].parent_vertex = elimination_edge[{v_k, L_[v_k][j].parent_vertex}];
            }
        }
    }

    auto end6 = std::chrono::high_resolution_clock::now();
    case_info.time_post = std::chrono::duration_cast<std::chrono::nanoseconds>(end6 - begin6).count() / 1e9;
    //---------------------------------------------------------------------------------------------------------------------------------

    clear_gloval_values_CT_extension();

    case_info.time_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end6 - begin1).count() / 1e9;
}

/**
 * Get the path using extension parent_vertex
 */
vector<pair<int, int>> query_path_sp1(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_CT_v2_case_info& case_info, int source, int terminal) {
    if (source == terminal) {
        return vector<pair<int, int>>();
    }

    auto& L = case_info.two_hop_case_info.L;
    int s_vertex, t_vertex;
    int s_size = L[source].size(), t_size = L[terminal].size();
    int w = -2, sw_, tw_;
    double distance = std::numeric_limits<double>::max(), new_dis;

    for (int i = 0; i < s_size; i++) {
        s_vertex = L[source][i].vertex;
        for (int j = 0; j < t_size; j++) {
            t_vertex = L[terminal][j].vertex;
            if (s_vertex == t_vertex) {
                new_dis = L[source][i].distance + L[terminal][j].distance;
                if (new_dis < distance) {
                    distance = new_dis;
                    w = s_vertex;
                    sw_ = L[source][i].parent_vertex;
                    tw_ = L[terminal][j].parent_vertex;
                }
            }
        }
    }

    //    cout <<"\t sp1 w=" << w << " sw_=" << sw_ << " tw_=" << tw_ <<
    //        " s=" << source << " t=" << terminal <<endl;

    if (w == -2) {
        vector<pair<int, int>> rslt;
        rslt.push_back({INT_MAX, INT_MAX});
        return rslt;
    }

    // path(s,w)
    vector<pair<int, int>> path_sw;
    if (source == w) {
        // nothing
    }
    else {
        double edge_sw = graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, w);
        if (edge_sw != std::numeric_limits<double>::max() && (sw_ == -1 || sw_ == w)) {
            path_sw.push_back({source, w});
        }
        else {
            vector<pair<int, int>> p1 = query_path_sp1(instance_graph, case_info, source, sw_);
            vector<pair<int, int>> p2 = query_path_sp1(instance_graph, case_info, sw_, w);
            path_sw.insert(path_sw.end(), p1.begin(), p1.end());
            path_sw.insert(path_sw.end(), p2.begin(), p2.end());
        }
    }

    // path(t,w)
    vector<pair<int, int>> path_tw;
    if (terminal == w) {
        // nothing
    }
    else {
        double edge_tw = graph_hash_of_mixed_weighted_edge_weight(instance_graph, terminal, w);
        if (edge_tw != std::numeric_limits<double>::max() && (tw_ == -1 || tw_ == w)) {
            path_tw.push_back({terminal, w});
        }
        else {
            vector<pair<int, int>> p1 = query_path_sp1(instance_graph, case_info, terminal, tw_);
            vector<pair<int, int>> p2 = query_path_sp1(instance_graph, case_info, tw_, w);
            path_tw.insert(path_tw.end(), p1.begin(), p1.end());
            path_tw.insert(path_tw.end(), p2.begin(), p2.end());
        }
    }

    path_sw.insert(path_sw.end(), path_tw.begin(), path_tw.end());
    return path_sw;
}

vector<pair<int, int>> query_path_sp2(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_CT_v2_case_info& case_info, int source, int terminal, double distance = 0) {
    // s is in tree, t is in core, and t is a hub of s

    if (source == terminal) {
        return vector<pair<int, int>>();
    }

    double edge_st = graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, terminal);
    if (edge_st != std::numeric_limits<double>::max() && abs(edge_st - distance) < 1e-5) {
        vector<pair<int, int>> rslt;
        rslt.push_back({source, terminal});
        return rslt;
    }

    auto& L = case_info.two_hop_case_info.L;
    int s_ptr = 0, s_size = L[source].size();
    int w = -2;

    while (s_ptr < s_size) {
        if (L[source][s_ptr].vertex == terminal) {
            w = L[source][s_ptr].parent_vertex;
            break;
        }
        s_ptr++;
    }

//        cout << "\t sp2 w=" << w << endl;

    if (w == -2) {
        return vector<pair<int, int>>();
    }

    // path(s,t)
    vector<pair<int, int>> path;
    double dist_st = graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, terminal);
    if ((w == -1 || w == terminal) && dist_st != std::numeric_limits<double>::max()) {
        path.push_back({source, terminal});
    }
    else {
        vector<pair<int, int>> p1 = query_path_sp1(instance_graph, case_info, source, w);
        vector<pair<int, int>> p2 = query_path_sp1(instance_graph, case_info, w, terminal);
        path.insert(path.end(), p1.begin(), p1.end());
        path.insert(path.end(), p2.begin(), p2.end());
    }

    return path;
}

void CT_extension_extract_path(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_CT_v2_case_info& case_info, int source, int terminal, vector<pair<int, int>>& path) {
    /* return empty vector */
    if (source == terminal) {
        return;
    }

    auto& root = case_info.root;
    auto& isIntree = case_info.isIntree;
    auto& Bags = case_info.Bags;
    auto& L = case_info.two_hop_case_info.L;

    int s_size = L[source].size();
    int t_size = L[terminal].size();
    double distance = std::numeric_limits<double>::max();

    if (!isIntree[source] && !isIntree[terminal]) {
        /* both in core */

        if (print) { cout << "case 1" << endl; }
        pair<int, int> two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, source, terminal);

        if (two_predecessors.first == source && two_predecessors.second == terminal) {  // disconnected
            // disconnect
            path.push_back({INT_MAX, INT_MAX});
            return;
        }

        int w = -2;
        if (source != two_predecessors.first) {
            int s_next = two_predecessors.first;
            double edge = graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, s_next);
            if (edge != std::numeric_limits<double>::max()) {
                path.push_back({source, s_next});
            }
            else {
                // need unfold the edge
                for (int i = 0; i < s_size; i++) {
                    if (L[source][i].vertex == s_next) {
                        w = L[source][i].parent_vertex;
                        break;
                    }
                }
                if (w == -2) {
                    int s_next_size = L[s_next].size();
                    for (int i = 0; i < s_next_size; i++) {
                        if (L[s_next][i].vertex == source) {
                            w = L[s_next][i].parent_vertex;
                            break;
                        }
                    }
                }
                vector<pair<int, int>> unfold_path1 = query_path_sp2(instance_graph, case_info, w, source);
                vector<pair<int, int>> unfold_path2 = query_path_sp2(instance_graph, case_info, w, s_next);
                path.insert(path.end(), unfold_path1.begin(), unfold_path1.end());
                path.insert(path.end(), unfold_path2.begin(), unfold_path2.end());
            }
        }

        if (terminal != two_predecessors.second) {
            int t_next = two_predecessors.second;
            double edge = graph_hash_of_mixed_weighted_edge_weight(instance_graph, terminal, t_next);
            if (edge != std::numeric_limits<double>::max()) {
                path.push_back({terminal, t_next});
            }
            else {
                // need unfold the edge
                for (int i = 0; i < t_size; i++) {
                    if (L[terminal][i].vertex == t_next) {
                        w = L[terminal][i].parent_vertex;
                        break;
                    }
                }
                if (w == -2) {
                    int t_next_size = L[t_next].size();
                    for (int i = 0; i < t_next_size; i++) {
                        if (L[t_next][i].vertex == source) {
                            w = L[t_next][i].parent_vertex;
                            break;
                        }
                    }
                }
                vector<pair<int, int>> unfold_path1 = query_path_sp2(instance_graph, case_info, w, terminal);
                vector<pair<int, int>> unfold_path2 = query_path_sp2(instance_graph, case_info, w, t_next);
                path.insert(path.end(), unfold_path1.begin(), unfold_path1.end());
                path.insert(path.end(), unfold_path2.begin(), unfold_path2.end());
            }
        }

        CT_extension_extract_path(instance_graph, case_info, two_predecessors.first, two_predecessors.second, path);
        return;
    }
    else if (!isIntree[source] && isIntree[terminal]) {
        /* source is in core, terminal is in tree */

        if (print){ cout << "case 2" << endl; }
        int r = root[terminal];
        int r_size = Bags[r].size();
        int interface_c = INT_MAX;
        double interface_distance = 0;

        for (int i = 0; i < r_size; i++) {
            if (L[terminal][i].distance > distance)
                continue;

            // x is the interface vertex
            int x = L[terminal][i].vertex;
            double x_dis = L[terminal][i].distance;

            double dis = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, source, x);
            if (distance > dis + x_dis) {
                distance = dis + x_dis;
                interface_c = x;
                interface_distance = x_dis;
            }
        }

        if (print) { cout << "\t interface c: " << interface_c << " distance: " << distance << endl; }
        if (interface_c == INT_MAX) {
            path.push_back({INT_MAX, INT_MAX});
            return;
        }

        // path(t,c)
        vector<pair<int, int>> path_tc = query_path_sp2(instance_graph, case_info, terminal, interface_c, interface_distance);
        // path(s,c)
        CT_extension_extract_path(instance_graph, case_info, source, interface_c, path);

        path.insert(path.end(), path_tc.begin(), path_tc.end());

        return;
    }
    else if (isIntree[source] && !isIntree[terminal]) {
        // source is in tree

        if (print) { cout << "case 3" << endl; }
        return CT_extension_extract_path(instance_graph, case_info, terminal, source, path);
    }
    else {
        // both in tree

        if (print) { cout << "case 4" << endl; }
        for (int i = 0; i < s_size; i++) {
            for (int j = 0; j < t_size; j++) {
                if (L[source][i].vertex == L[terminal][j].vertex) {
                    double new_dis = L[source][i].distance + L[terminal][j].distance;
                    if (distance > new_dis) {
                        distance = new_dis;
                    }
                }
            }
        }

        distance += 1e-5;
        int c = -2, d = -2;
        double dist_sc = 0, dist_cd = 0, dist_dt = 0;
        for (int i = 0; i < s_size; i++) {
            int c_ = L[source][i].vertex;
            if (!isIntree[c_]) {
                dist_sc = L[source][i].distance;
                for (int j = 0; j < t_size; j++) {
                    int d_ = L[terminal][j].vertex;
                    if (!isIntree[d_]) {
                        dist_dt = L[terminal][j].distance;
                        dist_cd = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, c_, d_);
                        if (dist_cd == std::numeric_limits<double>::max()) {
                            continue;
                        }
                        if (distance >= (dist_cd + dist_sc + dist_dt)) {
                            distance = (dist_cd + dist_sc + dist_dt);
                            c = c_;
                            d = d_;
                        }
                    }
                }
            }
        }

        if (c == -2 && d == -2) {
            vector<pair<int, int>> rslt = query_path_sp1(instance_graph, case_info, source, terminal);
            path.insert(path.end(), rslt.begin(), rslt.end());
        }
        else {
            CT_extension_extract_path(instance_graph, case_info, source, c, path);
            CT_extension_extract_path(instance_graph, case_info, c, d, path);
            CT_extension_extract_path(instance_graph, case_info, d, terminal, path);
        }
        return;
    }
}
