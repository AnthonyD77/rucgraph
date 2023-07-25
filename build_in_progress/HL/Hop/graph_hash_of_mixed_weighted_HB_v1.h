#pragma once
#include "graph_v_of_v_idealID/graph_v_of_v_idealID.h"
#include <iostream>
#include <limits>
#include <numeric>
#include <tool_functions/ThreadPool.h>
#include <shared_mutex>
#include <chrono>
#include <map>
#include <boost/heap/fibonacci_heap.hpp>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <graph_hash_of_mixed_weighted/two_graphs_operations/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID.h>
#include <build_in_progress/HL/Hop/graph_hash_of_mixed_weighted_two_hop_labels_v1.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_update_vertexIDs.h>

using namespace std;

struct HBBPS_v1_node
{
  public:
    int vertex, parent_vertex;
};

struct HBPLL_v1_node {
public:
	int vertex, parent_vertex, hop;
	double priority_value;
}; // define the node in the queue
bool operator<(HBPLL_v1_node const& x, HBPLL_v1_node const& y) {
	return x.priority_value > y.priority_value; // < is the max-heap; > is the min heap
}
typedef typename boost::heap::fibonacci_heap<HBPLL_v1_node>::handle_type graph_hash_of_mixed_weighted_HL_PLL_v1_handle_t_for_sp;

void graph_hash_of_mixed_weighted_HL_HB_v1_thread_function_HBBFS(int v_k, int N, int upper_k)
{
    /* Pruned Dijkstra from vertex v_k; v_k = u in HBPLL */

    mtx_599[max_N_ID_for_mtx_599 - 1].lock();
    int used_id = Qid_599.front();
    Qid_599.pop();
    mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

    queue<int> T_changed_vertices;
    vector<vector<double>> d;

    /* line 3 */
    HBBPS_v1_node node;
    vector<vector<HBBPS_v1_node>> Q;
    Q.resize(1);
    node.vertex = v_k;
    node.parent_vertex = v_k;
    Q[0].push_back(node);

    /* line 4, 5 */
    vector<double> d_0;
    d_0.resize(N);
    for (int i = 0; i < N; i++)
        d_0[i] = (i == v_k) ? 0 : std::numeric_limits<double>::max();
    d.push_back(d_0);

    /* T数组加速 */
    mtx_599[v_k].lock();
    int L_v_k_size = L_temp_599[v_k].size();
    for (int i = 0; i < L_v_k_size; i++)
    {
        int L_v_k_i_vertex = L_temp_599[v_k][i].vertex;
        T_bfs_599[used_id][L_v_k_i_vertex].first = L_temp_599[v_k][i].distance;
        T_bfs_599[used_id][L_v_k_i_vertex].second = L_temp_599[v_k][i].hop;
        T_changed_vertices.push(L_v_k_i_vertex);
    }
    mtx_599[v_k].unlock();

    long long int new_label_num = 0;
    two_hop_label_v1 xx;

    int h = 0;
    int upper_h = upper_k == 0 ? std::numeric_limits<int>::max() : upper_k;

    while (1) // line 6
    {
        if (h > upper_h || Q[h].empty()) // line 7, 8
            break;

        Q.push_back(vector<HBBPS_v1_node>()); // line 9
        d.push_back(d[h]);                    // line 10

        int Q_h_size = Q[h].size();
        for (int ii = 0; ii < Q_h_size; ii++) // line 11
        {
            node = Q[h][ii];
            int v = node.vertex;
            int v_parent = node.parent_vertex;

            if (v_k > v) // line 12, 13: rank pruning
                continue;

            /* line 14: query(u,v,h) */
            /* Under Windows, lock should be put in the for loop */
            double query_vk_v = std::numeric_limits<double>::max();
            mtx_599[v].lock();
            auto L_u_size1 = L_temp_599[v].size();
            for (int i = 0; i < L_u_size1; i++)
            {
                if (L_temp_599[v][i].hop + T_bfs_599[used_id][L_temp_599[v][i].vertex].second <= h)
                {
                    double dis = L_temp_599[v][i].distance + T_bfs_599[used_id][L_temp_599[v][i].vertex].first;
                    if (query_vk_v > dis)
                        query_vk_v = dis;
                }
            }
            mtx_599[v].unlock();

            /* line 14, 15 */
            if (d[h][v] + 1e-5 >= query_vk_v)
                continue;

            /* line 16 */
            xx.vertex = v_k;
            xx.distance = d[h][v];
            xx.parent_vertex = v_parent;
            xx.hop = h;
            mtx_599[v].lock();
            L_temp_599[v].push_back(xx);
            mtx_599[v].unlock();
            new_label_num++;

            int v_adj_size = ideal_graph_599[v].size();
            for (int i = 0; i < v_adj_size; i++) // line 17
            {
                int adj_v = ideal_graph_599[v][i].first;
                double ec = ideal_graph_599[v][i].second;
                if (d[h][v] + ec < d[h + 1][adj_v]) // line 18
                {
                    node.vertex = adj_v;
                    node.parent_vertex = v;
                    Q[h + 1].push_back(node);       // line 19
                    d[h + 1][adj_v] = d[h][v] + ec; // line 20
                }
            }
        }
        h++;
    }

    while (T_changed_vertices.size() > 0)
    {
        T_bfs_599[used_id][T_changed_vertices.front()].first = std::numeric_limits<double>::max(); // reverse-allocate T values
        T_changed_vertices.pop();
    }

    mtx_599[v_k].lock();
    vector<two_hop_label_v1>(L_temp_599[v_k]).swap(L_temp_599[v_k]);
    mtx_599[v_k].unlock();

    mtx_599[max_N_ID_for_mtx_599 - 1].lock();
    Qid_599.push(used_id);
    labal_size_599 = labal_size_599 + new_label_num;
    mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

    if (labal_size_599 > max_labal_size_599)
    {
        throw reach_limit_error_string_MB; // after catching error, must call terminate_procedures_599(), otherwise this PLL cannot be reused
    }

    if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time_599)
            .count() > max_run_time_nanoseconds_599)
    {
        throw reach_limit_error_string_time; // after catching error, must call terminate_procedures_599(), otherwise this PLL cannot be reused
    }
}

void graph_hash_of_mixed_weighted_HL_HB_v1_thread_function_HBDIJ(int v_k, int N, int upper_k, bool use_hb)
{
	mtx_599[max_N_ID_for_mtx_599 - 1].lock();
	int used_id = Qid_599.front();
	Qid_599.pop();
	mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

	queue<int> P_changed_vertices, T_changed_vertices;
	vector<graph_hash_of_mixed_weighted_HL_PLL_v1_handle_t_for_sp> Q_handles(N);

	HBPLL_v1_node node;
	boost::heap::fibonacci_heap<HBPLL_v1_node> Q;
	two_hop_label_v1 xx;

	node.vertex = v_k;
	node.parent_vertex = v_k;
    node.hop = 0;
	node.priority_value = 0;
	Q_handles[v_k] = Q.push(node);
	P_dij_599[used_id][v_k] = 0;
	P_changed_vertices.push(v_k);

	mtx_599[v_k].lock();
	int L_v_k_size = L_temp_599[v_k].size();
	for (int i = 0; i < L_v_k_size; i++) {
		int L_v_k_i_vertex = L_temp_599[v_k][i].vertex;
		T_dij_599[used_id][L_v_k_i_vertex] = L_temp_599[v_k][i].distance;
		T_changed_vertices.push(L_v_k_i_vertex);
	}
	mtx_599[v_k].unlock();

	long long int new_label_num = 0;

	while (Q.size() > 0) {

		node = Q.top();
		Q.pop();
		int u = node.vertex;

		if (v_k <= u) {
			int u_parent = node.parent_vertex;
            int u_hop = node.hop;
			double P_u = node.priority_value;
			double P_u_with_error = P_u + 1e-5;
			double query_v_k_u = std::numeric_limits<double>::max();

    #ifdef _WIN32
			mtx_599[u].lock();
			auto L_u_size = L_temp_599[u].size(); // a vector<PLL_with_non_adj_reduction_sorted_label>
			mtx_599[u].unlock();
			for (int i = 0; i < L_u_size; i++) {
				mtx_599[u].lock();      // put lock in for loop is very slow, but it may be the only way under Windows
				double dis = L_temp_599[u][i].distance + T_bfs_599[used_id][L_temp_599[u][i].vertex];
				mtx_599[u].unlock();
				if (query_v_k_u > dis) { query_v_k_u = dis; }
			} //求query的值		
    #else
			mtx_599[u].lock();
			auto L_u_size1 = L_temp_599[u].size(); // a vector<PLL_with_non_adj_reduction_sorted_label>
			for (int i = 0; i < L_u_size1; i++) {
				if (L_temp_599[u][i].hop + T_dij_599[used_id][L_temp_599[u][i].vertex] <= u_hop || !use_hb)
                {
                    double dis = L_temp_599[u][i].distance + T_dij_599[used_id][L_temp_599[u][i].vertex];   // dont know why this code does not work under Windows
				    if (query_v_k_u > dis) { query_v_k_u = dis; }
                }      
			}
			mtx_599[u].unlock();
    #endif

			if (P_u_with_error < query_v_k_u) { // this is pruning
				xx.vertex = v_k;
				xx.distance = P_u;
				xx.parent_vertex = u_parent;
                xx.hop = u_hop;

				mtx_599[u].lock();
				L_temp_599[u].push_back(xx); //新增标签，并行时L_temp_599[u]里面的标签不一定是按照vertex ID排好序的，但是因为什么query时用了T_bfs_599的trick，没必要让L_temp_599[u]里面的标签排好序
				mtx_599[u].unlock();
				new_label_num++;

				/*下面是dij更新邻接点的过程，同时更新优先队列和距离*/
				int u_adj_size = ideal_graph_599[u].size();
				for (int i = 0; i < u_adj_size; i++) {
					int adj_v = ideal_graph_599[u][i].first; // this needs to be locked
					double ec = ideal_graph_599[u][i].second;
					if (P_dij_599[used_id][adj_v] == std::numeric_limits<double>::max()) { //尚未到达的点
						node.vertex = adj_v;
						node.parent_vertex = u;
						node.priority_value = P_u + ec;
                        node.hop = u_hop + 1;
                        /* R2添加的边: hop+1 */
                        if (new_edges_with_middle_v.find({u,adj_v}) != new_edges_with_middle_v.end() || new_edges_with_middle_v.find({adj_v, u}) != new_edges_with_middle_v.end())
                            node.hop ++;
                        Q_handles[adj_v] = Q.push(node);
						P_dij_599[used_id][adj_v] = node.priority_value;
						P_changed_vertices.push(adj_v);
					}
					else {
						if (P_dij_599[used_id][adj_v] > P_u + ec) {
							node.vertex = adj_v;
							node.parent_vertex = u;
							node.priority_value = P_u + ec;
                            node.hop = u_hop + 1;
							Q.update(Q_handles[adj_v], node);
							P_dij_599[used_id][adj_v] = node.priority_value;
						}
					}
				}
			}	
		}
	}

	while (P_changed_vertices.size() > 0) {
		P_dij_599[used_id][P_changed_vertices.front()] = std::numeric_limits<double>::max(); // reverse-allocate P values
		P_changed_vertices.pop();
	}
	while (T_changed_vertices.size() > 0) {
		T_dij_599[used_id][T_changed_vertices.front()] = std::numeric_limits<double>::max(); // reverse-allocate T values
		T_changed_vertices.pop();
	}

	mtx_599[v_k].lock();
	vector<two_hop_label_v1>(L_temp_599[v_k]).swap(L_temp_599[v_k]); // swap释放vector中多余空间： https://blog.csdn.net/qq_41929943/article/details/103190891 
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


void graph_hash_of_mixed_weighted_HB_v1_sort_labels_thread(vector<vector<two_hop_label_v1>> *output_L, int v_k)
{
    sort(L_temp_599[v_k].begin(), L_temp_599[v_k].end(), compare_two_hop_label_small_to_large);
    (*output_L)[v_k] = L_temp_599[v_k];
    vector<two_hop_label_v1>().swap(L_temp_599[v_k]); // clear new labels for RAM efficiency
}

vector<vector<two_hop_label_v1>> graph_hash_of_mixed_weighted_HB_v1_sort_labels(int N, int max_N_ID, int num_of_threads)
{
    vector<vector<two_hop_label_v1>> output_L(max_N_ID);
    vector<vector<two_hop_label_v1>> *p = &output_L;

    ThreadPool pool(num_of_threads);
    std::vector<std::future<int>> results;
    for (int v_k = 0; v_k < N; v_k++)
    {
        results.emplace_back(pool.enqueue([p, v_k] {
            graph_hash_of_mixed_weighted_HB_v1_sort_labels_thread(p, v_k);
            return 1;
        }));
    }
    for (auto &&result : results)
        result.get(); // all threads finish here

    return output_L;
}

/*the following parallel PLL_with_non_adj_reduction code cannot be run parallelly, due to the above globel values*/

void graph_hash_of_mixed_weighted_HB_v1(graph_v_of_v_idealID &input_graph, int max_N_ID, bool weighted, int num_of_threads, graph_hash_of_mixed_weighted_two_hop_case_info_v1 &case_info)
{
    //----------------------------------- step 1: initialization -----------------------------------
    cout << "step 1: initialization" << endl;

    auto begin = std::chrono::high_resolution_clock::now();
    /*information prepare*/
    begin_time_599 = std::chrono::high_resolution_clock::now();
    max_run_time_nanoseconds_599 = case_info.max_run_time_seconds * 1e9;
    labal_size_599 = 0;
    max_labal_size_599 = case_info.max_labal_size;

    full_two_hop_labels = !case_info.use_dummy_dij_search_in_PLL;

    if (max_N_ID > max_N_ID_for_mtx_599)
    {
        cout << "max_N_ID > max_N_ID_for_mtx_599; max_N_ID_for_mtx_599 is too small!" << endl;
        exit(1);
    }

    mtx_599[max_N_ID_for_mtx_599 - 1].lock();
    if (this_parallel_PLL_PSL_is_running_599 == true)
    {
        cout << "the following parallel PLL_with_non_adj_reduction code cannot be run parallelly, due to the above (static) globel values" << endl;
        exit(1);
    }
    this_parallel_PLL_PSL_is_running_599 = true;
    mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

    L_temp_599.resize(max_N_ID);
    int N = input_graph.size();

    /*sorting vertices*/
    // vector<vector<pair<int, double>>>().swap(adjs);
    // adjs.resize(max_N_ID);
    // vector<pair<int, double>>().swap(min_adjs);
    // min_adjs.resize(max_N_ID);
    // vector<pair<int, int>> sorted_vertices;
    // for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++)
    // {
    //     sorted_vertices.push_back({it->first, input_graph.degree(it->first)});
    //     adjs[it->first] = input_graph.adj_v_and_ec(it->first);
    //     min_adjs[it->first] = input_graph.min_adj(it->first);
    // }
    // sort(sorted_vertices.begin(), sorted_vertices.end(), compare_pair_second_large_to_small);
    // unordered_map<int, int> vertexID_old_to_new;
    // vertexID_new_to_old_599.resize(N);
    // for (int i = 0; i < N; i++)
    // {
    //     vertexID_old_to_new[sorted_vertices[i].first] = i;
    //     vertexID_new_to_old_599[i] = sorted_vertices[i].first;
    // }
    // vector<pair<int, int>>().swap(sorted_vertices);
    // ideal_graph_599 = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(input_graph, vertexID_old_to_new);
    
    ideal_graph_599 = input_graph;
    

    auto end = std::chrono::high_resolution_clock::now();
    case_info.time_initialization = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

    //---------------------------------------------------------------------------------------------------------------------------------------

    //----------------------------------------------- step 2: reduction ---------------------------------------------------------------
    cout << "step 2: reduction" << endl;

    /*redcution: add and remove certain edges*/
    case_info.reduction_measures_2019R2.clear(); // for using this function multiple times
    case_info.reduction_measures_2019R2.resize(max_N_ID);

    /*reduction 2; 用了2019 R2 enhance之后的图就是weighted，不能使用Unweighted bfs了！*/
    if (weighted == 0 && case_info.use_2019R2 + case_info.use_enhanced2019R2 > 0)
    {
        cout << "weighted = 1; // 用了2019 R2 enhance之后的图就是weighted，不能使用Unweighted bfs了！" << endl;
        weighted = 1;
    }
    begin = std::chrono::high_resolution_clock::now();
    if (case_info.use_2019R2)
    {
        case_info.MG_num = 0;
        for (int x = 0; x < N; x++)
        {
            if (ideal_graph_599[x].size() > 0)
            {
                if (x > ideal_graph_599[x][ideal_graph_599[x].size() - 1].first)
                {
                    case_info.reduction_measures_2019R2[x] = 2;
                    cout << "reduce " << x << endl;
                }
            }
        }
        for (int x = N - 1; x >= 0; x--)
        {
            if (case_info.reduction_measures_2019R2[x] == 2)
            {
                /*add edge*/
                auto it1 = ideal_graph_599[x].begin();
                for (int m = ideal_graph_599[x].size() - 1; m > 0; m--)
                {
                    for (int n = m - 1; n >= 0; n--)
                    {
                        double s_vPLUSv_t = (it1 + m)->second + (it1 + n)->second;
                        int e1 = (it1 + m)->first;
                        int e2 = (it1 + n)->first;
                        if (s_vPLUSv_t < graph_v_of_v_idealID_edge_weight(ideal_graph_599, e1, e2))
                        { // (s,v)+(v,t) is the shorter path or there is no edge between s and t
                            graph_v_of_v_idealID_add_edge(ideal_graph_599, e1, e2, s_vPLUSv_t);
                            new_edges_with_middle_v[{e1, e2}] = x;
                        }
                    }
                }
                /*remove edge*/
                case_info.MG_num++;
                graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_599, x);
            }
        }
    }
    else if (case_info.use_enhanced2019R2)
    { // e.g., 2019R2enhance_10
        case_info.MG_num = 0;
        for (int x = 0; x < N; x++)
        {
            if (ideal_graph_599[x].size() > 0)
            { // Prevent memory overflow
                if (x > ideal_graph_599[x][ideal_graph_599[x].size() - 1].first)
                { // Here, only one comparison is needed. A little trick.
                    case_info.reduction_measures_2019R2[x] = 2;
                    // cout << "reduce " << x << endl;
                }
            }
        }
        int bound = case_info.max_degree_MG_enhanced2019R2;
        for (int x = N - 1; x >= 0; x--)
        { // from low ranking to high ranking
            if (case_info.reduction_measures_2019R2[x] == 0 &&
                ideal_graph_599[x].size() <= bound)
            { // bound is the max degree for reduction
                bool no_adj_MG_vertices = true;
                for (auto it = ideal_graph_599[x].begin(); it != ideal_graph_599[x].end(); it++)
                {
                    if (case_info.reduction_measures_2019R2[it->first] == 2)
                    {
                        no_adj_MG_vertices = false;
                        break;
                    }
                }
                if (no_adj_MG_vertices)
                {
                    case_info.reduction_measures_2019R2[x] = 2; // new reduction
                    // cout << "new reduce " << x << endl;
                }
            }
        }
        for (int x = N - 1; x >= 0; x--)
        {
            if (case_info.reduction_measures_2019R2[x] == 2)
            {
                /*add edge*/
                auto it1 = ideal_graph_599[x].begin();
                for (int m = ideal_graph_599[x].size() - 1; m > 0; m--)
                {
                    for (int n = m - 1; n >= 0; n--)
                    {
                        double s_vPLUSv_t = (it1 + m)->second + (it1 + n)->second;
                        int e1 = (it1 + m)->first;
                        int e2 = (it1 + n)->first;
                        if (s_vPLUSv_t < graph_v_of_v_idealID_edge_weight(ideal_graph_599, e1, e2))
                        {
                            // (s,v)+(v,t) is the shorter path or there is no edge between s and t
                            graph_v_of_v_idealID_add_edge(ideal_graph_599, e1, e2, s_vPLUSv_t);
                            new_edges_with_middle_v[{e1, e2}] = x;
                        }
                    }
                }
                /*remove edge*/
                case_info.MG_num++;
                graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_599, x);
            }
        }
    }
    else if (case_info.use_non_adj_reduc_degree)
    {
        case_info.MG_num = 0;
        int bound = case_info.max_degree_MG_enhanced2019R2;
        for (int x = N - 1; x >= 0; x--)
        { // from low ranking to high ranking
            if (case_info.reduction_measures_2019R2[vertexID_new_to_old_599[x]] == 0 &&
                ideal_graph_599[x].size() <= bound)
            { // bound is the max degree for reduction
                bool no_adj_MG_vertices = true;
                for (auto it = ideal_graph_599[x].begin(); it != ideal_graph_599[x].end(); it++)
                {
                    if (case_info.reduction_measures_2019R2[vertexID_new_to_old_599[it->first]] == 2)
                    {
                        no_adj_MG_vertices = false;
                        break;
                    }
                }
                if (no_adj_MG_vertices)
                {
                    case_info.reduction_measures_2019R2[vertexID_new_to_old_599[x]] = 2; // new reduction
                    // cout << "new reduce " << vertexID_new_to_old_599[x] << endl;
                }
            }
        }
        for (int x = N - 1; x >= 0; x--)
        {
            if (case_info.reduction_measures_2019R2[vertexID_new_to_old_599[x]] == 2)
            {
                /*add edge*/
                auto it1 = ideal_graph_599[x].begin();
                for (int m = ideal_graph_599[x].size() - 1; m > 0; m--)
                {
                    for (int n = m - 1; n >= 0; n--)
                    {
                        double s_vPLUSv_t = (it1 + m)->second + (it1 + n)->second;
                        int e1 = (it1 + m)->first;
                        int e2 = (it1 + n)->first;
                        if (s_vPLUSv_t < graph_v_of_v_idealID_edge_weight(ideal_graph_599, e1, e2))
                        {
                            // (s,v)+(v,t) is the shorter path or there is no edge between s and t
                            graph_v_of_v_idealID_add_edge(ideal_graph_599, e1, e2, s_vPLUSv_t);
                            new_edges_with_middle_v[{e1, e2}] = x;
                        }
                    }
                }
                /*remove edge*/
                case_info.MG_num++;
                graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_599, x);
            }
        }
    }
    end = std::chrono::high_resolution_clock::now();
    case_info.time_2019R2_or_enhanced_pre = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

    //---------------------------------------------------------------------------------------------------------------------------------------

    //----------------------------------------------- step 3: generate labels ---------------------------------------------------------------
    cout << "step 3: generate labels" << endl;
    begin = std::chrono::high_resolution_clock::now();

    /*seaching shortest paths*/
    bool use_dij = case_info.use_dij;
    bool use_hb = case_info.use_hb;
    int upper_k = use_hb ? case_info.upper_k : 0;
    ThreadPool pool(num_of_threads);
    std::vector<std::future<int>> results;              
    int num_of_threads_per_push = num_of_threads * 100;
    if (!use_dij)
    {
        T_bfs_599.resize(num_of_threads);
        for (int i = 0; i < num_of_threads; i++)
        {
            T_bfs_599[i].resize(N);
            for (int j = 0; j < N; j++)
            {
                T_bfs_599[i][j].first = std::numeric_limits<double>::max();
                T_bfs_599[i][j].second = 0;
            }
            Qid_599.push(i);
        }
        int push_num = 0;
        for (int v_k = 0; v_k < N; v_k++)
        {
            if (ideal_graph_599[v_k].size() > 0)
            { // not from isolated vertices

                results.emplace_back(pool.enqueue([v_k, N, upper_k] { // pass const type value j to thread; [] can be empty
                    graph_hash_of_mixed_weighted_HL_HB_v1_thread_function_HBBFS(v_k, N, upper_k);
                    return 1; // return to results; the return type must be the same with results
                }));

                push_num++;
            }
            if (push_num % num_of_threads_per_push == 0)
            {
                for (auto &&result : results)
                    result.get(); // all threads finish here
                results.clear();
            }
        }
    }
    else
    {
        P_dij_599.resize(num_of_threads);
		T_dij_599.resize(num_of_threads);
		for (int i = 0; i < num_of_threads; i++)
		{
			P_dij_599[i].resize(N);
			T_dij_599[i].resize(N);
			for (int j = 0; j < N; j++)
			{
				P_dij_599[i][j] = std::numeric_limits<double>::max();
				T_dij_599[i][j] = std::numeric_limits<double>::max();
			}
			Qid_599.push(i);
		}
		int push_num = 0;
		for (int v_k = 0; v_k < N; v_k++) {
			if (ideal_graph_599[v_k].size() > 0) {  // not from isolated vertices
                results.emplace_back(
                    pool.enqueue([v_k, N, upper_k, use_hb] { // pass const type value j to thread; [] can be empty
                        graph_hash_of_mixed_weighted_HL_HB_v1_thread_function_HBDIJ(v_k, N, upper_k, use_hb);
                        return 1; // return to results; the return type must be the same with results
                        })
                );
				push_num++;
			}
			if (push_num % num_of_threads_per_push == 0) {
				for (auto&& result : results)
					result.get(); //all threads finish here
				results.clear();
			}
		}
    }

    for (auto &&result : results)
        result.get(); // all threads finish here

    end = std::chrono::high_resolution_clock::now();
    case_info.time_generate_labels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

    //---------------------------------------------------------------------------------------------------------------------------------------

    /*
        update predecessors for this non_adj_reduction,
        this update is for correct recursive direction.
    */
    begin = std::chrono::high_resolution_clock::now();
    for (auto it = new_edges_with_middle_v.begin() ; it != new_edges_with_middle_v.end() ; it++)
    {
        int e1 = it->first.first;
        int e2 = it->first.second;
        int middle_k = it->second;
        /*
            why just change the labels of e1 and e2 ?
            because 'parent_vertex' stands for 'next vertex on the shortest path', so it can only be shown in e1 and e2's labels
        */
        for (int j = L_temp_599[e1].size() - 1; j >= 0; j--)
        {
            if (L_temp_599[e1][j].parent_vertex == e2)
            {
                L_temp_599[e1][j].parent_vertex = middle_k;
            }
        }
        for (int j = L_temp_599[e2].size() - 1; j >= 0; j--)
        {
            if (L_temp_599[e2][j].parent_vertex == e1)
            {
                L_temp_599[e2][j].parent_vertex = middle_k;
            }
        }
    }
    end = std::chrono::high_resolution_clock::now();
    case_info.time_2019R2_or_enhanced_fixlabels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin) .count() / 1e9; // s

    //----------------------------------------------- step 4: canonical_repair---------------------------------------------------------------
    cout << "step 4: canonical_repair" << endl;

    /*canonical_repair based on the sorted new ID order, not the original ID order!*/
    if (case_info.use_canonical_repair)
    {
        begin = std::chrono::high_resolution_clock::now();
        // reduction_measures_2019R1_new_ID.resize(max_N_ID, 0);
        // reduction_measures_2019R2_new_ID.resize(max_N_ID, 0);
        f_2019R1_new_ID.resize(max_N_ID, 0);
        for (int i = 0; i < max_N_ID; i++)
        {
            sort(L_temp_599[i].begin(), L_temp_599[i].end(), compare_two_hop_label_small_to_large); // sort is necessary
        }
        end = std::chrono::high_resolution_clock::now();
        case_info.time_canonical_repair1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

        if (case_info.print_label_before_canonical_fix)
        {
            cout << "label_before_canonical_fix:" << endl;
            for (int i = 0; i < L_temp_599.size(); i++)
            {
                cout << "L[" << i << "]:\t";
                int ii = i;
                for (int j = 0; j < L_temp_599[ii].size(); j++)
                {
                    cout << "{" << L_temp_599[ii][j].vertex << "," << L_temp_599[ii][j].distance << "," << L_temp_599[ii][j].parent_vertex << "," << L_temp_599[ii][j].hop << "}";
                }
                cout << endl;
            }
        }

        begin = std::chrono::high_resolution_clock::now();
        canonical_repair_multi_threads(case_info.label_size_before_canonical_repair,
                                       case_info.label_size_after_canonical_repair,
                                       case_info.canonical_repair_remove_label_ratio, num_of_threads);
        end = std::chrono::high_resolution_clock::now();
        case_info.time_canonical_repair2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
    }

    //---------------------------------------------------------------------------------------------------------------------------------------

    //----------------------------------------------- step 5: update_old_IDs_in_label ---------------------------------------------------------------
    cout << "step 5: update_old_IDs_in_labels" << endl;
    begin = std::chrono::high_resolution_clock::now();
    
    /* sort lables */
    case_info.L = graph_hash_of_mixed_weighted_HB_v1_sort_labels(N, max_N_ID, num_of_threads);
    /* restore input_graph */
    ideal_graph_599 = input_graph;

    end = std::chrono::high_resolution_clock::now();
    case_info.time_update_old_IDs_in_labels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
    //---------------------------------------------------------------------------------------------------------------------------------------

    graph_hash_of_mixed_weighted_two_hop_clear_global_values();
}