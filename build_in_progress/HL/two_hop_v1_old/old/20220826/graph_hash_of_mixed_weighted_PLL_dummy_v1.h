#pragma once
#include <graph_hash_of_mixed_weighted/HL/two_hop_v1/graph_hash_of_mixed_weighted_PLL_v1.h>

/*
suppose that 0,...,max_non_dummy_ID are non-dummy vertices, max_non_dummy_ID+1,...,max_N_ID-1 are dummy vertices

PLL_dummy_v1:
end the search process after searching from all dummy vertices;
max dij priority value: 2*c(dummy edge) ---- stop search passing dummy vertex;
*/

double max_dij_priority_value = 2 * 1e6; // suppose that dummy edge has a weight of 1e6

void graph_hash_of_mixed_weighted_HL_PLL_v1_thread_function_dij_mixed_dummy(int v_k, int N)
{
	/*Pruned Dijkstra from vertex v_k; see Algorithm 1 in 2013 Japan SIGMOD paper*/

	mtx_595[max_N_ID_for_mtx_595 - 1].lock();
	int used_id = Qid_595.front();
	Qid_595.pop();
	mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

	queue<int> P_changed_vertices, T_changed_vertices;
	vector<graph_hash_of_mixed_weighted_HL_PLL_v1_handle_t_for_sp> Q_handles(N);

	PLL_v1_node_for_sp node;
	boost::heap::fibonacci_heap<PLL_v1_node_for_sp> Q;
	two_hop_label_v1 xx;

	node.vertex = v_k;
	node.parent_vertex = v_k;
	node.priority_value = 0;
	Q_handles[v_k] = Q.push(node);
	P_dij_595[used_id][v_k] = 0;
	P_changed_vertices.push(v_k);

	mtx_595[v_k].lock();
	int L_v_k_size = L_temp_595[v_k].size();
	for (int i = 0; i < L_v_k_size; i++) {
		int L_v_k_i_vertex = L_temp_595[v_k][i].vertex;
		T_dij_595[used_id][L_v_k_i_vertex] = L_temp_595[v_k][i].distance; //allocate T values for L_temp_595[v_k]
		T_changed_vertices.push(L_v_k_i_vertex);
	}
	mtx_595[v_k].unlock();
	//因为v-k的标签在从自己出发的过程中不会发生改变，并且在求query的过程中每次都会用到，所以可以提前取出来放在T数组，节省后面查找的时间

	long long int new_label_num = 0;

	while (Q.size() > 0) {

		node = Q.top();
		Q.pop();
		int u = node.vertex;

		if (v_k <= u) { // this pruning condition is not in the original 2013 PLL paper
			int u_parent = node.parent_vertex;
			float P_u = node.priority_value;
			float P_u_with_error = P_u + 1e-5;
			float query_v_k_u = std::numeric_limits<float>::max();

#ifdef _WIN32
			mtx_595[u].lock();
			auto L_u_size = L_temp_595[u].size(); // a vector<PLL_with_non_adj_reduction_sorted_label>
			mtx_595[u].unlock();
			for (int i = 0; i < L_u_size; i++) {
				mtx_595[u].lock();      // put lock in for loop is very slow, but it may be the only way under Windows
				float dis = L_temp_595[u][i].distance + T_dij_595[used_id][L_temp_595[u][i].vertex];
				mtx_595[u].unlock();
				if (query_v_k_u > dis) { query_v_k_u = dis; }
			} //求query的值		
#else
			mtx_595[u].lock();
			auto L_u_size1 = L_temp_595[u].size(); // a vector<PLL_with_non_adj_reduction_sorted_label>
			for (int i = 0; i < L_u_size1; i++) {
				float dis = L_temp_595[u][i].distance + T_dij_595[used_id][L_temp_595[u][i].vertex];   // dont know why this code does not work under Windows
				if (query_v_k_u > dis) { query_v_k_u = dis; }
			} //求query的值
			mtx_595[u].unlock();
#endif

			if (P_u_with_error < query_v_k_u) { // this is pruning
				xx.vertex = v_k;
				xx.distance = P_u;
				xx.parent_vertex = u_parent;

				mtx_595[u].lock();
				L_temp_595[u].push_back(xx); //新增标签，并行时L_temp_595[u]里面的标签不一定是按照vertex ID排好序的，但是因为什么query时用了T_dij_595的trick，没必要让L_temp_595[u]里面的标签排好序
				mtx_595[u].unlock();
				new_label_num++;

				/*下面是dij更新邻接点的过程，同时更新优先队列和距离*/
				int u_adj_size = ideal_graph_595[u].size();
				for (int i = 0; i < u_adj_size; i++) {
					int adj_v = ideal_graph_595[u][i].first; // this needs to be locked
					float ec = ideal_graph_595[u][i].second;
					if (P_dij_595[used_id][adj_v] == std::numeric_limits<float>::max()) { //尚未到达的点
						node.vertex = adj_v;
						node.parent_vertex = u;
						node.priority_value = P_u + ec;
						if (node.priority_value < max_dij_priority_value) {
							Q_handles[adj_v] = Q.push(node);
							P_dij_595[used_id][adj_v] = node.priority_value;
							P_changed_vertices.push(adj_v);
						}	
					}
					else {
						if (P_dij_595[used_id][adj_v] > P_u + ec) {
							node.vertex = adj_v;
							node.parent_vertex = u;
							node.priority_value = P_u + ec;
							Q.update(Q_handles[adj_v], node);
							P_dij_595[used_id][adj_v] = node.priority_value;
						}
					}
				}
			}
		}
	}

	while (P_changed_vertices.size() > 0) {
		P_dij_595[used_id][P_changed_vertices.front()] = std::numeric_limits<float>::max(); // reverse-allocate P values
		P_changed_vertices.pop();
	}
	while (T_changed_vertices.size() > 0) {
		T_dij_595[used_id][T_changed_vertices.front()] = std::numeric_limits<float>::max(); // reverse-allocate T values
		T_changed_vertices.pop();
	}

	mtx_595[v_k].lock();
	vector<two_hop_label_v1>(L_temp_595[v_k]).swap(L_temp_595[v_k]); // swap释放vector中多余空间： https://blog.csdn.net/qq_41929943/article/details/103190891 
	mtx_595[v_k].unlock();

	mtx_595[max_N_ID_for_mtx_595 - 1].lock();
	Qid_595.push(used_id);
	labal_size_595 = labal_size_595 + new_label_num;
	mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

	if (labal_size_595 > max_labal_size_595) {
		throw reach_limit_error_string;  // after catching error, must call terminate_procedures_595(), otherwise this PLL cannot be reused
	}

	if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time_595).count() > max_run_time_nanoseconds_595) {
		throw reach_limit_error_string;  // after catching error, must call terminate_procedures_595(), otherwise this PLL cannot be reused
	}
}













/*the following is wrong when dummy vertices may be in the tree index*/

//void PLL_dummy_v1(graph_hash_of_mixed_weighted& input_graph, int max_non_dummy_ID, int max_N_ID, bool weighted, int num_of_threads, graph_hash_of_mixed_weighted_two_hop_case_info_v1& case_info)
//{
//	//----------------------------------- step 1: initialization ------------------------------------------------------------------
//
//	auto begin = std::chrono::high_resolution_clock::now();
//	/*information prepare*/
//	begin_time_595 = std::chrono::high_resolution_clock::now();
//	max_run_time_nanoseconds_595 = case_info.max_run_time_seconds * 1e9;
//	labal_size_595 = 0;
//	max_labal_size_595 = case_info.max_labal_size;
//
//	if (max_N_ID > max_N_ID_for_mtx_595) {
//		cout << "max_N_ID > max_N_ID_for_mtx_595; max_N_ID_for_mtx_595 is too small!" << endl;
//		exit(1);
//	}
//
//	mtx_595[max_N_ID_for_mtx_595 - 1].lock();
//	if (this_parallel_PLL_PSL_is_running_595 == true) {
//		cout << "the following parallel PLL_with_non_adj_reduction code cannot be run parallelly, due to the above (static) globel values" << endl;
//		exit(1);
//	}
//	this_parallel_PLL_PSL_is_running_595 = true;
//	mtx_595[max_N_ID_for_mtx_595 - 1].unlock();
//
//	L_temp_595.resize(max_N_ID);
//	int N = input_graph.hash_of_vectors.size();
//
//
//	int not_searched_dummy_v_num = 0;
//	/*sorting vertices*/
//	vector <vector<pair<int, double>>>().swap(adjs);
//	adjs.resize(max_N_ID);
//	vector<pair<int, double>>().swap(min_adjs);
//	min_adjs.resize(max_N_ID);
//	vector<pair<int, int>> sorted_vertices;
//	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
//		int v_id = it->first;
//		int degree = input_graph.degree(it->first);
//		if (v_id > max_non_dummy_ID) {
//			not_searched_dummy_v_num++;
//		}
//		sorted_vertices.push_back({ v_id, degree });
//		adjs[it->first] = input_graph.adj_v_and_ec(it->first);
//		min_adjs[it->first] = input_graph.min_adj(it->first);
//	}
//	sort(sorted_vertices.begin(), sorted_vertices.end(), compare_pair_second_large_to_small);
//	unordered_map<int, int> vertexID_old_to_new;
//	vertexID_new_to_old_595.resize(N);
//	for (int i = 0; i < N; i++) {
//		vertexID_old_to_new[sorted_vertices[i].first] = i;
//		vertexID_new_to_old_595[i] = sorted_vertices[i].first;
//	}
//	vector<pair<int, int>>().swap(sorted_vertices);
//	ideal_graph_595 = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(input_graph, vertexID_old_to_new);
//
//
//	vector<pair<pair<int, int>, int>> new_edges_with_middle_v;	//Record newly added edges
//	/*redcution: add and remove certain edges*/
//	case_info.reduction_measures_2019R2.clear(); // for using this function multiple times
//	case_info.reduction_measures_2019R2.resize(max_N_ID);
//	/*clear graph_hash_of_mixed_weighted_HL_PLL_v1_f_2019R1*/
//	case_info.reduction_measures_2019R1.clear(); // for using this function multiple times
//	case_info.reduction_measures_2019R1.resize(max_N_ID);
//	case_info.f_2019R1.resize(max_N_ID);
//	std::iota(std::begin(case_info.f_2019R1), std::end(case_info.f_2019R1), 0); // Fill with 0, 1, ...
//
//	auto end = std::chrono::high_resolution_clock::now();
//	case_info.time_initialization = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
//
//	//---------------------------------------------------------------------------------------------------------------------------------------
//
//
//
//	//----------------------------------------------- step 2: reduction ---------------------------------------------------------------
//
//
//	/* redcution1: equivalence between vertices; we assume that vIDs in ideal_graph_595 are sorted by degree from large to small*/
//	if (case_info.use_2019R1)
//	{
//		auto begin = std::chrono::high_resolution_clock::now();
//		int ideal_graph_size = ideal_graph_595.size();
//
//		ThreadPool pool(num_of_threads);
//		std::vector< std::future<int> > results; // return typename: xxx
//		for (int i = 0; i < ideal_graph_size; i++)
//		{
//			vector<int>* xx = &(case_info.reduction_measures_2019R1);
//			vector<int>* yy = &(case_info.f_2019R1);
//			results.emplace_back(
//				pool.enqueue([i, ideal_graph_size, xx, yy] { // pass const type value j to thread; [] can be empty
//					update_2019R1_condition_PLL_with_non_adj_reduction(i, ideal_graph_size, xx, yy);
//					return 1; // return to results; the return type must be the same with results
//					})
//			);
//		}
//		for (auto&& result : results)
//			result.get(); // all threads finish here
//		/* remove edges */
//		case_info.reduce_V_num_2019R1 = 0;
//		for (int i = 0; i < max_N_ID; i++)
//		{
//			if (case_info.f_2019R1[i] != i)
//			{
//				case_info.reduce_V_num_2019R1++;
//				graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_595, vertexID_old_to_new[i]);
//			}
//		}
//		auto end = std::chrono::high_resolution_clock::now();
//		case_info.time_2019R1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
//	}
//
//	/*reduction 2; 用了2019 R2 enhance之后的图就是weighted，不能使用Unweighted bfs了！*/
//	if (weighted == 0 && case_info.use_2019R2 + case_info.use_enhanced2019R2 > 0) {
//		cout << "weighted = 1; // 用了2019 R2 enhance之后的图就是weighted，不能使用Unweighted bfs了！" << endl;
//		weighted = 1;
//	}
//	begin = std::chrono::high_resolution_clock::now();
//	if (case_info.use_2019R2) {
//		case_info.MG_num = 0;
//		for (int x = 0; x < N; x++) {
//			if (ideal_graph_595[x].size() > 0) {										//Prevent memory overflow
//				if (x > ideal_graph_595[x][ideal_graph_595[x].size() - 1].first) {		//Here, only one comparison is needed. A little trick.
//					case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] = 2;
//					//cout << "reduce " << vertexID_new_to_old_595[x] << endl;
//				}
//			}
//		}
//		for (int x = N - 1; x >= 0; x--) {
//			if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] == 2) {
//				/*add edge*/
//				auto it1 = ideal_graph_595[x].begin();
//				for (int m = ideal_graph_595[x].size() - 1; m > 0; m--)
//				{
//					for (int n = m - 1; n >= 0; n--)
//					{
//						double s_vPLUSv_t = (it1 + m)->second + (it1 + n)->second;
//						int e1 = (it1 + m)->first;
//						int e2 = (it1 + n)->first;
//						if (s_vPLUSv_t < graph_v_of_v_idealID_edge_weight(ideal_graph_595, e1, e2)) {
//							// (s,v)+(v,t) is the shorter path or there is no edge between s and t
//							graph_v_of_v_idealID_add_edge(ideal_graph_595, e1, e2, s_vPLUSv_t);
//							new_edges_with_middle_v.push_back({ {e1, e2}, x });
//						}
//					}
//				}
//				/*remove edge*/
//				case_info.MG_num++;
//				graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_595, x);
//			}
//		}
//	}
//	else if (case_info.use_enhanced2019R2) { // e.g., 2019R2enhance_10
//		case_info.MG_num = 0;
//		for (int x = 0; x < N; x++) {
//			if (ideal_graph_595[x].size() > 0) {										//Prevent memory overflow
//				if (x > ideal_graph_595[x][ideal_graph_595[x].size() - 1].first) {		//Here, only one comparison is needed. A little trick.
//					case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] = 2;
//					//cout << "reduce " << vertexID_new_to_old_595[x] << endl;
//				}
//			}
//		}
//		int bound = case_info.max_degree_MG_enhanced2019R2;
//		for (int x = N - 1; x >= 0; x--) { // from low ranking to high ranking
//			if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] == 0 && ideal_graph_595[x].size() <= bound) { // bound is the max degree for reduction
//				bool no_adj_MG_vertices = true;
//				for (auto it = ideal_graph_595[x].begin(); it != ideal_graph_595[x].end(); it++) {
//					if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[it->first]] == 2) {
//						no_adj_MG_vertices = false;
//						break;
//					}
//				}
//				if (no_adj_MG_vertices) {
//					case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] = 2; // new reduction
//					//cout << "new reduce " << vertexID_new_to_old_595[x] << endl;
//				}
//			}
//		}
//		for (int x = N - 1; x >= 0; x--) {
//			if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] == 2) {
//				/*add edge*/
//				auto it1 = ideal_graph_595[x].begin();
//				for (int m = ideal_graph_595[x].size() - 1; m > 0; m--)
//				{
//					for (int n = m - 1; n >= 0; n--)
//					{
//						double s_vPLUSv_t = (it1 + m)->second + (it1 + n)->second;
//						int e1 = (it1 + m)->first;
//						int e2 = (it1 + n)->first;
//						if (s_vPLUSv_t < graph_v_of_v_idealID_edge_weight(ideal_graph_595, e1, e2)) {
//							// (s,v)+(v,t) is the shorter path or there is no edge between s and t
//							graph_v_of_v_idealID_add_edge(ideal_graph_595, e1, e2, s_vPLUSv_t);
//							new_edges_with_middle_v.push_back({ {e1, e2}, x });
//						}
//					}
//				}
//				/*remove edge*/
//				case_info.MG_num++;
//				graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_595, x);
//			}
//		}
//	}
//	else if (case_info.use_non_adj_reduc_degree) {
//		case_info.MG_num = 0;
//		int bound = case_info.max_degree_MG_enhanced2019R2;
//		for (int x = N - 1; x >= 0; x--) { // from low ranking to high ranking
//			if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] == 0 && ideal_graph_595[x].size() <= bound) { // bound is the max degree for reduction
//				bool no_adj_MG_vertices = true;
//				for (auto it = ideal_graph_595[x].begin(); it != ideal_graph_595[x].end(); it++) {
//					if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[it->first]] == 2) {
//						no_adj_MG_vertices = false;
//						break;
//					}
//				}
//				if (no_adj_MG_vertices) {
//					case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] = 2; // new reduction
//					//cout << "new reduce " << vertexID_new_to_old_595[x] << endl;
//				}
//			}
//		}
//		for (int x = N - 1; x >= 0; x--) {
//			if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] == 2) {
//				/*add edge*/
//				auto it1 = ideal_graph_595[x].begin();
//				for (int m = ideal_graph_595[x].size() - 1; m > 0; m--)
//				{
//					for (int n = m - 1; n >= 0; n--)
//					{
//						double s_vPLUSv_t = (it1 + m)->second + (it1 + n)->second;
//						int e1 = (it1 + m)->first;
//						int e2 = (it1 + n)->first;
//						if (s_vPLUSv_t < graph_v_of_v_idealID_edge_weight(ideal_graph_595, e1, e2)) {
//							// (s,v)+(v,t) is the shorter path or there is no edge between s and t
//							graph_v_of_v_idealID_add_edge(ideal_graph_595, e1, e2, s_vPLUSv_t);
//							new_edges_with_middle_v.push_back({ {e1, e2}, x });
//						}
//					}
//				}
//				/*remove edge*/
//				case_info.MG_num++;
//				graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_595, x);
//			}
//		}
//	}
//	end = std::chrono::high_resolution_clock::now();
//	case_info.time_2019R2_or_enhanced_pre = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
//
//
//	//---------------------------------------------------------------------------------------------------------------------------------------
//
//
//	//----------------------------------------------- step 3: generate labels ---------------------------------------------------------------
//	begin = std::chrono::high_resolution_clock::now();
//
//	/*seaching shortest paths*/
//	ThreadPool pool(num_of_threads);
//	std::vector< std::future<int> > results; // return typename: xxx
//	int num_of_threads_per_push = num_of_threads * 100; // 每次push进去 num_of_threads_per_push 线程，如果没有异常，继续push进去num_of_threads_per_push线程；如果全都一起push进去必须全部线程都结束才能catch异常
//	if (weighted) {
//		P_dij_595.resize(num_of_threads);
//		T_dij_595.resize(num_of_threads);
//		for (int i = 0; i < num_of_threads; i++)
//		{
//			P_dij_595[i].resize(N);
//			T_dij_595[i].resize(N);
//			for (int j = 0; j < N; j++)
//			{
//				P_dij_595[i][j] = std::numeric_limits<float>::max();
//				T_dij_595[i][j] = std::numeric_limits<float>::max();
//			}
//			Qid_595.push(i);
//		}
//		int push_num = 0;
//		for (int v_k = 0; v_k < N; v_k++) {
//
//			if (vertexID_new_to_old_595[v_k] > max_non_dummy_ID) {
//				not_searched_dummy_v_num--;
//				if (not_searched_dummy_v_num == -1) { // all dummy vertices have been pushed into thread
//					for (auto&& result : results)
//						result.get(); //all threads finish here
//					results.clear();
//					break;
//				}
//			}
//
//			if (ideal_graph_595[v_k].size() > 0) {  // not from isolated vertices
//				results.emplace_back(
//					pool.enqueue([v_k, N] { // pass const type value j to thread; [] can be empty
//						graph_hash_of_mixed_weighted_HL_PLL_v1_thread_function_dij_mixed_dummy(v_k, N);
//						return 1; // return to results; the return type must be the same with results
//						})
//				);
//				push_num++;
//			}
//			if (push_num % num_of_threads_per_push == 0) {
//				for (auto&& result : results)
//					result.get(); //all threads finish here
//				results.clear();
//			}
//		}
//	}
//	else {
//		P_bfs_595.resize(num_of_threads);
//		T_bfs_595.resize(num_of_threads);
//		for (int i = 0; i < num_of_threads; i++)
//		{
//			P_bfs_595[i].resize(N);
//			T_bfs_595[i].resize(N);
//			for (int j = 0; j < N; j++)
//			{
//				P_bfs_595[i][j] = INT_MAX;
//				T_bfs_595[i][j] = INT_MAX;
//			}
//			Qid_595.push(i);
//
//		}
//		int push_num = 0;
//		for (int v_k = 0; v_k < N; v_k++) {
//
//			if (vertexID_new_to_old_595[v_k] > max_non_dummy_ID) {
//				not_searched_dummy_v_num--;
//				if (not_searched_dummy_v_num == -1) { // all dummy vertices have been pushed into thread
//					for (auto&& result : results)
//						result.get(); //all threads finish here
//					results.clear();
//					break;
//				}
//			}
//
//			if (ideal_graph_595[v_k].size() > 0) {  // not from isolated vertices
//				results.emplace_back(
//					pool.enqueue([v_k, N] { // pass const type value j to thread; [] can be empty
//						graph_hash_of_mixed_weighted_HL_PLL_v1_thread_function_bfs_mixed(v_k, N);
//						return 1; // return to results; the return type must be the same with results
//						})
//				);
//				push_num++;
//				if (push_num % num_of_threads_per_push == 0) {
//					for (auto&& result : results)
//						result.get(); // all threads finish here
//					results.clear();
//				}
//			}
//		}
//	}
//
//	for (auto&& result : results)
//		result.get(); //all threads finish here
//
//	end = std::chrono::high_resolution_clock::now();
//	case_info.time_generate_labels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
//	//---------------------------------------------------------------------------------------------------------------------------------------
//
//
//
//
///*
//	update predecessors for this non_adj_reduction,
//	this update is for correct recursive direction.
//*/
//	begin = std::chrono::high_resolution_clock::now();
//	for (int i = new_edges_with_middle_v.size() - 1; i >= 0; i--) {
//		int e1 = new_edges_with_middle_v[i].first.first;
//		int e2 = new_edges_with_middle_v[i].first.second;
//		int middle_k = new_edges_with_middle_v[i].second;
//		/*
//			why just change the labels of e1 and e2 ?
//			because 'parent_vertex' stands for 'next vertex on the shortest path', so it can only be shown in e1 and e2's labels
//		*/
//		for (int j = L_temp_595[e1].size() - 1; j >= 0; j--) {
//			if (L_temp_595[e1][j].parent_vertex == e2) {
//				L_temp_595[e1][j].parent_vertex = middle_k;
//			}
//		}
//		for (int j = L_temp_595[e2].size() - 1; j >= 0; j--) {
//			if (L_temp_595[e2][j].parent_vertex == e1) {
//				L_temp_595[e2][j].parent_vertex = middle_k;
//			}
//		}
//	}
//	end = std::chrono::high_resolution_clock::now();
//	case_info.time_2019R2_or_enhanced_fixlabels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
//
//
//
//
//
//
//
//	//----------------------------------------------- step 4: canonical_repair ---------------------------------------------------------------
//
//
//
//	//cout << "print_L_temp_595:" << endl;
//	//for (int i = 0; i < L_temp_595.size(); i++) {
//	//	cout << "L[" << i << "]=";
//	//	for (int j = 0; j < L_temp_595[i].size(); j++) {
//	//		cout << "{" << L_temp_595[i][j].vertex << "," << L_temp_595[i][j].distance << "," << L_temp_595[i][j].parent_vertex << "}";
//	//	}
//	//	cout << endl;
//	//}
//
//	/*canonical_repair based on the sorted new ID order, not the original ID order!*/
//	if (case_info.use_canonical_repair) {
//		begin = std::chrono::high_resolution_clock::now();
//		reduction_measures_2019R1_new_ID.resize(max_N_ID, 0);
//		reduction_measures_2019R2_new_ID.resize(max_N_ID, 0);
//		f_2019R1_new_ID.resize(max_N_ID, 0);
//		for (int i = 0; i < max_N_ID; i++) {
//			reduction_measures_2019R1_new_ID[vertexID_old_to_new[i]] = case_info.reduction_measures_2019R1[i];
//			reduction_measures_2019R2_new_ID[vertexID_old_to_new[i]] = case_info.reduction_measures_2019R2[i];
//			f_2019R1_new_ID[vertexID_old_to_new[i]] = vertexID_old_to_new[case_info.f_2019R1[i]];
//			sort(L_temp_595[i].begin(), L_temp_595[i].end(), compare_two_hop_label_small_to_large); // sort is necessary
//		}
//		graph_hash_of_mixed_weighted new_ID_g = graph_hash_of_mixed_weighted_update_vertexIDs(input_graph, vertexID_old_to_new);
// vector <vector<pair<int, double>>>().swap(adjs_new_IDs);
//adjs_new_IDs.resize(max_N_ID);
//vector<pair<int, double>>().swap(min_adjs_new_IDs);
//min_adjs_new_IDs.resize(max_N_ID);
//for (auto it = new_ID_g.hash_of_vectors.begin(); it != new_ID_g.hash_of_vectors.end(); it++) {
//	adjs_new_IDs[it->first] = new_ID_g.adj_v_and_ec(it->first);
//	min_adjs_new_IDs[it->first] = new_ID_g.min_adj(it->first);
//}
//		end = std::chrono::high_resolution_clock::now();
//		case_info.time_canonical_repair1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
//
//
//		begin = std::chrono::high_resolution_clock::now();
//		canonical_repair_multi_threads(new_ID_g, case_info.label_size_before_canonical_repair, case_info.label_size_after_canonical_repair, case_info.canonical_repair_remove_label_ratio, num_of_threads);
//		end = std::chrono::high_resolution_clock::now();
//		case_info.time_canonical_repair2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
//	}
//
//	//cout << "print_L_temp_595:" << endl;
//	//for (int i = 0; i < L_temp_595.size(); i++) {
//	//	cout << "L[" << i << "]=";
//	//	for (int j = 0; j < L_temp_595[i].size(); j++) {
//	//		cout << "{" << L_temp_595[i][j].vertex << "," << L_temp_595[i][j].distance << "," << L_temp_595[i][j].parent_vertex << "}";
//	//	}
//	//	cout << endl;
//	//}
//
//	//---------------------------------------------------------------------------------------------------------------------------------------
//
//
//	//----------------------------------------------- step 5: update_old_IDs_in_labels ---------------------------------------------------------------
//	begin = std::chrono::high_resolution_clock::now();
//
//	/*return L for old_IDs*/
//	case_info.L = graph_hash_of_mixed_weighted_HL_PLL_v1_transform_labels_to_old_vertex_IDs(N, max_N_ID, num_of_threads);
//
//	end = std::chrono::high_resolution_clock::now();
//	case_info.time_update_old_IDs_in_labels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
//	//---------------------------------------------------------------------------------------------------------------------------------------
//
//
//	graph_hash_of_mixed_weighted_two_hop_clear_global_values();
//}

