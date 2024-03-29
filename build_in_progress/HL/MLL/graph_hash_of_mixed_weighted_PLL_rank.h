#pragma once

/*this PLL record ranks*/

#include <build_in_progress/HL/two_hop_v1/graph_hash_of_mixed_weighted_PLL_v1.h>

void graph_hash_of_mixed_weighted_PLL_rank(graph_hash_of_mixed_weighted& input_graph,
	int max_N_ID, bool weighted, int num_of_threads,
	graph_hash_of_mixed_weighted_two_hop_case_info_v1& case_info,
	vector<int>& rank, int lambda)
{
	//----------------------------------- step 1: initialization ------------------------------------------------------------------
	//cout << "step 1: initialization" << endl;

	auto begin = std::chrono::high_resolution_clock::now();
	/*information prepare*/
	begin_time_595 = std::chrono::high_resolution_clock::now();
	max_run_time_nanoseconds_595 = case_info.max_run_time_seconds * 1e9;
	labal_size_595 = 0;
	max_labal_size_595 = case_info.max_labal_size;

	full_two_hop_labels = !case_info.use_dummy_dij_search_in_PLL;

	if (max_N_ID > max_N_ID_for_mtx_595) {
		cout << "max_N_ID > max_N_ID_for_mtx_595; max_N_ID_for_mtx_595 is too small!" << endl;
		exit(1);
	}

	mtx_595[max_N_ID_for_mtx_595 - 1].lock();
	if (this_parallel_PLL_PSL_is_running_595 == true) {
		cout << "the following parallel PLL_with_non_adj_reduction code cannot be run parallelly, due to the above (static) globel values" << endl;
		exit(1);
	}
	this_parallel_PLL_PSL_is_running_595 = true;
	mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

	L_temp_595.resize(max_N_ID);
	int N = input_graph.hash_of_vectors.size();


	/*sorting vertices*/
	vector <vector<pair<int, double>>>().swap(adjs);
	adjs.resize(max_N_ID);
	vector<pair<int, double>>().swap(min_adjs);
	min_adjs.resize(max_N_ID);
	vector<pair<int, int>> sorted_vertices;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		sorted_vertices.push_back({ it->first, input_graph.degree(it->first) });
		adjs[it->first] = input_graph.adj_v_and_ec(it->first);
		min_adjs[it->first] = input_graph.min_adj(it->first);
	}
	// sort(sorted_vertices.begin(), sorted_vertices.end(), compare_pair_second_large_to_small);
	// unordered_map<int, int> vertexID_old_to_new;
	// vertexID_new_to_old_595.resize(N);
	// for (int i = 0; i < N; i++) {
	// 	vertexID_old_to_new[sorted_vertices[i].first] = i;
	// 	vertexID_new_to_old_595[i] = sorted_vertices[i].first;
	// }
	// vector<pair<int, int>>().swap(sorted_vertices);
	// ideal_graph_595 = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(input_graph, vertexID_old_to_new);

	unordered_map<int, int> vertexID_old_to_new;
	vertexID_new_to_old_595.resize(N);
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++)
	{

		int new_id = rank[it->first] - lambda - 1; // new_id = 0 is the vertex with lowest degree, then turn over
		new_id = N - 1 - new_id;
		// cout << it->first << " " << new_id << endl;
		vertexID_old_to_new[it->first] = new_id;
		vertexID_new_to_old_595[new_id] = it->first;
	}
	ideal_graph_595 = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(input_graph, vertexID_old_to_new);


	vector<pair<pair<int, int>, int>> new_edges_with_middle_v;	//Record newly added edges
	/*redcution: add and remove certain edges*/
	case_info.reduction_measures_2019R2.clear(); // for using this function multiple times
	case_info.reduction_measures_2019R2.resize(max_N_ID);
	/*clear graph_hash_of_mixed_weighted_HL_PLL_v1_f_2019R1*/
	case_info.reduction_measures_2019R1.clear(); // for using this function multiple times
	case_info.reduction_measures_2019R1.resize(max_N_ID);
	case_info.f_2019R1.resize(max_N_ID);
	std::iota(std::begin(case_info.f_2019R1), std::end(case_info.f_2019R1), 0); // Fill with 0, 1, ...

	auto end = std::chrono::high_resolution_clock::now();
	case_info.time_initialization = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	//---------------------------------------------------------------------------------------------------------------------------------------



	//----------------------------------------------- step 2: reduction ---------------------------------------------------------------
	//cout << "step 2: reduction" << endl;

	/* redcution1: equivalence between vertices; we assume that vIDs in ideal_graph_595 are sorted by degree from large to small*/
	if (case_info.use_2019R1)
	{
		auto begin = std::chrono::high_resolution_clock::now();
		int ideal_graph_size = ideal_graph_595.size();

		ThreadPool pool(num_of_threads);
		std::vector< std::future<int> > results; // return typename: xxx
		for (int i = 0; i < ideal_graph_size; i++)
		{
			vector<int>* xx = &(case_info.reduction_measures_2019R1);
			vector<int>* yy = &(case_info.f_2019R1);
			results.emplace_back(
				pool.enqueue([i, ideal_graph_size, xx, yy] { // pass const type value j to thread; [] can be empty
					update_2019R1_condition_PLL_with_non_adj_reduction(i, ideal_graph_size, xx, yy);
					return 1; // return to results; the return type must be the same with results
					})
			);
		}
		for (auto&& result : results)
			result.get(); // all threads finish here
		/* remove edges */
		case_info.reduce_V_num_2019R1 = 0;
		for (int i = 0; i < max_N_ID; i++)
		{
			if (case_info.f_2019R1[i] != i)
			{
				if (case_info.f_2019R1[case_info.f_2019R1[i]] != case_info.f_2019R1[i]) {
					cout << "f error due to the above parallelly updating f" << endl;
					//getchar();
				}
				case_info.reduce_V_num_2019R1++;
				graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_595, vertexID_old_to_new[i]);
			}
		}
		auto end = std::chrono::high_resolution_clock::now();
		case_info.time_2019R1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	}

	/*reduction 2; 用了2019 R2 enhance之后的图就是weighted，不能使用Unweighted bfs了！*/
	if (weighted == 0 && case_info.use_2019R2 + case_info.use_enhanced2019R2 > 0) {
		cout << "weighted = 1; // 用了2019 R2 enhance之后的图就是weighted，不能使用Unweighted bfs了！" << endl;
		weighted = 1;
	}
	begin = std::chrono::high_resolution_clock::now();
	if (case_info.use_2019R2) {
		case_info.MG_num = 0;
		for (int x = 0; x < N; x++) {
			if (ideal_graph_595[x].size() > 0) {										//Prevent memory overflow
				if (x > ideal_graph_595[x][ideal_graph_595[x].size() - 1].first) {		//Here, only one comparison is needed. A little trick.
					case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] = 2;
					//cout << "reduce " << vertexID_new_to_old_595[x] << endl;
				}
			}
		}
		for (int x = N - 1; x >= 0; x--) {
			if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] == 2) {
				/*add edge*/
				auto it1 = ideal_graph_595[x].begin();
				for (int m = ideal_graph_595[x].size() - 1; m > 0; m--)
				{
					for (int n = m - 1; n >= 0; n--)
					{
						double s_vPLUSv_t = (it1 + m)->second + (it1 + n)->second;
						int e1 = (it1 + m)->first;
						int e2 = (it1 + n)->first;
						if (s_vPLUSv_t < graph_v_of_v_idealID_edge_weight(ideal_graph_595, e1, e2)) {
							// (s,v)+(v,t) is the shorter path or there is no edge between s and t
							graph_v_of_v_idealID_add_edge(ideal_graph_595, e1, e2, s_vPLUSv_t);
							new_edges_with_middle_v.push_back({ {e1, e2}, x });
						}
					}
				}
				/*remove edge*/
				case_info.MG_num++;
				graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_595, x);
			}
		}
	}
	else if (case_info.use_enhanced2019R2) { // e.g., 2019R2enhance_10
		case_info.MG_num = 0;
		for (int x = 0; x < N; x++) {
			if (ideal_graph_595[x].size() > 0) {										//Prevent memory overflow
				if (x > ideal_graph_595[x][ideal_graph_595[x].size() - 1].first) {		//Here, only one comparison is needed. A little trick.
					case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] = 2;
					//cout << "reduce " << vertexID_new_to_old_595[x] << endl;
				}
			}
		}
		int bound = case_info.max_degree_MG_enhanced2019R2;
		for (int x = N - 1; x >= 0; x--) { // from low ranking to high ranking
			if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] == 0 && ideal_graph_595[x].size() <= bound) { // bound is the max degree for reduction
				bool no_adj_MG_vertices = true;
				for (auto it = ideal_graph_595[x].begin(); it != ideal_graph_595[x].end(); it++) {
					if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[it->first]] == 2) {
						no_adj_MG_vertices = false;
						break;
					}
				}
				if (no_adj_MG_vertices) {
					case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] = 2; // new reduction
					//cout << "new reduce " << vertexID_new_to_old_595[x] << endl;
				}
			}
		}
		for (int x = N - 1; x >= 0; x--) {
			if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] == 2) {
				/*add edge*/
				auto it1 = ideal_graph_595[x].begin();
				for (int m = ideal_graph_595[x].size() - 1; m > 0; m--)
				{
					for (int n = m - 1; n >= 0; n--)
					{
						double s_vPLUSv_t = (it1 + m)->second + (it1 + n)->second;
						int e1 = (it1 + m)->first;
						int e2 = (it1 + n)->first;
						if (s_vPLUSv_t < graph_v_of_v_idealID_edge_weight(ideal_graph_595, e1, e2)) {
							// (s,v)+(v,t) is the shorter path or there is no edge between s and t
							graph_v_of_v_idealID_add_edge(ideal_graph_595, e1, e2, s_vPLUSv_t);
							new_edges_with_middle_v.push_back({ {e1, e2}, x });
						}
					}
				}
				/*remove edge*/
				case_info.MG_num++;
				graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_595, x);
			}
		}
	}
	else if (case_info.use_non_adj_reduc_degree) {
		case_info.MG_num = 0;
		int bound = case_info.max_degree_MG_enhanced2019R2;
		for (int x = N - 1; x >= 0; x--) { // from low ranking to high ranking
			if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] == 0 && ideal_graph_595[x].size() <= bound) { // bound is the max degree for reduction
				bool no_adj_MG_vertices = true;
				for (auto it = ideal_graph_595[x].begin(); it != ideal_graph_595[x].end(); it++) {
					if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[it->first]] == 2) {
						no_adj_MG_vertices = false;
						break;
					}
				}
				if (no_adj_MG_vertices) {
					case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] = 2; // new reduction
					//cout << "new reduce " << vertexID_new_to_old_595[x] << endl;
				}
			}
		}
		for (int x = N - 1; x >= 0; x--) {
			if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] == 2) {
				/*add edge*/
				auto it1 = ideal_graph_595[x].begin();
				for (int m = ideal_graph_595[x].size() - 1; m > 0; m--)
				{
					for (int n = m - 1; n >= 0; n--)
					{
						double s_vPLUSv_t = (it1 + m)->second + (it1 + n)->second;
						int e1 = (it1 + m)->first;
						int e2 = (it1 + n)->first;
						if (s_vPLUSv_t < graph_v_of_v_idealID_edge_weight(ideal_graph_595, e1, e2)) {
							// (s,v)+(v,t) is the shorter path or there is no edge between s and t
							graph_v_of_v_idealID_add_edge(ideal_graph_595, e1, e2, s_vPLUSv_t);
							new_edges_with_middle_v.push_back({ {e1, e2}, x });
						}
					}
				}
				/*remove edge*/
				case_info.MG_num++;
				graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_595, x);
			}
		}
	}
	end = std::chrono::high_resolution_clock::now();
	case_info.time_2019R2_or_enhanced_pre = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s


	//---------------------------------------------------------------------------------------------------------------------------------------


	//----------------------------------------------- step 3: generate labels ---------------------------------------------------------------
	//cout << "step 3: generate labels" << endl;
	begin = std::chrono::high_resolution_clock::now();

	/*seaching shortest paths*/
	ThreadPool pool(num_of_threads);
	std::vector< std::future<int> > results; // return typename: xxx
	int num_of_threads_per_push = num_of_threads * 100; // 每次push进去 num_of_threads_per_push 线程，如果没有异常，继续push进去num_of_threads_per_push线程；如果全都一起push进去必须全部线程都结束才能catch异常
	if (weighted) {
		P_dij_595.resize(num_of_threads);
		T_dij_595.resize(num_of_threads);
		for (int i = 0; i < num_of_threads; i++)
		{
			P_dij_595[i].resize(N);
			T_dij_595[i].resize(N);
			for (int j = 0; j < N; j++)
			{
				P_dij_595[i][j] = std::numeric_limits<double>::max();
				T_dij_595[i][j] = std::numeric_limits<double>::max();
			}
			Qid_595.push(i);
		}
		int push_num = 0;
		for (int v_k = 0; v_k < N; v_k++) {
			if (ideal_graph_595[v_k].size() > 0) {  // not from isolated vertices
				if (case_info.use_dummy_dij_search_in_PLL) {
					results.emplace_back(
						pool.enqueue([v_k, N] { // pass const type value j to thread; [] can be empty
							graph_hash_of_mixed_weighted_HL_PLL_v1_thread_function_dij_mixed_dummy(v_k, N);
							return 1; // return to results; the return type must be the same with results
							})
					);
				}
				else {
					results.emplace_back(
						pool.enqueue([v_k, N] { // pass const type value j to thread; [] can be empty
							graph_hash_of_mixed_weighted_HL_PLL_v1_thread_function_dij_mixed(v_k, N);
							return 1; // return to results; the return type must be the same with results
							})
					);
				}
				push_num++;
			}
			if (push_num % num_of_threads_per_push == 0) {
				for (auto&& result : results)
					result.get(); //all threads finish here
				results.clear();
			}
		}
	}
	else {
		P_bfs_595.resize(num_of_threads);
		T_bfs_595.resize(num_of_threads);
		for (int i = 0; i < num_of_threads; i++)
		{
			P_bfs_595[i].resize(N);
			T_bfs_595[i].resize(N);
			for (int j = 0; j < N; j++)
			{
				P_bfs_595[i][j] = INT_MAX;
				T_bfs_595[i][j] = INT_MAX;
			}
			Qid_595.push(i);

		}
		int push_num = 0;
		for (int v_k = 0; v_k < N; v_k++) {
			if (ideal_graph_595[v_k].size() > 0) {  // not from isolated vertices
				results.emplace_back(
					pool.enqueue([v_k, N] { // pass const type value j to thread; [] can be empty
						graph_hash_of_mixed_weighted_HL_PLL_v1_thread_function_bfs_mixed(v_k, N);
						return 1; // return to results; the return type must be the same with results
						})
				);
				push_num++;
				if (push_num % num_of_threads_per_push == 0) {
					for (auto&& result : results)
						result.get(); // all threads finish here
					results.clear();
				}
			}
		}
	}

	for (auto&& result : results)
		result.get(); //all threads finish here

	end = std::chrono::high_resolution_clock::now();
	case_info.time_generate_labels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	//---------------------------------------------------------------------------------------------------------------------------------------




/*
	update predecessors for this non_adj_reduction,
	this update is for correct recursive direction.
*/
	begin = std::chrono::high_resolution_clock::now();
	for (int i = new_edges_with_middle_v.size() - 1; i >= 0; i--) {
		int e1 = new_edges_with_middle_v[i].first.first;
		int e2 = new_edges_with_middle_v[i].first.second;
		int middle_k = new_edges_with_middle_v[i].second;
		/*
			why just change the labels of e1 and e2 ?
			because 'parent_vertex' stands for 'next vertex on the shortest path', so it can only be shown in e1 and e2's labels
		*/
		for (int j = L_temp_595[e1].size() - 1; j >= 0; j--) {
			if (L_temp_595[e1][j].parent_vertex == e2) {
				L_temp_595[e1][j].parent_vertex = middle_k;
			}
		}
		for (int j = L_temp_595[e2].size() - 1; j >= 0; j--) {
			if (L_temp_595[e2][j].parent_vertex == e1) {
				L_temp_595[e2][j].parent_vertex = middle_k;
			}
		}
	}
	end = std::chrono::high_resolution_clock::now();
	case_info.time_2019R2_or_enhanced_fixlabels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s







	//----------------------------------------------- step 4: canonical_repair ---------------------------------------------------------------
	//cout << "step 4: canonical_repair" << endl;


	//cout << "print_L_temp_595:" << endl;
	//for (int i = 0; i < L_temp_595.size(); i++) {
	//	cout << "L[" << i << "]=";
	//	for (int j = 0; j < L_temp_595[i].size(); j++) {
	//		cout << "{" << L_temp_595[i][j].vertex << "," << L_temp_595[i][j].distance << "," << L_temp_595[i][j].parent_vertex << "}";
	//	}
	//	cout << endl;
	//}

	/*canonical_repair based on the sorted new ID order, not the original ID order!*/
	if (case_info.use_canonical_repair) {
		begin = std::chrono::high_resolution_clock::now();
		reduction_measures_2019R1_new_ID.resize(max_N_ID, 0);
		reduction_measures_2019R2_new_ID.resize(max_N_ID, 0);
		f_2019R1_new_ID.resize(max_N_ID, 0);
		for (int i = 0; i < max_N_ID; i++) {
			reduction_measures_2019R1_new_ID[vertexID_old_to_new[i]] = case_info.reduction_measures_2019R1[i];
			reduction_measures_2019R2_new_ID[vertexID_old_to_new[i]] = case_info.reduction_measures_2019R2[i];
			f_2019R1_new_ID[vertexID_old_to_new[i]] = vertexID_old_to_new[case_info.f_2019R1[i]];
			sort(L_temp_595[i].begin(), L_temp_595[i].end(), compare_two_hop_label_small_to_large); // sort is necessary
		}
		graph_hash_of_mixed_weighted new_ID_g = graph_hash_of_mixed_weighted_update_vertexIDs(input_graph, vertexID_old_to_new);
		vector <vector<pair<int, double>>>().swap(adjs_new_IDs);
		adjs_new_IDs.resize(max_N_ID);
		vector<pair<int, double>>().swap(min_adjs_new_IDs);
		min_adjs_new_IDs.resize(max_N_ID);
		for (auto it = new_ID_g.hash_of_vectors.begin(); it != new_ID_g.hash_of_vectors.end(); it++) {
			adjs_new_IDs[it->first] = new_ID_g.adj_v_and_ec(it->first);
			min_adjs_new_IDs[it->first] = new_ID_g.min_adj(it->first);
		}
		end = std::chrono::high_resolution_clock::now();
		case_info.time_canonical_repair1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s


		begin = std::chrono::high_resolution_clock::now();
		canonical_repair_multi_threads(new_ID_g, case_info.label_size_before_canonical_repair, case_info.label_size_after_canonical_repair, case_info.canonical_repair_remove_label_ratio, num_of_threads);
		end = std::chrono::high_resolution_clock::now();
		case_info.time_canonical_repair2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	}

	//cout << "print_L_temp_595:" << endl;
	//for (int i = 0; i < L_temp_595.size(); i++) {
	//	cout << "L[" << i << "]=";
	//	for (int j = 0; j < L_temp_595[i].size(); j++) {
	//		cout << "{" << L_temp_595[i][j].vertex << "," << L_temp_595[i][j].distance << "," << L_temp_595[i][j].parent_vertex << "}";
	//	}
	//	cout << endl;
	//}

	//---------------------------------------------------------------------------------------------------------------------------------------


	//----------------------------------------------- step 5: update_old_IDs_in_labels ---------------------------------------------------------------
	//cout << "step 5: update_old_IDs_in_labels" << endl;
	begin = std::chrono::high_resolution_clock::now();

	/*return L for old_IDs*/
	case_info.L = graph_hash_of_mixed_weighted_HL_PLL_v1_transform_labels_to_old_vertex_IDs(N, max_N_ID, num_of_threads);

	end = std::chrono::high_resolution_clock::now();
	case_info.time_update_old_IDs_in_labels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	//---------------------------------------------------------------------------------------------------------------------------------------


	graph_hash_of_mixed_weighted_two_hop_clear_global_values();
}