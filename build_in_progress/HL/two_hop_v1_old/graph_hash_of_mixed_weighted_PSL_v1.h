#pragma once
#include <iostream>
#include <tool_functions/ThreadPool.h>
#include <shared_mutex>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <graph_hash_of_mixed_weighted/two_graphs_operations/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID.h>
#include <build_in_progress/HL/two_hop_v1/graph_hash_of_mixed_weighted_two_hop_labels_v1.h>

 


void update_2019R1_condition_PSL_enhancedoriginalR2(int v1, int ideal_graph_size, vector<int>* reduction_measures_2, vector<int>* f_2019R1) {
	/*here, we assume v1 and v2 have the same number of adjs*/


	for (int v2 = v1 + 1; v2 < ideal_graph_size; v2++)
	{
		/*here is a little trick. it's certain that i has adjs no less than j*/
		if (ideal_graph_595[v1].size() > ideal_graph_595[v2].size())
			break; // no need to j++ any more


		int condition;

		if (graph_v_of_v_idealID_contain_edge(ideal_graph_595, v1, v2)) { // may be equivalent_2
			bool is_equivalent_2 = true;
			int size = ideal_graph_595[v1].size();
			auto it1 = ideal_graph_595[v1].begin();
			auto it2 = ideal_graph_595[v2].begin();
			while (it1 != ideal_graph_595[v1].end() && it2 != ideal_graph_595[v2].end()) {
				if (it1->first == v2) {
					it1++;
					continue;
				}
				if (it2->first == v1) {
					it2++;
					continue;
				}
				if (it1->first == it2->first) {
					if (abs(it1->second - it2->second) > 1e-5) {
						is_equivalent_2 = false;
						break;
					}
				}
				else {
					is_equivalent_2 = false;
					break;
				}
				it1++;
				it2++;
			}
			if (is_equivalent_2)
				condition = 2;
			else
				condition = 0;
		}
		else { // may be equivalent_1
			bool is_equivalent_1 = true;
			int size = ideal_graph_595[v1].size();
			auto it1 = ideal_graph_595[v1].begin();
			auto it2 = ideal_graph_595[v2].begin();
			while (it1 != ideal_graph_595[v1].end()) {
				if (it1->first == it2->first) {
					if (abs(it1->second - it2->second) > 1e-5) {
						is_equivalent_1 = false;
						break;
					}
				}
				else {
					is_equivalent_1 = false;
					break;
				}
				it1++;
				it2++;
			}
			if (is_equivalent_1)
				condition = 1;
			else
				condition = 0;
		}


		/*there are bugs below; 1. there should be locks below 2. even with locks, parallelly updating f function could cause wrong f mapping;
		but it seems that such bugs never occur, since v1 is increasing pushing into threads*/
		if (condition == 1)
		{
			(*reduction_measures_2)[vertexID_new_to_old_595[v1]] = 11;
			(*reduction_measures_2)[vertexID_new_to_old_595[v2]] = 11;
			(*f_2019R1)[vertexID_new_to_old_595[v2]] = (*f_2019R1)[vertexID_new_to_old_595[v1]];
		}
		else if (condition == 2)
		{
			(*reduction_measures_2)[vertexID_new_to_old_595[v1]] = 12;
			(*reduction_measures_2)[vertexID_new_to_old_595[v2]] = 12;
			(*f_2019R1)[vertexID_new_to_old_595[v2]] = (*f_2019R1)[vertexID_new_to_old_595[v1]];
		}


	}

}

void graph_hash_of_mixed_weighted_PSL_v1_thread_function_dij(int u, graph_hash_of_mixed_weighted_two_hop_case_info_v1* case_info)
{
	/*the id of T and dirty_flag used here */
	int used_id;
	mtx_595_1.lock();
	used_id = Qid_595.front(); // ID of thread
	Qid_595.pop();
	mtx_595_1.unlock();

	auto& reduction_measures_2019R2 = (*case_info).reduction_measures_2019R2;

	/* save labels of u in T */ //有了这个T之后,后面的query就会快一点,这个应该不用改
	int u_label_size = L_595[u].size();
	for (int i = 0; i < u_label_size; i++)
	{
		int w = L_595[u][i].vertex;
		double dis = L_595[u][i].distance;
		if (dirty_tag_595[used_id][w]) // the same hub may have redundancy, record the shortest distance
		{
			dirty_tag_595[used_id][w] = false; // first time w as hub
			T_dij_595[used_id][w] = dis;
		}
		else if (dis < T_dij_595[used_id][w]) // 
			T_dij_595[used_id][w] = dis; // T_dij_595 record the minimum dis for a hub
	}

	/*
	if (reduction_measures[vertexID_new_to_old[u]] == 2) //如果u是lms中的点,那么相当于它被删除掉了,不需要再给它做label了
	{
		for (int i = 0; i < u_label_size; i++)
			dirty_tag_595[used_id][L_595[u][i].vertex] = true; // recover dirty tag
	}
	*/

	int u_adj_size = ideal_graph_595[u].size(); // u的邻居的数量
	for (int i = 0; i < u_adj_size; i++)        // 遍历u的邻居v
	{
		int v = ideal_graph_595[u][i].first;
		double ec = ideal_graph_595[u][i].second;        //(u,v)距离
		int v_label_size = pos_595[v] + increment_595[v]; // increment_595[v] is label size in the last iteration
		if (reduction_measures_2019R2[vertexID_new_to_old_595[v]] == 2) //如果点v是lms集合中的点,那么进入这个循环,找到v除了u以外的邻接点(N2部分)
		{		
			int v_adj_size = ideal_graph_595[v].size();
			for (int j = 0; j < v_adj_size; j++)
			{
				if (ideal_graph_595[v][j].first == u) //需要排除是u点的情况
					continue;
				int v_1 = ideal_graph_595[v][j].first;     //这是点v1的index
				double ec_2 = ideal_graph_595[v][j].second; // v到v1的边的长度
				//遍历v_1的d-2级别label,用上pos_2
				for (int k = pos_2_595[v_1]; k < pos_595[v_1]; k++) 
				{
					int w = L_595[v_1][k].vertex; //d-2 hubs of v1, which is adj of v

					if (w > u)
						continue; // cannot be the hub of u

					/* query pruning */
					double dis = L_595[v_1][k].distance + ec + ec_2;
					int w_label_size = L_595[w].size();
					bool flag = false;
					for (int k = 0; k < w_label_size; k++)
						if (!dirty_tag_595[used_id][L_595[w][k].vertex]) // L_595[w][k].vertex is in the <d hub
						{
							double query_dis = T_dij_595[used_id][L_595[w][k].vertex] + L_595[w][k].distance;
							if (query_dis - 1e-6 <= dis)
							{
								flag = true;
								break;
							}
						}
					if (flag)
						continue;

					/*add a new label*/
					two_hop_label_v1 xx;
					// two_hop_label_v1 xx;
					xx.vertex = w;
					xx.distance = dis;
					xx.parent_vertex = v;
					L_temp_595[u].push_back(xx);
					if_continue_595 = true;
				}
			}
		}
		else // v不是lms集合中的点,那么进入这个循环,找到v的d-1label(N1部分)
		{
			for (int j = pos_595[v]; j < v_label_size; j++) // Ld-1(v); v_label_size is the position of the last element in Ld-1
			{
				int w = L_595[v][j].vertex; // visit each hub of v; d-1 hub of v
				if (w > u)
					continue; // cannot be the hub of u

				/* query pruning */
				double dis = L_595[v][j].distance + ec;
				int w_label_size = L_595[w].size();
				bool flag = false;
				for (int k = 0; k < w_label_size; k++)
					if (!dirty_tag_595[used_id][L_595[w][k].vertex]) // L_595[w][k].vertex is a common hub of u and w
					{
						double query_dis = T_dij_595[used_id][L_595[w][k].vertex] + L_595[w][k].distance; // T_dij_595[used_id][L_595[w][k].vertex] is u to common hub dis
						if (query_dis - 1e-6 <= dis)
						{
							flag = true;
							break;
						}
					}
				if (flag)
					continue;

				/*add a new label*/
				two_hop_label_v1 xx;
				// two_hop_label_v1 xx;
				xx.vertex = w;
				xx.distance = dis;
				xx.parent_vertex = v;
				L_temp_595[u].push_back(xx);
				if_continue_595 = true;

				//cout << "add (" << w << "," << dis << "," << v << ") to " << u << endl;
			}
		}
	}

	for (int i = 0; i < u_label_size; i++)
		dirty_tag_595[used_id][L_595[u][i].vertex] = true; // recover dirty tag



	mtx_595_1.lock();
	Qid_595.push(used_id);
	mtx_595_1.unlock();
	mtx_595_2.lock();
	labal_size_595 = labal_size_595 + L_temp_595[u].size(); // notably, there are redundent labels that may be removed in the end, so max_labal_size_595 is a limit for middle process, not only for the end
	mtx_595_2.unlock();

	if (labal_size_595 > max_labal_size_595) {
		throw reach_limit_error_string_MB;  // after catching error, must call terminate_procedures_595(), otherwise this PSL cannot be reused
	}
	if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time_595).count() > max_run_time_nanoseconds_595) {
		throw reach_limit_error_string_time;  // after catching error, must call terminate_procedures_595(), otherwise this PSL cannot be reused
	}
}

void graph_hash_of_mixed_weighted_PSL_v1_thread_function_add_new_labels(int u)
{
	pos_2_595[u] = pos_595[u];
	pos_595[u] += increment_595[u]; 
	increment_595[u] = L_temp_595[u].size();
	L_595[u].insert(L_595[u].end(), L_temp_595[u].begin(), L_temp_595[u].end());

	vector<two_hop_label_v1>().swap(L_temp_595[u]);
	vector<two_hop_label_v1>(L_595[u]).swap(L_595[u]);  // swap释放vector中多余空间： https://blog.csdn.net/qq_41929943/article/details/103190891 
}

void graph_hash_of_mixed_weighted_PSL_v1_transform_labels_to_old_vertex_IDs_element(vector<vector<two_hop_label_v1>>* output_L, int v_k) {

	int L_v_k_size = L_temp_595[v_k].size();
	for (int i = 0; i < L_v_k_size; i++) {
		auto it = &L_temp_595[v_k][i];
		it->vertex = vertexID_new_to_old_595[it->vertex];
		it->parent_vertex = vertexID_new_to_old_595[it->parent_vertex];
	}
	sort(L_temp_595[v_k].begin(), L_temp_595[v_k].end(), compare_two_hop_label_small_to_large);

	(*output_L)[vertexID_new_to_old_595[v_k]] = L_temp_595[v_k];
	vector<two_hop_label_v1>().swap(L_temp_595[v_k]); // clear new labels for RAM efficiency

}

vector<vector<two_hop_label_v1>> graph_hash_of_mixed_weighted_PSL_v1_transform_labels_to_old_vertex_IDs
(int N, int max_N, int num_of_threads) {

	/*time complexity: O(V*L*logL), where L is average number of labels per vertex*/

	vector<vector<two_hop_label_v1>> output_L(max_N);

	/*time complexity: O(V*L*logL), where L is average number of labels per vertex*/
	ThreadPool pool(num_of_threads);
	std::vector< std::future<int> > results; // return typename: xxx
	vector<vector<two_hop_label_v1>>* p = &output_L;
	for (int v_k = 0; v_k < N; v_k++) {
		results.emplace_back(
			pool.enqueue([p, v_k] { // pass const type value j to thread; [] can be empty
			graph_hash_of_mixed_weighted_PSL_v1_transform_labels_to_old_vertex_IDs_element(p, v_k);
			return 1; // return to results; the return type must be the same with results
		})
		);
	}
	for (auto&& result : results)
		result.get(); // all threads finish here

	return output_L;
}

void clean_L_element(int u) {

	/*the id of T and dirty_flag and pre_595 used here */
	int used_id;
	mtx_595_1.lock();
	used_id = Qid_595.front(); // ID of thread
	Qid_595.pop();
	mtx_595_1.unlock();

	/* save the smallest label for each hub of u in T */
	int u_label_size = L_595[u].size();
	for (int i = 0; i < u_label_size; i++)
	{
		int w = L_595[u][i].vertex; double dis = L_595[u][i].distance;
		if (dirty_tag_595[used_id][w]) // the same hub may have redundancy, record the shortest distance
		{
			dirty_tag_595[used_id][w] = false;
			T_dij_595[used_id][w] = dis;
			pre_595[used_id][w] = L_595[u][i].parent_vertex;
		}
		else if (dis < T_dij_595[used_id][w])
		{
			T_dij_595[used_id][w] = dis;
			pre_595[used_id][w] = L_595[u][i].parent_vertex;
		}
	}
	for (int i = 0; i < u_label_size; i++)
	{
		int w = L_595[u][i].vertex;
		if (!dirty_tag_595[used_id][w]) // push w only once
		{
			dirty_tag_595[used_id][w] = true; // push w only once
			two_hop_label_v1 xx;
			xx.vertex = w;
			xx.distance = T_dij_595[used_id][w];
			xx.parent_vertex = pre_595[used_id][w];
			L_temp_595[u].push_back(xx);
		}
	}
	vector<two_hop_label_v1>().swap(L_595[u]); // do not have two L in RAM simultaneously

	mtx_595_1.lock();
	Qid_595.push(used_id);
	mtx_595_1.unlock();

}

void graph_hash_of_mixed_weighted_PSL_v1
(graph_hash_of_mixed_weighted& input_graph, int max_N_ID, int num_of_threads, graph_hash_of_mixed_weighted_two_hop_case_info_v1& case_info) {

	//----------------------------------- step 1: initialization ------------------------------------------------------------------

	auto begin = std::chrono::high_resolution_clock::now();

	begin_time_595 = std::chrono::high_resolution_clock::now();
	max_run_time_nanoseconds_595 = case_info.max_run_time_seconds * 1e9;
	labal_size_595 = 0;
	max_labal_size_595 = case_info.max_labal_size;


	mtx_595_1.lock();
	if (this_parallel_PLL_PSL_is_running_595 == true) {
		cout << "the following parallel PSL_enhancedoriginalR2 code cannot be run parallelly, due to the above (static) globel values" << endl;
		exit(1);
	}
	this_parallel_PLL_PSL_is_running_595 = true;
	mtx_595_1.unlock();

	int N = input_graph.hash_of_vectors.size();

	/*sort vertices by degrees*/
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
	sort(sorted_vertices.begin(), sorted_vertices.end(), compare_pair_second_large_to_small);


	/*graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID*/
	unordered_map<int, int> vertexID_old_to_new;
	vertexID_new_to_old_595.resize(N);
	for (int i = 0; i < N; i++) {
		vertexID_old_to_new[sorted_vertices[i].first] = i;
		vertexID_new_to_old_595[i] = sorted_vertices[i].first;
	}
	vector<pair<int, int>>().swap(sorted_vertices);
	ideal_graph_595 = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(input_graph, vertexID_old_to_new);

	/*redcution: add and remove certain edges*/
	case_info.reduction_measures_2019R2.clear(); // for using this function multiple times
	case_info.reduction_measures_2019R2.resize(max_N_ID);
	/*clear case_info.f_2019R1*/
	case_info.reduction_measures_2019R1.clear(); // for using this function multiple times
	case_info.reduction_measures_2019R1.resize(max_N_ID);
	case_info.f_2019R1.resize(max_N_ID);
	std::iota(std::begin(case_info.f_2019R1), std::end(case_info.f_2019R1), 0); // Fill with 0, 1, ...


	auto end = std::chrono::high_resolution_clock::now();
	case_info.time_initialization = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	//---------------------------------------------------------------------------------------------------------------------------------------

	//graph_v_of_v_idealID_print(ideal_graph_595);



	//----------------------------------------------- step 2: reduction ---------------------------------------------------------------

	/* redcution1: equivalence between vertices */
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
				update_2019R1_condition_PSL_enhancedoriginalR2(i, ideal_graph_size, xx, yy);
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

	begin = std::chrono::high_resolution_clock::now();
	if (case_info.use_2019R2) { // no edge is removed here
		case_info.MG_num = 0;
		for (int x = 0; x < N; x++) {
			if (ideal_graph_595[x].size() > 0) {										//Prevent memory overflow
				if (x > ideal_graph_595[x][ideal_graph_595[x].size() - 1].first) {		//Here, only one comparison is needed. A little trick.
					case_info.MG_num++;
					case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] = 2;
					//cout << "reduce " << vertexID_new_to_old_595[x] << endl;
				}
			}
		}
	}
	else if (case_info.use_enhanced2019R2) { // e.g., 2019R2enhance_10
		case_info.MG_num = 0;
		for (int x = 0; x < N; x++) {
			if (ideal_graph_595[x].size() > 0) {										//Prevent memory overflow
				if (x > ideal_graph_595[x][ideal_graph_595[x].size() - 1].first) {		//Here, only one comparison is needed. A little trick.
					case_info.MG_num++;
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
					case_info.MG_num++;
					case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] = 2; // new reduction
					//cout << "new reduce " << vertexID_new_to_old_595[x] << endl;
				}
			}
		}
	}
	else if (case_info.use_non_adj_reduc_degree) { // e.g., 2019R2enhance_10
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
					case_info.MG_num++;
					case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] = 2; // new reduction
					//cout << "new reduce " << vertexID_new_to_old_595[x] << endl;
				}
			}
		}
	}
	end = std::chrono::high_resolution_clock::now();
	case_info.time_2019R2_or_enhanced_pre = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	//---------------------------------------------------------------------------------------------------------------------------------------





	//----------------------------------------------- step 3: generate labels ---------------------------------------------------------------
	begin = std::chrono::high_resolution_clock::now();

	L_595.resize(N); pos_595.resize(N); pos_2_595.resize(N); increment_595.resize(N);
	for (int v_k = 0; v_k < N; v_k++)  //initialization
	{
		if (ideal_graph_595[v_k].size() == 0) continue; // do not search isolated vertices
		pos_595[v_k] = 0; // record the position of the begining of the label in last iteration (d-1 iteration)
		pos_2_595[v_k] = 0; // for d-2 iteration
		increment_595[v_k] = 1; // record the size of the labels in last iteration (d-1 iteration)
		two_hop_label_v1 xx;
		xx.vertex = v_k;
		xx.distance = 0;
		xx.parent_vertex = v_k;
		L_595[v_k].push_back(xx);
	}

	T_dij_595.resize(num_of_threads);
	pre_595.resize(num_of_threads);
	dirty_tag_595.resize(num_of_threads);
	for (int i = 0; i < num_of_threads; i++)
	{
		pre_595[i].resize(N);
		T_dij_595[i].resize(N);
		dirty_tag_595[i].resize(N);
		for (int j = 0; j < N; j++) dirty_tag_595[i][j] = true;
		Qid_595.push(i);
	}

	L_temp_595.resize(N);
	if_continue_595 = true;
	int if_continue_595_false_time = 0;

	ThreadPool pool(num_of_threads);
	std::vector< std::future<int> > results;
	int num_of_threads_per_push = num_of_threads * 1000; // 每次push进去 num_of_threads_per_push 线程，如果没有异常，继续push进去num_of_threads_per_push线程；如果全都一起push进去必须全部线程都结束才能catch异常
	while (if_continue_595 || if_continue_595_false_time < 2) // since R2 skip some vertices, some new labels can only be generated when d increases 2, not 1, thus terminate the loop only when if_continue_595==false twice
	{
		//cout << "here" << endl;

		if (if_continue_595 == true) {
			if_continue_595_false_time = 0;
		}
		if_continue_595 = false;	

		int push_num = 0;
		for (int u = 0; u < N; u++)
		{
			if (ideal_graph_595[u].size() == 0 || case_info.reduction_measures_2019R2[vertexID_new_to_old_595[u]] == 2) continue; // do not search isolated vertices
			auto* xx = &case_info;
			results.emplace_back(
				pool.enqueue([u, xx] { // pass const type value j to thread; [] can be empty
				graph_hash_of_mixed_weighted_PSL_v1_thread_function_dij(u, xx); // new labels add into L_temp[u], but also read L in the process
				return 1; // return to results; the return type must be the same with results
			})
			);
			push_num++;
			if (push_num % num_of_threads_per_push == 0) {
				for (auto&& result : results)
					result.get(); //all threads finish here
				results.clear();
			}
		}

		for (auto&& result : results)
			result.get();
		results.clear();

		for (int u = 0; u < N; u++)
		{
			if (ideal_graph_595[u].size() == 0 || case_info.reduction_measures_2019R2[vertexID_new_to_old_595[u]] == 2) continue; // do not search isolated vertices
			results.emplace_back(
				pool.enqueue([u] { // pass const type value j to thread; [] can be empty
				graph_hash_of_mixed_weighted_PSL_v1_thread_function_add_new_labels(u); // new labels in L_temp[u] add into L[u], to avoid locking L in dij process
				return 1; // return to results; the return type must be the same with results
			})
			);
		}

		for (auto&& result : results)
			result.get();
		results.clear(); // 如果本轮没有异常则继续

		if (if_continue_595 == false) {
			if_continue_595_false_time++; // if if_continue_595_false_time==2, then if_continue_595==false twice
		}
	}

	/* clean the labels in L_595 to L_temp_595*/
	for (int u = 0; u < N; u++)
	{
		results.emplace_back(
			pool.enqueue([u] { // pass const type value j to thread; [] can be empty
				clean_L_element(u);
				return 1; // return to results; the return type must be the same with results
				})
		);
	}
	for (auto&& result : results)
		result.get();
	results.clear();

	end = std::chrono::high_resolution_clock::now();
	case_info.time_generate_labels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	//---------------------------------------------------------------------------------------------------------------------------------------




	//----------------------------------------------- step 4: canonical_repair ---------------------------------------------------------------

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
		for (int i = 0; i < N; i++) {
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
	begin = std::chrono::high_resolution_clock::now();

	/*return unordered_map_L for old_IDs*/
	case_info.L = graph_hash_of_mixed_weighted_PSL_v1_transform_labels_to_old_vertex_IDs(N, max_N_ID, num_of_threads);

	end = std::chrono::high_resolution_clock::now();
	case_info.time_update_old_IDs_in_labels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	//---------------------------------------------------------------------------------------------------------------------------------------



	graph_hash_of_mixed_weighted_two_hop_clear_global_values();
}
