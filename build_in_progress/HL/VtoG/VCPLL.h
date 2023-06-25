#pragma once
#include <iostream>
#include <tool_functions/ThreadPool.h>
#include <shared_mutex>
#include <chrono>
#include <boost/heap/fibonacci_heap.hpp>
#include <set>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <graph_hash_of_mixed_weighted/two_graphs_operations/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID.h>
#include <build_in_progress/HL/VtoG/graph_hash_of_mixed_weighted_two_hop_labels_v1.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_update_vertexIDs.h>
#include <build_in_progress/HL/VtoG/graph_hash_of_mixed_weighted_PLL_v1.h>

/*global values*/
vector<vector<two_hop_label_v1>> Messages;
vector<unordered_map<int, pair<double, int>>> Hash;

void Scatter(int a) {
	int a_adj_size = ideal_graph_595[a].size();
	for (int i = 0; i < a_adj_size; i++) {
		int v = ideal_graph_595[a][i].first; // this needs to be locked
		double ec = ideal_graph_595[a][i].second;
	
		for (auto& it : L_595[a]) {
			if (it.vertex < v) {
				two_hop_label_v1 x;
				x.vertex = it.vertex;
				x.distance = it.distance + ec;
				x.parent_vertex = a;
				mtx_595[v].lock();
				Messages[v].push_back(x);
				mtx_595[v].unlock();
			}
		}	
	}
}

void Gather(int v, vector<int>& ActiveVertices, int used_id) {

	auto& TT = T_dij_595[used_id];
	
	vector<int> T_changed_vertices;

	/*new idea: remove redundant labels in Messages[v]   (there may be multiple labels in Messages[v] with the same hub, you only need to keep the minimum distance one)*/
	vector<two_hop_label_v1> new_Messages_v;
	for (auto& it : Messages[v]) {
		int u = it.vertex;
		if (TT[u] > it.distance) {
			TT[u] = it.distance;
		}
	}
	for (auto& it : Messages[v]) {
		int u = it.vertex;
		if (TT[u] == it.distance) {
			new_Messages_v.push_back(it);
			TT[u] = std::numeric_limits<double>::max();
		}
	}

	/*L_temp_595 is not needed when there is only 1 batch*/
	//int L_v_size = L_temp_595[v].size(); 
	//for (int i = 0; i < L_v_size; i++) {
	//	int L_v_i_vertex = L_temp_595[v][i].vertex;
	//	TT[L_v_i_vertex] = L_temp_595[v][i].distance; //allocate T values for L_temp_595[v_k]
	//	T_changed_vertices.push_back(L_v_i_vertex);
	//}
	mtx_595[v].lock();
	for (auto& label : Hash[v]) {
		TT[label.first] = label.second.first;
		T_changed_vertices.push_back(label.first);
	}
	mtx_595[v].unlock();

	vector<two_hop_label_v1>().swap(L_595[v]);

	for (auto& it : new_Messages_v) {

		if (TT[it.vertex] < it.distance + 1e-5) { // new idea, can accelerate
			continue;
		}

		int u = it.vertex;
		double query_v_u = std::numeric_limits<double>::max();
//#ifdef _WIN32
//		auto L_u_size = L_temp_595[u].size(); // a vector<PLL_with_non_adj_reduction_sorted_label>
//		for (int i = 0; i < L_u_size; i++) {
//			double dis = L_temp_595[u][i].distance + TT[L_temp_595[u][i].vertex];
//			if (query_v_u > dis) { query_v_u = dis; }
//		} //求query的值
//#else
//		auto L_u_size1 = L_temp_595[u].size(); // a vector<PLL_with_non_adj_reduction_sorted_label>
//		for (int i = 0; i < L_u_size1; i++) {
//			double dis = L_temp_595[u][i].distance + T_dij_595[used_id][L_temp_595[u][i].vertex];   // dont know why this code does not work under Windows
//			if (query_v_u > dis) { query_v_u = dis; }
//		} //求query的值
//#endif
		/*it is too slow to add the following codes, since it makes it necessary to lock Hash; but it's necessary when there is just 1 batch, as L_temp_595 is empty when processing this batch*/
		mtx_595[u].lock();
		for (auto& label : Hash[u]) {
			double dis = label.second.first + TT[label.first];
			if (query_v_u > dis) { query_v_u = dis; }
		}
		mtx_595[u].unlock();

		if (it.distance + 1e-5 < query_v_u) {
			L_595[v].push_back(it);
		}
	}
	if (L_595[v].size()) {
		mtx_595[v].lock();
		auto& it2 = Hash[v];
		for (auto& label : L_595[v]) { // since redundant labels in Messages[v] have been removed, it's guarantted that it2[label.vertex].first < label.distance
			it2[label.vertex] = { label.distance,label.parent_vertex };
		}
		mtx_595[v].unlock();
		mtx_595[max_N_ID_for_mtx_595 - 1].lock();
		ActiveVertices.push_back(v);
		mtx_595[max_N_ID_for_mtx_595 - 1].unlock();
	}

	vector<two_hop_label_v1>().swap(Messages[v]);

	for (auto v : T_changed_vertices) {
		TT[v] = std::numeric_limits<double>::max(); // reverse-allocate T values
	}
	vector<int>().swap(T_changed_vertices);
}

void batch_process(int N, vector<int>& batch_V, int num_of_threads, ThreadPool& pool, std::vector<std::future<int>>& results, bool use_hash_clean) {

	auto ActiveVertices = batch_V;
	for (auto u : ActiveVertices) {
		two_hop_label_v1 x;
		x.vertex = u;
		x.distance = 0;
		x.parent_vertex = u;
		L_595[u] = { x };
		L_temp_595[u].push_back(x);
		Hash[u][u] = { 0,u };
	}

	vector<vector<int>> Vs(num_of_threads);

	int mmm = 0;

	while (ActiveVertices.size()) {
		cout << "while Scatter " << ++mmm << endl;

		int xx1 = 0;
		for (auto a : ActiveVertices) {
			Vs[xx1++ % num_of_threads].push_back(a);
		}
		for (int i = 0; i < num_of_threads; i++) {
			auto& p = Vs[i];
			results.emplace_back(
				pool.enqueue([&p] { // pass const type value j to thread; [] can be empty
					for (auto v : p) {
						Scatter(v);
					}
					vector<int>().swap(p);
					return 1; // return to results; the return type must be the same with results
					})
			);
		}
		for (auto&& result : results)
			result.get(); //all threads finish here
		results.clear();

		cout << "while Gather " << mmm << endl;

		vector<int>().swap(ActiveVertices);
		auto& xx = ActiveVertices;

		for (int v = 0; v < N; v++) {
			if (Messages[v].size()) {
				Vs[xx1++ % num_of_threads].push_back(v);
			}
		}
		for (int i = 0; i < num_of_threads; i++) {
			auto& p = Vs[i];
			results.emplace_back(
				pool.enqueue([&p, &xx, i] { // pass const type value j to thread; [] can be empty
					for (auto v : p) {
						Gather(v, xx, i);
					}
					vector<int>().swap(p);
					return 1; // return to results; the return type must be the same with results
					})
			);
		}
		for (auto&& result : results)
			result.get(); //all threads finish here
		results.clear();
	}

	/*hash cleaning and insert L_temp_595*/
	for (int v = 0; v < N; v++) {
		if (use_hash_clean) {
			results.emplace_back(
				pool.enqueue([v] { // pass const type value j to thread; [] can be empty					
					for (auto& label : Hash[v]) {
						int u = label.first;
						/*query d_uv*/
						double query_v_u = std::numeric_limits<double>::max();
						for (auto& label2 : Hash[u]) {
							if (Hash[v].count(label2.first) == 1) {
								double query_d = Hash[v][label2.first].first + label2.second.first;
								if (query_d < query_v_u) {
									query_v_u = query_d;
								}
							}
						}
						/*clean*/
						if (query_v_u + 1e-5 > label.second.first) {
							two_hop_label_v1 x;
							x.vertex = label.first;
							x.distance = label.second.first;
							x.parent_vertex = label.second.second;
							L_temp_595[v].push_back(x);
						}
						//else {
						//	cout << "clean" << endl;
						//}
					}
					return 1; // return to results; the return type must be the same with results
					})
			);
		}
		else {
			results.emplace_back(
				pool.enqueue([v] { // pass const type value j to thread; [] can be empty					
					for (auto& label : Hash[v]) {
						two_hop_label_v1 x;
						x.vertex = label.first;
						x.distance = label.second.first;
						x.parent_vertex = label.second.second;
						L_temp_595[v].push_back(x);
					}
					return 1; // return to results; the return type must be the same with results
					})
			);
		}
		vector<two_hop_label_v1>().swap(L_595[v]);
	}
	for (auto&& result : results)
		result.get(); //all threads finish here
	results.clear();

	for (int v = 0; v < N; v++) {
		unordered_map<int, pair<double, int>>().swap(Hash[v]);
	}
}

void clean_labels(int N, ThreadPool& pool, std::vector<std::future<int>>& results) {

	for (int v = 0; v < N; v++) {
		results.emplace_back(
			pool.enqueue([v] { // pass const type value j to thread; [] can be empty

				mtx_595[max_N_ID_for_mtx_595 - 1].lock();
				int used_id = Qid_595.front();
				Qid_595.pop();
				mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

				queue<int> T_changed_vertices;

				for (auto& it : L_temp_595[v]) {
					if (T_dij_595[used_id][it.vertex] > it.distance) {
						T_dij_595[used_id][it.vertex] = it.distance;
					}
				}
				vector<two_hop_label_v1> new_Lv;
				for (auto& it : L_temp_595[v]) {
					if (T_dij_595[used_id][it.vertex] != std::numeric_limits<double>::max() && T_dij_595[used_id][it.vertex] == it.distance) {
						new_Lv.push_back(it);
						T_dij_595[used_id][it.vertex] = std::numeric_limits<double>::max();
					}
				}

				L_temp_595[v] = new_Lv;

				mtx_595[max_N_ID_for_mtx_595 - 1].lock();
				Qid_595.push(used_id);
				mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

				return 1; // return to results; the return type must be the same with results
				})
		);
	}
	for (auto&& result : results)
		result.get(); //all threads finish here
	results.clear();
}

void VCPLL(graph_hash_of_mixed_weighted& input_graph, int max_N_ID, bool use_hash_clean, int num_of_threads, graph_hash_of_mixed_weighted_two_hop_case_info_v1& case_info)
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

	L_595.resize(max_N_ID); // deltaL
	L_temp_595.resize(max_N_ID); // L
	Messages.resize(max_N_ID);
	Hash.resize(max_N_ID);
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
	sort(sorted_vertices.begin(), sorted_vertices.end(), compare_pair_second_large_to_small);
	unordered_map<int, int> vertexID_old_to_new;
	vertexID_new_to_old_595.resize(N);
	for (int i = 0; i < N; i++) {
		vertexID_old_to_new[sorted_vertices[i].first] = i;
		vertexID_new_to_old_595[i] = sorted_vertices[i].first;
		//vertexID_old_to_new[i] = i;
		//vertexID_new_to_old_595[i] = i;
	}
	//cout << "rank above" << endl;
	vector<pair<int, int>>().swap(sorted_vertices);
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

	//----------------------------------------------- step 3: generate labels ---------------------------------------------------------------
	//cout << "step 3: generate labels" << endl;
	begin = std::chrono::high_resolution_clock::now();

	/*seaching shortest paths*/
	ThreadPool pool(num_of_threads);
	std::vector< std::future<int> > results; // return typename: xxx
	P_dij_595.resize(num_of_threads);
	T_dij_595.resize(num_of_threads);
	for (int i = 0; i < num_of_threads; i++)
	{
		P_dij_595[i].resize(N, std::numeric_limits<double>::max());
		T_dij_595[i].resize(N, std::numeric_limits<double>::max());
		Qid_595.push(i);
	}
	int batch_size = N; // 512 is very slow
	int batch_ID = 0;
	vector<int> batch_V;
	for (int v_k = 0; v_k < N; v_k++) {
		if (ideal_graph_595[v_k].size() > 0) {  // not from isolated vertices
			batch_V.push_back(v_k);
		}
		if (batch_V.size() == batch_size || v_k == N - 1) {
			cout << "batch_process " << ++batch_ID << endl;
			batch_process(N, batch_V, num_of_threads, pool, results, use_hash_clean);
			vector<int>().swap(batch_V);
		}
	}
	end = std::chrono::high_resolution_clock::now();
	case_info.time_generate_labels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	//---------------------------------------------------------------------------------------------------------------------------------------

	clean_labels(N, pool, results);/*remove redundant labels*/


	//graph_v_of_v_idealID_print(ideal_graph_595);

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

	//----------------------------------------------- step 5: update_old_IDs_in_labels ---------------------------------------------------------------
	//cout << "step 5: update_old_IDs_in_labels" << endl;
	begin = std::chrono::high_resolution_clock::now();

	/*return L for old_IDs*/
	case_info.L = graph_hash_of_mixed_weighted_HL_PLL_v1_transform_labels_to_old_vertex_IDs(N, max_N_ID, num_of_threads);

	end = std::chrono::high_resolution_clock::now();
	case_info.time_update_old_IDs_in_labels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	//---------------------------------------------------------------------------------------------------------------------------------------


	graph_hash_of_mixed_weighted_two_hop_clear_global_values();
	vector<vector<two_hop_label_v1>>().swap(Messages);
	vector<unordered_map<int, pair<double, int>>>().swap(Hash);
}