#pragma once
#include <iostream>
#include <chrono>
#include <boost/heap/fibonacci_heap.hpp>
#include <build_in_progress/HL/dgraph/dgraph_two_hop_label.h>

void propagate(dgraph_v_of_v<two_hop_weight_type>* input_graph, int k, int u) {

	if (PSL_throw_error) {
		return;
	}
	if (labal_size_PSL > max_labal_size_PSL) { // for PSL, this is not to limit the final label size, but to limit the middle process label size
		PSL_throw_error = true;
		throw reach_limit_error_string_MB;  // after catching error, must call dgraph_clear_global_values_PLL_PSL, otherwise this PLL cannot be reused
	}
	if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time_PSL).count() > max_run_time_nanoseconds_PSL) {
		PSL_throw_error = true;
		throw reach_limit_error_string_time;  // after catching error, must call dgraph_clear_global_values_PLL_PSL, otherwise this PLL cannot be reused
	}

	mtx.lock();
	int current_tid = Qid_595.front();
	Qid_595.pop();
	mtx.unlock();

	int N = input_graph->INs.size();
	vector<int> dist_dirt_change;
	vector<vector<two_hop_label>>* Lit; vector<vector<two_hop_label>>* Lit2;
	std::vector<std::pair<int, two_hop_weight_type>>* adj_list;
	if (k == 1) {
		Lit = &L_temp_in;
		Lit2 = &L_temp_out;
		adj_list = &(input_graph->INs)[u];
	}
	else if (k == 0) {
		Lit = &L_temp_out;
		Lit2 = &L_temp_in;
		adj_list = &(input_graph->OUTs)[u];
	}

	for (auto& cur : (*Lit)[u]) {
		int v = cur.vertex;
		two_hop_weight_type d = cur.distance;
		dist[current_tid][v] = std::min(dist[current_tid][v], d);
		if (dirt[current_tid][v]) { // using dirt is slightly faster than not using dirt
			dirt[current_tid][v] = false;
			dist[current_tid][v] = d;
			dist_dirt_change.push_back(v);
		}
		else {
			dist[current_tid][v] = std::min(dist[current_tid][v], d);
		}
	}

	for (auto label : *adj_list) {
		int v = label.first;
		two_hop_weight_type dd = label.second;
		for (int i = pos_595[k][v]; i < pos_595[k][v] + increment[k][v]; i++) {
			int x = (*Lit)[v][i].vertex;
			if (x <= u) { // lower ID has higher rank
				two_hop_weight_type d_temp = dd + (*Lit)[v][i].distance;
				bool flag = false;
				for (auto it : (*Lit2)[x]) {
					int y = it.vertex;
					if (dist[current_tid][y] + it.distance <= d_temp) {
						flag = true;
						break;
					}
				}
				if (flag) {
					continue;
				}
				L_PSL_temp[k][u].push_back({ x, d_temp });
			}
		}
	}

	for (int v : dist_dirt_change) {
		dist[current_tid][v] = std::numeric_limits<two_hop_weight_type>::max();
		dirt[current_tid][v] = true;
	}

	mtx.lock();
	Qid_595.push(current_tid);
	labal_size_PSL += L_PSL_temp[k][u].size();
	mtx.unlock();
}

void append(int k, int u) {
	pos_2_595[k][u] = pos_595[k][u];
	pos_595[k][u] += increment[k][u];
	increment[k][u] = L_PSL_temp[k][u].size();

	if (k == 1) {
		L_temp_in[u].insert(L_temp_in[u].end(), L_PSL_temp[k][u].begin(), L_PSL_temp[k][u].end());
	}
	else {
		L_temp_out[u].insert(L_temp_out[u].end(), L_PSL_temp[k][u].begin(), L_PSL_temp[k][u].end());
	}

	if (increment[k][u]) {
		no_new_label_PSL = false;
	}
	vector<two_hop_label>().swap(L_PSL_temp[k][u]);
}

void clean_incorrect_labels(int N, ThreadPool& pool, std::vector<std::future<int>>& results) {

	for (int u = 0; u < N; u++) {
		results.emplace_back(pool.enqueue([u, N] {

			mtx.lock();
			int current_tid = Qid_595.front();
			Qid_595.pop();
			mtx.unlock();
			vector<two_hop_label>* Lit;
			for (int k = 0; k < 2; k++) {
				if (k == 1) {
					Lit = &L_temp_in[u];
				}
				else if (k == 0) {
					Lit = &L_temp_out[u];
				}
				/* save the smallest label for each hub of u in T */
				vector<two_hop_label> L_temp;
				vector<int> dirt_change;
				int u_label_size = Lit->size();
				for (int i = 0; i < u_label_size; i++) {
					int w = (*Lit)[i].vertex;
					two_hop_weight_type dis = (*Lit)[i].distance;
					if (dirt[current_tid][w]) { // the same hub may have redundancy, record the shortest distance            
						dirt[current_tid][w] = false;
						dirt_change.push_back(w);
						dist[current_tid][w] = dis;
					}
					else if (dis < dist[current_tid][w]) {
						dist[current_tid][w] = dis;
					}
				}
				sort(dirt_change.begin(), dirt_change.end());
				for (int i : dirt_change) {
					dirt[current_tid][i] = true;
					two_hop_label xx;
					xx.vertex = i;
					xx.distance = dist[current_tid][i];
					L_temp.push_back(xx);
				}
				*Lit = L_temp;
			}
			mtx.lock();
			Qid_595.push(current_tid);
			mtx.unlock();

			return 1; }));
	}

	for (auto&& result : results) {
		result.get();
	}
	results.clear();
}

void dgraph_PSL(dgraph_v_of_v<two_hop_weight_type>& input_graph, int num_of_threads, dgraph_case_info_v1& case_info) {

	mtx_595[max_N_ID_for_mtx_595 - 1].lock();
	if (this_parallel_PSL_is_running == true) {
		cout << "dgraph_PSL cannot be run parallelly, due to the above (static) globel values" << endl;
		exit(1);
	}
	this_parallel_PSL_is_running = true;
	mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

	begin_time_PSL = std::chrono::high_resolution_clock::now();
	max_run_time_nanoseconds_PSL = case_info.max_run_time_seconds * 1e9;
	labal_size_PSL = 0;
	max_labal_size_PSL = case_info.max_labal_bit_size / sizeof(two_hop_label);

	int N = input_graph.INs.size();
	L_temp_in.resize(N);
	L_temp_out.resize(N);

	// prepare for thread pool
	dist.resize(num_of_threads);
	dist2.resize(num_of_threads);
	dirt.resize(num_of_threads);
	for (int i = 0; i < num_of_threads; i++) {
		dist[i].resize(N, std::numeric_limits<two_hop_weight_type>::max());
		dist2[i].resize(N, std::numeric_limits<two_hop_weight_type>::max());
		dirt[i].resize(N, true);
		Qid_595.push(i);
	}

	for (int k = 0; k < 2; k++) {
		L_PSL_temp[k].resize(N);
		pos_595[k].resize(N, 0);
		pos_2_595[k].resize(N, 0);
		increment[k].resize(N, 1);
	}

	case_info.time1_PLL_PSL_initialization = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time_PSL).count() / 1e9;
	auto begin1 = std::chrono::high_resolution_clock::now();
	//---------------------------------------------------------------------------------

	/*label generation*/
	ThreadPool pool(num_of_threads);
	std::vector<std::future<int>> results;
	for (int i = 0; i < N; i++) {
		L_temp_in[i].push_back({ i, 0 });
		L_temp_out[i].push_back({ i, 0 });
	}
	dgraph_v_of_v<two_hop_weight_type>* g_it = &input_graph;
	no_new_label_PSL = false;
	while (!no_new_label_PSL) {
		no_new_label_PSL = true;
		for (int u = 0; u < N; u++) {
			results.emplace_back(pool.enqueue([u, g_it] {
				propagate(g_it, 0, u);
				propagate(g_it, 1, u);
				return 1; }));
		}
		for (auto&& result : results) {
			result.get();
		}
		results.clear();
		for (int u = 0; u < N; u++) {
			results.emplace_back(pool.enqueue([u] {
				append(0, u);
				append(1, u);
				return 1; }));
		}
		for (auto&& result : results) {
			result.get();
		}
		results.clear();
	} 

	case_info.time2_PLL_PSL_label_generation = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin1).count() / 1e9;
	auto begin2 = std::chrono::high_resolution_clock::now();
	//---------------------------------------------------------------------------------

	clean_incorrect_labels(N, pool, results); /*clean out incorrect labels and sort*/

	case_info.time3_PLL_PSL_label_postprocess = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin2).count() / 1e9;
	auto begin3 = std::chrono::high_resolution_clock::now();
	//---------------------------------------------------------------------------------

	case_info.label_size_before_canonical_repair = compute_label_bit_size(L_temp_in, L_temp_out);
	if (case_info.use_canonical_repair) {
		canonical_repair_multi_threads(pool, results, &case_info.L_in, &case_info.L_out);
	}
	else {
		case_info.L_in = L_temp_in;
		case_info.L_out = L_temp_out;
	}
	case_info.time4_PLL_PSL_label_canonical = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin3).count() / 1e9;
	//---------------------------------------------------------------------------------

	case_info.label_size_after_canonical_repair = compute_label_bit_size(case_info.L_in, case_info.L_out);
	dgraph_clear_global_values_PLL_PSL();
	case_info.time5_PLL_PSL_total = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time_PSL).count() / 1e9;
}