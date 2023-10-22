#include <vector>
#include <queue>
#include <algorithm>
#include <math.h>
using namespace std;

#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

void Distance_Dijsktra(graph_v_of_v_idealID& instance_graph, int s, vector<weightTYPE>& d) {
	int n = instance_graph.size();
	vector<bool> mark(n, false);
	d[s] = 0;
	priority_queue<pair<weightTYPE, int>, vector<pair<weightTYPE, int> >, greater<pair<weightTYPE, int> > > Q;
	Q.push(pair<weightTYPE, int>(0, s));
	while (!Q.empty()) {
		int v = Q.top().second;
		Q.pop();
		if (mark[v]) continue;
		mark[v] = true;

		for (auto nei : instance_graph[v]) {
			if (d[nei.first] > d[v] + nei.second) {
				d[nei.first] = d[v] + nei.second;
				Q.push(pair<weightTYPE, int>(d[nei.first], nei.first));
			}
		}
	}
}

void FindAffectedNode(graph_v_of_v_idealID& instance_graph,
	graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int x, int y, weightTYPE w_old, vector<int>& A, vector<bool>& a) {

	int n = instance_graph.size();
	vector<bool> mark(n, false);

	vector<weightTYPE> dy(n, MAX_VALUE);
	Distance_Dijsktra(instance_graph, y, dy);
	queue<int> Q;
	Q.push(x);
	while (!Q.empty()) {

		if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
			throw reach_limit_time_string;
		}

		int v = Q.front();
		Q.pop();
		A.push_back(v);
		a[v] = true;
		for (auto ux : instance_graph[v]) {
			int u = ux.first;
			if (mark[u]) continue;
			auto query_dy_old = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(mm.L, u, y);
			weightTYPE dy_old = query_dy_old.first;
			weightTYPE dx_old = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(mm.L, u, x);
			if (abs(dy_old - dy[u]) > 1e-5) {
				mark[u] = true;
				Q.push(u);
			}
			else {
				int h = query_dy_old.second;
				if (a[h] || (h == u || h == y) && abs(dy_old - (dx_old + w_old)) <= 1e-5) {
					mark[u] = true;
					Q.push(u);
				}
			}
		}
	}
}

void RemoveAffectedHub(graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, vector<int>& AFF_x, vector<int>& AFF_y, vector<bool>& ax, vector<bool>& ay,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {
	for (auto v : AFF_x) {
		results_dynamic.emplace_back(pool_dynamic.enqueue([v, &mm, &ay] {

			if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
				return 1;
			}

			vector<two_hop_label_v1> Lv;
			for (auto it = mm.L[v].begin(); it != mm.L[v].end(); it++) {
				if (ay[it->vertex] == 0) {
					Lv.push_back(*it);
				}
			}
			mm.L[v] = Lv;


			//for (auto it = mm.L[v].begin(); it != mm.L[v].end();) {
			//	if (ay[it->vertex]) {
			//		it = mm.L[v].erase(it); // this is likely to be faster than building a new vector with not erased elements
			//	}
			//	else {
			//		it++;
			//	}
			//}


			return 1; }));
	}
	for (auto&& result : results_dynamic) {
		result.get();
	}
	results_dynamic.clear();
	for (auto v : AFF_y) {
		results_dynamic.emplace_back(pool_dynamic.enqueue([v, &mm, &ax] {

			if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
				return 1;
			}

			vector<two_hop_label_v1> Lv;
			for (auto it = mm.L[v].begin(); it != mm.L[v].end(); it++) {
				if (ax[it->vertex] == 0) {
					Lv.push_back(*it);
				}
			}
			mm.L[v] = Lv;


			//for (auto it = mm.L[v].begin(); it != mm.L[v].end();) {
			//	if (ax[it->vertex]) {
			//		it = mm.L[v].erase(it);
			//	}
			//	else {
			//		it++;
			//	}
			//}


			return 1; }));
	}
	for (auto&& result : results_dynamic) {
		result.get();
	}
	results_dynamic.clear();
}

void GreedyRestore(graph_v_of_v_idealID& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, vector<int>& AFF_x,
	vector<bool>& ay, ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {
	int n = instance_graph.size();
	vector<int>& SA = AFF_x;
	for (auto a : SA) {
		results_dynamic.emplace_back(pool_dynamic.enqueue([n, &instance_graph, a, &mm, &ay] {

			if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
				return 1;
			}

			mtx_595_1.lock();
			int current_tid = Qid_595.front();
			Qid_595.pop();
			mtx_595_1.unlock();

			auto& dist = Q_value[current_tid];
			vector<int> dist_changed;

			priority_queue<pair<weightTYPE, int>, vector<pair<weightTYPE, int> >, greater<pair<weightTYPE, int> > > Q;
			dist[a] = 0;
			dist_changed.push_back(a);
			Q.push(pair<weightTYPE, int>(0, a));
			while (!Q.empty()) {

				if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
					for (int i : dist_changed) {
						dist[i] = MAX_VALUE;
					}
					mtx_595_1.lock();
					Qid_595.push(current_tid);
					mtx_595_1.unlock();
					return 1;
				}

				int v = Q.top().second;
				Q.pop();
				if (ay[v]) {
					mtx_595[a].lock();
					auto L_a = mm.L[a];
					mtx_595[a].unlock();
					mtx_595[v].lock();
					auto L_v = mm.L[v];
					mtx_595[v].unlock();
					weightTYPE qdist = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc3(L_v, L_a);
					if (dist[v] - 1e-5 < qdist) {
						if (v < a) {
							mtx_595[a].lock();
							insert_sorted_two_hop_label(mm.L[a], v, dist[v]);
							mtx_595[a].unlock();
						}
						else {
							mtx_595[v].lock();
							insert_sorted_two_hop_label(mm.L[v], a, dist[v]);
							mtx_595[v].unlock();
						}
					}
				}

				for (auto u : instance_graph[v]) {
					if (dist[u.first] - 1e-5 < dist[v] + u.second) continue;
					dist[u.first] = dist[v] + u.second;
					dist_changed.push_back(u.first);
					Q.push(pair<weightTYPE, int>(dist[u.first], u.first));
				}
			}

			for (int i : dist_changed) {
				dist[i] = MAX_VALUE;
			}
			mtx_595_1.lock();
			Qid_595.push(current_tid);
			mtx_595_1.unlock();

			return 1; }));
	}
	for (auto&& result : results_dynamic) {
		result.get();
	}
	results_dynamic.clear();
}

void OrderRestore(graph_v_of_v_idealID& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, vector<int>& AFF_x, vector<int>& AFF_y,
	vector<bool>& ax, vector<bool>& ay, ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {
	int n = instance_graph.size();
	vector<int>& FA = AFF_x;
	FA.insert(FA.end(), AFF_y.begin(), AFF_y.end());
	sort(FA.begin(), FA.end());
	for (auto a : FA) {
		results_dynamic.emplace_back(pool_dynamic.enqueue([n, &instance_graph, a, &mm, &ax, &ay] {

			if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
				return 1;
			}

			mtx_595_1.lock();
			int current_tid = Qid_595.front();
			Qid_595.pop();
			mtx_595_1.unlock();

			auto& dist = Q_value[current_tid];
			vector<int> dist_changed;

			priority_queue<pair<weightTYPE, int>, vector<pair<weightTYPE, int> >, greater<pair<weightTYPE, int> > > Q;
			dist[a] = 0;
			dist_changed.push_back(a);
			Q.push(pair<weightTYPE, int>(0, a));
			while (!Q.empty()) {

				if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
					for (int i : dist_changed) {
						dist[i] = MAX_VALUE;
					}
					mtx_595_1.lock();
					Qid_595.push(current_tid);
					mtx_595_1.unlock();
					return 1;
				}

				int v = Q.top().second;
				Q.pop();
				// if (v < a) continue;
				if ((ay[v] && ax[a]) || (ay[a] && ax[v])) {
					mtx_595[a].lock();
					auto L_a = mm.L[a];
					mtx_595[a].unlock();
					mtx_595[v].lock();
					auto L_v = mm.L[v];
					mtx_595[v].unlock();
					weightTYPE qdist = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc3(L_v, L_a);
					if (dist[v] - 1e-5 < qdist) {
						mtx_595[v].lock();
						insert_sorted_two_hop_label(mm.L[v], a, dist[v]);
						mtx_595[v].unlock();
					}
				}

				for (auto u : instance_graph[v]) {
					if (dist[u.first] - 1e-5 < dist[v] + u.second) continue;
					dist[u.first] = dist[v] + u.second;
					dist_changed.push_back(u.first);
					if (u.first < a) continue;
					Q.push(pair<weightTYPE, int>(dist[u.first], u.first));
				}
			}

			for (int i : dist_changed) {
				dist[i] = MAX_VALUE;
			}

			mtx_595_1.lock();
			Qid_595.push(current_tid);
			mtx_595_1.unlock();

			return 1; }));
	}
	for (auto&& result : results_dynamic) {
		result.get();
	}
	results_dynamic.clear();
}

void WeightIncrease2019(graph_v_of_v_idealID& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int x, int y, weightTYPE w_old, weightTYPE w_new,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic, double max_run_second) {

	begin_time = std::chrono::high_resolution_clock::now();
	max_run_time_nanosec = max_run_second * 1e9;

	double n = instance_graph.size();

	vector<int> AFF_x, AFF_y;
	vector<bool> ax(n, false);
	vector<bool> ay(n, false);

	FindAffectedNode(instance_graph, mm, x, y, w_old, AFF_x, ax);
	FindAffectedNode(instance_graph, mm, y, x, w_old, AFF_y, ay);
	sort(AFF_x.begin(), AFF_x.end());
	AFF_x.erase(unique(AFF_x.begin(), AFF_x.end()), AFF_x.end());
	sort(AFF_y.begin(), AFF_y.end());
	AFF_y.erase(unique(AFF_y.begin(), AFF_y.end()), AFF_y.end());

	if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
		throw reach_limit_time_string;
	}

	RemoveAffectedHub(mm, AFF_x, AFF_y, ax, ay, pool_dynamic, results_dynamic);

	if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
		throw reach_limit_time_string;
	}

	double small_size = min(AFF_x.size(), AFF_y.size());

	if (1) { // new strategy
		double all_size = AFF_x.size() + AFF_y.size();
		if (all_size / small_size > 500) {
			if (AFF_x.size() < AFF_y.size())
				GreedyRestore(instance_graph, mm, AFF_x, ay, pool_dynamic, results_dynamic);
			else
				GreedyRestore(instance_graph, mm, AFF_y, ax, pool_dynamic, results_dynamic);
		}
		else {
			if (small_size < n / log(n)) {
				if (AFF_x.size() < AFF_y.size())
					GreedyRestore(instance_graph, mm, AFF_x, ay, pool_dynamic, results_dynamic);
				else
					GreedyRestore(instance_graph, mm, AFF_y, ax, pool_dynamic, results_dynamic);
			}
			else {
				OrderRestore(instance_graph, mm, AFF_x, AFF_y, ax, ay, pool_dynamic, results_dynamic);
			}
		}
	}
	else { // 2019 strategy
		if (small_size < n / log(n)) {
			if (AFF_x.size() < AFF_y.size())
				GreedyRestore(instance_graph, mm, AFF_x, ay, pool_dynamic, results_dynamic);
			else
				GreedyRestore(instance_graph, mm, AFF_y, ax, pool_dynamic, results_dynamic);
		}
		else {
			OrderRestore(instance_graph, mm, AFF_x, AFF_y, ax, ay, pool_dynamic, results_dynamic);
		}
	}

	if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
		throw reach_limit_time_string;
	}
}