#include <vector>
#include <queue>
#include <algorithm>
#include <math.h>
using namespace std;

#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

void Distance_Dijsktra(graph_hash_of_mixed_weighted& instance_graph, int s, vector<weightTYPE>& d) {
	int n = instance_graph.hash_of_vectors.size();
	vector<bool> mark(n, false);
	d[s] = 0;
	priority_queue<pair<weightTYPE, int>, vector<pair<weightTYPE, int> >, greater<pair<weightTYPE, int> > > Q;
	Q.push(pair<weightTYPE, int>(0, s));
	while (!Q.empty()) {
		int v = Q.top().second;
		Q.pop();
		if (mark[v]) continue;
		mark[v] = true;
		for (auto nei : instance_graph.adj_v_and_ec(v)) {
			if (d[nei.first] > d[v] + nei.second) {
				d[nei.first] = d[v] + nei.second;
				Q.push(pair<weightTYPE, int>(d[nei.first], nei.first));
			}
		}
	}
}

void FindAffectedNode(graph_hash_of_mixed_weighted& instance_graph,
	graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int x, int y, weightTYPE w_old, vector<int>& A) {

	int n = instance_graph.hash_of_vectors.size();
	vector<bool> mark(n, false);

	vector<weightTYPE> dy(n, MAX_VALUE);
	Distance_Dijsktra(instance_graph, y, dy);
	queue<int> Q;
	Q.push(x);
	while (!Q.empty()) {
		int v = Q.front();
		Q.pop();
		A.push_back(v);
		for (auto u : instance_graph.adj_v(v)) {
			if (mark[u]) continue;
			weightTYPE dy_old = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(mm.L, u, y);
			weightTYPE dy_new = dy[u];
			weightTYPE dx_old = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(mm.L, u, x);
			if (abs(dy_old - dy_new) > 1e-5) {
				mark[u] = true;
				Q.push(u);
			}
			else {
				int h = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(mm.L, u, y).second;
				if (find(A.begin(), A.end(), h) != A.end() || (h == u || h == y) && abs(dy_old - (dx_old + w_old)) <= 1e-5) {
					mark[u] = true;
					Q.push(u);
				}
			}
		}
	}
}

void RemoveAffectedHub(graph_hash_of_mixed_weighted& instance_graph,
	graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, vector<int>& AFF_x, vector<int>& AFF_y,
	vector<bool>& ax, vector<bool>& ay) {
	for (auto v : AFF_x) {
		for (auto it = mm.L[v].begin(); it != mm.L[v].end();) {
			if (ay[it->vertex]) {
				it = mm.L[v].erase(it);
			}
			else {
				it++;
			}
		}
	}
	for (auto v : AFF_y) {
		for (auto it = mm.L[v].begin(); it != mm.L[v].end();) {
			if (ax[it->vertex]) {
				it = mm.L[v].erase(it);
			}
			else {
				it++;
			}
		}
	}
}

void GreedyRestore(graph_hash_of_mixed_weighted& instance_graph,
	graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, vector<int>& AFF_x, vector<int>& AFF_y,
	vector<bool>& ax, vector<bool>& ay) {
	int n = instance_graph.hash_of_vectors.size();
	vector<int>& SA = AFF_x;
	vector<int>& LA = AFF_y;
	for (auto a : SA) {
		vector<weightTYPE> dist(n, MAX_VALUE);
		priority_queue<pair<weightTYPE, int>, vector<pair<weightTYPE, int> >, greater<pair<weightTYPE, int> > > Q;
		dist[a] = 0;
		Q.push(pair<weightTYPE, int>(0, a));
		while (!Q.empty()) {
			int v = Q.top().second;
			Q.pop();
			if (ay[v]) {
				weightTYPE qdist = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(mm.L, v, a);
				if (dist[v] - 1e-5 < qdist) {
					if (v < a) insert_sorted_two_hop_label(mm.L[a], v, dist[v]);
					else insert_sorted_two_hop_label(mm.L[v], a, dist[v]);
				}
			}
			for (auto u : instance_graph.adj_v_and_ec(v)) {
				if (dist[u.first] - 1e-5 < dist[v] + u.second) continue;
				dist[u.first] = dist[v] + u.second;
				Q.push(pair<weightTYPE, int>(dist[u.first], u.first));
			}
		}
	}
}

void OrderRestore(graph_hash_of_mixed_weighted& instance_graph,
	graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, vector<int>& AFF_x, vector<int>& AFF_y,
	vector<bool>& ax, vector<bool>& ay) {
	int n = instance_graph.hash_of_vectors.size();
	vector<int>& FA = AFF_x;
	FA.insert(FA.end(), AFF_y.begin(), AFF_y.end());
	sort(FA.begin(), FA.end());
	for (auto a : FA) {
		vector<weightTYPE> dist(n, MAX_VALUE);
		priority_queue<pair<weightTYPE, int>, vector<pair<weightTYPE, int> >, greater<pair<weightTYPE, int> > > Q;
		dist[a] = 0;
		Q.push(pair<weightTYPE, int>(0, a));
		while (!Q.empty()) {
			int v = Q.top().second;
			Q.pop();
			if (v < a) continue;
			if ((ay[v] && ax[a]) || (ay[a] && ax[v])) {
				weightTYPE qdist = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(mm.L, v, a);
				if (dist[v] - 1e-5 < qdist) {
					insert_sorted_two_hop_label(mm.L[v], a, dist[v]);
				}
			}
			for (auto u : instance_graph.adj_v_and_ec(v)) {
				if (dist[u.first] - 1e-5 < dist[v] + u.second) continue;
				dist[u.first] = dist[v] + u.second;
				Q.push(pair<weightTYPE, int>(dist[u.first], u.first));
			}
		}
	}
}

void WeightIncrease2019(graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int x, int y, weightTYPE w_old) {

	vector<int> AFF_x, AFF_y;
	vector<bool> ax(instance_graph.hash_of_vectors.size(), false);
	vector<bool> ay(instance_graph.hash_of_vectors.size(), false);

	FindAffectedNode(instance_graph, mm, x, y, w_old, AFF_x);
	FindAffectedNode(instance_graph, mm, y, x, w_old, AFF_y);
	sort(AFF_x.begin(), AFF_x.end());
	AFF_x.erase(unique(AFF_x.begin(), AFF_x.end()), AFF_x.end());
	sort(AFF_y.begin(), AFF_y.end());
	AFF_y.erase(unique(AFF_y.begin(), AFF_y.end()), AFF_y.end());

	for (int i = 0; i < AFF_x.size(); i++) {
		ax[AFF_x[i]] = true;
	}
	for (int i = 0; i < AFF_y.size(); i++) {
		ay[AFF_y[i]] = true;
	}

	RemoveAffectedHub(instance_graph, mm, AFF_x, AFF_y, ax, ay);

	double small_size = min(AFF_x.size(), AFF_y.size());
	double n = instance_graph.hash_of_vectors.size();

	if (small_size > n / log(n)) {
		if (AFF_x.size() < AFF_y.size())
			GreedyRestore(instance_graph, mm, AFF_x, AFF_y, ax, ay);
		else
			GreedyRestore(instance_graph, mm, AFF_y, AFF_x, ay, ax);
	}
	else {
		OrderRestore(instance_graph, mm, AFF_x, AFF_y, ax, ay);
	}

}