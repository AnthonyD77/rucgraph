#pragma once
#include <build_in_progress/HL/dynamic/PLL_dynamic.h>
#include <build_in_progress/HL/batch/WeightChangeMaintenance_batch_base.h>
#include <queue>

void SPREAD1(graph_v_of_v_idealID& instance_graph, vector<vector<two_hop_label_v1>>* L,
	std::vector<affected_label>& al1, std::vector<pair_label>* al2) {

	for (auto it : al1) {
        queue<pair<int, weightTYPE> > q; //(u,d)
        int v = it.second;
        q.push(pair<int, weightTYPE>(it.first, it.dis));
        while (!q.empty()) {
            int x = q.front().first;
            weightTYPE dx = q.front().second;
            q.pop();
            insert_sorted_two_hop_label((*L)[x], v, MAX_VALUE); // this does not change the size of L[x] here, so does not need to lock here
            al2->push_back(pair_label(x, v));
            for (auto nei : instance_graph[x]) {
                if (v < nei.first) {
                    weightTYPE search_weight = search_sorted_two_hop_label((*L)[nei.first], v);
                    if (abs(dx + nei.second - search_weight) < 1e-5) {
                        q.push(pair<int, weightTYPE>(nei.first, dx + nei.second));
                    }
                }
            }
        }
	}
}

void SPREAD2(graph_v_of_v_idealID& instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR,
	std::vector<pair_label>& al2, std::vector<affected_label>* al3) {

	for (auto it : al2) {
        int v = it.first, u = it.second;
        std::vector<int> temp = PPR_retrieve(*PPR, v, u);
        PPR_binary_operations_insert(temp, u);
        for (auto t : temp) {
            if (v < t) {
                weightTYPE d1 = MAX_VALUE;
                for (auto nei : instance_graph[t]) {
                    d1 = min(d1, search_sorted_two_hop_label((*L)[nei.first], v) + (weightTYPE)nei.second);
                }
                auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, t, v);
                if (query_result.first > d1 + 1e-5) { // only add new label when it's absolutely necessary
                    al3->push_back(affected_label(t, v, d1));
                    //cout<<"spread21: "<<t<<' '<<v << ' ' << d1 << ' ' << query_result.first << endl;
                }
                else {
                    if (query_result.second != v) {
                        PPR_insert(*PPR, t, query_result.second, v);
                    }
                    if (query_result.second != t) {
                        PPR_insert(*PPR, v, query_result.second, t);
                    }
                }
            }
            else if (t < v) {
                weightTYPE d1 = MAX_VALUE;
                for (auto nei : instance_graph[v]) {
                    d1 = min(d1, search_sorted_two_hop_label((*L)[nei.first], t) + (weightTYPE)nei.second);
                }
                auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, v, t);
                if (query_result.first > d1 + 1e-5) {
                    al3->push_back(affected_label(v, t, d1));
                    //cout << "spread22: " << v << ' ' << t << ' ' << d1 << ' ' << query_result.first << endl;
                }
                else {
                    if (query_result.second != v) {
                        PPR_insert(*PPR, t, query_result.second, v);
                    }
                    if (query_result.second != t) {
                        PPR_insert(*PPR, v, query_result.second, t);
                    }
                }
            }
        }
	}
}

void SPREAD3(graph_v_of_v_idealID& instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR, std::vector<affected_label>& al3) {

	for (auto it : al3) {

        int u = it.first, v = it.second;
        weightTYPE du = it.dis;

        auto Lv = (*L)[v]; // to avoid interlocking

        auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4((*L)[u], Lv);

        if (query_result.first + 1e-5 < du) {
            if (query_result.second != v) {
                PPR_insert(*PPR, u, query_result.second, v);
            }
            if (query_result.second != u) {
                PPR_insert(*PPR, v, query_result.second, u);
            }
        }

        vector<int> Dis_changed;
        auto& DIS = Dis[0];
        auto& Q_HANDLES = Q_handles[0];
        auto& Q_VALUE = Q_value[0];
        DIS[u] = { du, v }; // <distance, hub responsible for this distance>
        Dis_changed.push_back(u);

        boost::heap::fibonacci_heap<node_for_DIFFUSE> pq;
        Q_HANDLES[u] = pq.push(node_for_DIFFUSE(u, du));
        Q_VALUE[u] = du;

        //cout << "spread3 0: " << u << ' ' << v << " " << du << endl;
        while (!pq.empty()) {
            int x = pq.top().index;
            weightTYPE dx = pq.top().disx;
            pq.pop();
            Q_VALUE[x] = MAX_VALUE;

            if (dx < search_sorted_two_hop_label((*L)[x], v)) {
                insert_sorted_two_hop_label((*L)[x], v, dx);
            }

            for (auto nei : instance_graph[x]) {
                int xnei = nei.first;
                weightTYPE d_new = dx + nei.second;

                if (v < xnei) {
                    if (DIS[xnei].first == -1) {
                        DIS[xnei] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4((*L)[xnei], Lv);
                        Dis_changed.push_back(xnei);
                    }
                    if (DIS[xnei].first > d_new + 1e-5) {
                        DIS[xnei] = { d_new, v };
                        if (Q_VALUE[xnei] == MAX_VALUE) {
                            Q_HANDLES[xnei] = pq.push(node_for_DIFFUSE(xnei, d_new));
                        }
                        else {
                            pq.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
                        }
                        Q_VALUE[xnei] = d_new;
                    }
                    else {
                        if (DIS[xnei].second != v) {
                            PPR_insert(*PPR, xnei, DIS[xnei].second, v);
                        }
                        if (DIS[xnei].second != xnei) {
                            PPR_insert(*PPR, v, DIS[xnei].second, xnei);
                        }
                    }
                }
            }
        }

        for (int i : Dis_changed) {
            DIS[i] = { -1, -1 };
        }
	}
}

void WeightIncreaseMaintenance_improv(graph_v_of_v_idealID& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new) {

	std::vector<affected_label> al1, al3;
	std::vector<pair_label> al2;

	/*it's slow to paralize the following part*/
	for (auto it : mm.L[v1]) {
		if (it.vertex <= v2 && abs(search_sorted_two_hop_label(mm.L[v2], it.vertex) - it.distance - w_old) < 1e-5) {
			al1.push_back(affected_label(v2, it.vertex, it.distance + w_old));
		}
	}
	for (auto it : mm.L[v2]) {
		if (it.vertex <= v1 && abs(search_sorted_two_hop_label(mm.L[v1], it.vertex) - it.distance - w_old) < 1e-5) {
			al1.push_back(affected_label(v1, it.vertex, it.distance + w_old));
		}
	}

	//cout << "al1.size() " << al1.size() << endl;

	SPREAD1(instance_graph, &mm.L, al1, &al2);
	SPREAD2(instance_graph, &mm.L, &mm.PPR, al2, &al3);
	SPREAD3(instance_graph, &mm.L, &mm.PPR, al3);

	//for (auto it : al2) {
	//	cout << "al2 " << it.first << " " << it.second << endl;
	//}
	//for (auto it : al3) {
	//	cout << "al3 " << it.first << " " << it.second << " " << it.dis << endl;
	//}
}

void WeightIncrease_batch(graph_v_of_v_idealID& ideal_g, graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	std::vector<graph_edge_change>& edge_changes, double& avg_maintain_time) {
    int cnt = 0;
    for (auto it : edge_changes) {
        graph_hash_of_mixed_weighted_add_edge(instance_graph, it.source, it.terminal, it.new_weight); // increase weight
        graph_v_of_v_idealID_add_edge(ideal_g, it.source, it.terminal, it.new_weight);

        auto begin = std::chrono::high_resolution_clock::now();

        /*maintain labels*/
        //WeightIncrease2021(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight, pool_dynamic, results_dynamic);
        cnt ++;
        #ifdef CHECK_PROCESS
            cerr << cnt << endl;
        #endif
        WeightIncreaseMaintenance_improv(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight);
        #ifdef CHECK_PROCESS
            // cerr << cnt << endl;
        #endif
        //WeightIncrease2019(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight, pool_dynamic, results_dynamic, 1e1);
        auto end = std::chrono::high_resolution_clock::now();
        avg_maintain_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
    }
}

void WeightIncrease_batch_original(graph_v_of_v_idealID& ideal_g, graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	std::vector<graph_edge_change>& edge_changes, int thread_num, double& avg_maintain_time) {
    
    ThreadPool pool_dynamic(thread_num);
	std::vector<std::future<int>> results_dynamic; 
    int cnt = 0;
    for (auto it : edge_changes) {
        graph_hash_of_mixed_weighted_add_edge(instance_graph, it.source, it.terminal, it.new_weight); // increase weight
        graph_v_of_v_idealID_add_edge(ideal_g, it.source, it.terminal, it.new_weight);

        auto begin = std::chrono::high_resolution_clock::now();

        /*maintain labels*/
        //WeightIncrease2021(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight, pool_dynamic, results_dynamic);
        cnt ++;
        #ifdef CHECK_PROCESS
            cerr << cnt << endl;
        #endif
        WeightIncreaseMaintenance_improv(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight, pool_dynamic, results_dynamic);
        #ifdef CHECK_PROCESS
            // cerr << cnt << endl;
        #endif
        //WeightIncrease2019(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight, pool_dynamic, results_dynamic, 1e1);
        auto end = std::chrono::high_resolution_clock::now();
        avg_maintain_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
    }
}