#pragma once
#include <build_in_progress/HL/dynamic/PLL_dynamic.h>
#include <build_in_progress/HL/batch/WeightChangeMaintenance_batch_base.h>
#include <algorithm>
#include <build_in_progress/HL/dynamic/WeightDecreaseMaintenance_improv_multiThread.h>

bool Sort_Affected_Label(const affected_label& al1, const affected_label& al2) {
    if (al1.second != al2.second) {
        return al1.second > al2.second;
    }
    else {
        if (al1.dis != al2.dis) {
            return al1.dis < al2.dis;
        } 
        else {
            return true;
        }
    }
}

void UpdateLabels(graph_v_of_v_idealID& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, std::priority_queue<pair<weightTYPE,int> >& Q, int v) {
    vector<int> shortest_dist_changed, dist_changed;
    auto& shortest_dist = Dis_batch[0];
    auto& DIS = Dis_batch[0];
    auto& dist = Dis_batch[1]; // can use weightTYPE instead of pair
    auto& L = mm.L; 
    auto& PPR = mm.PPR;
    auto& Lv = mm.L[v];

    #ifdef DEBUG_MODE
        cerr << "label_v: " << v << endl;
    #endif

    while (!Q.empty()) {
        auto temp = Q.top();
        Q.pop();
        int x = temp.second;
        weightTYPE dx = -temp.first;
        if (shortest_dist[x].first == -1) {
            shortest_dist[x] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4(L[x], Lv);
            shortest_dist_changed.push_back(x);
            #ifdef DEBUG_MODE
                cerr << "source: " << x << " " << "shortest dis: " << DIS[x].first << endl;
            #endif

            auto search_result = search_sorted_two_hop_label2(L[x], v);
            if (search_result.second != -1) { 
                dist[x].first = search_result.first;
                dist_changed.push_back(x);
            }
        }

        #ifdef DEBUG_MODE
            cerr << "Q_element: " << "x: " << x << " " << "dx: " << dx << " DIS[x].first: " << DIS[x].first << endl;
        #endif

        if (dx <= shortest_dist[x].first + 1e-5) { //dx <= sd
            shortest_dist[x] = { dx, v };
            shortest_dist_changed.push_back(x);
            insert_sorted_two_hop_label(L[x], v, dx);
            dist[x].first = dx;
            dist_changed.push_back(x);
        }
        else if (dist[x].first != -1 && dx <= dist[x].first + 1e-5) {
            dist[x].first = dx;
            dist_changed.push_back(x);
            insert_sorted_two_hop_label(L[x], v, dx);
        }

        for (auto& nei : instance_graph[x]) {
            int xnei = nei.first;
            weightTYPE d_new = dx + nei.second;
            if (v >= xnei) continue;

            if(shortest_dist[xnei].first == -1){
                shortest_dist[xnei] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4(L[xnei], Lv);
                shortest_dist_changed.push_back(xnei);
                #ifdef DEBUG_MODE
                    cerr << "source: " << xnei << " " << "shortest dist: " << DIS[xnei].first << endl;
                #endif

                auto search_result = search_sorted_two_hop_label2(L[xnei], v);
                if (search_result.second != -1) {
                    dist[xnei].first = search_result.first;
                    dist_changed.push_back(xnei);
                }
            }
            
            if (shortest_dist[xnei].first > d_new + 1e-5) {
                shortest_dist[xnei] = { d_new, v };
                if (dist[xnei].first == -1) {
                    dist_changed.push_back(xnei);
                }
                dist[xnei].first = d_new;
                Q.push({-d_new, xnei});
            }
            else if (dist[xnei].first != -1) { //redundant label
                if (dist[xnei].first > d_new + 1e-5) {
                    dist[xnei].first = d_new;
                    Q.push({-d_new, xnei});
                }
            }
            else {
                //PPL
                if (shortest_dist[xnei].second != v) {
                    PPR_insert(PPR, xnei, shortest_dist[xnei].second, v);
                }
                if (shortest_dist[xnei].second != xnei) {
                    PPR_insert(PPR, v, shortest_dist[xnei].second, xnei);
                }
            }
        }
    }
    for (int i : shortest_dist_changed) {
        shortest_dist[i] = { -1, -1 };
    }
    for (int i : dist_changed) {
        dist[i] = {-1, -1};
    }
}

void WeightDecrease_batch_work(graph_v_of_v_idealID& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, std::vector<graph_edge_change>& edge_changes) {

    std::vector<affected_label> CL;
    auto& L = mm.L;

    for (auto edge_change : edge_changes) {
        int v1 = edge_change.source;
        int v2 = edge_change.terminal;
        weightTYPE w_new = edge_change.new_weight;

        #ifdef DEBUG_MODE
            cerr << v1 << " " << v2 << " " << w_new << endl;
        #endif
        
        for (int sl = 0; sl < 2; sl++) {
            if (sl == 1) {
                swap(v1, v2);
            }
            for (auto it : mm.L[v1]) {
                if (it.vertex <= v2) {
                    auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(L, it.vertex, v2); // query_result is {distance, common hub}
                    if (query_result.first > it.distance + w_new) {
                        CL.push_back(affected_label{ v2 , it.vertex, it.distance + w_new });
                    }
					else {
						auto search_result = search_sorted_two_hop_label2(L[v2], it.vertex);
						if (search_result.second != -1 && search_result.first > it.distance + w_new) {
							CL.push_back(affected_label{ v2, it.vertex, it.distance + w_new });
						}
					}
                }
            }
        }
    }
    if (CL.empty()) return;

    #ifdef DEBUG_MODE
        cerr << "CL.size(): " << CL.size() << endl;
    #endif
    
    sort(CL.begin(), CL.end(), Sort_Affected_Label);
    std::priority_queue<pair<weightTYPE,int> > Q;
    int last_label = -1;
    #ifdef DEBUG_MODE
        cerr << "----------Affected labels----------" << endl;
    #endif   
    for (auto al : CL) {
        if (al.second != last_label) {
            UpdateLabels(instance_graph, mm, Q, last_label);
            #ifdef DEBUG_MODE
                cerr << "--------------------------------" << endl;
            #endif
            last_label = al.second;
            #ifdef DEBUG_MODE
                if (!Q.empty()) cerr << "ERROR!" << endl;
            #endif
        }
        #ifdef DEBUG_MODE
            cerr << al.first << " " << al.second << " " << al.dis << endl;
        #endif
        Q.push({-al.dis, al.first});
    }
    UpdateLabels(instance_graph, mm, Q, last_label);
    #ifdef DEBUG_MODE
        cerr << "-----------------------------" << endl;
    #endif
}

void WeightDecrease_batch(graph_v_of_v_idealID& ideal_g, graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	std::vector<graph_edge_change>& edge_changes, double& avg_maintain_time) {
        
        for (auto it : edge_changes) {
            graph_hash_of_mixed_weighted_add_edge(instance_graph, it.source, it.terminal, it.new_weight);
            graph_v_of_v_idealID_add_edge(ideal_g, it.source, it.terminal, it.new_weight);
        }
        auto begin = std::chrono::high_resolution_clock::now();
        WeightDecrease_batch_work(ideal_g, mm, edge_changes);
        auto end = std::chrono::high_resolution_clock::now();
        avg_maintain_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
}

void WeightDecrease_batch_original(graph_v_of_v_idealID& ideal_g, graph_hash_of_mixed_weighted& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm,
	std::vector<graph_edge_change>& edge_changes, int thread_num, double& avg_maintain_time) {
        
        ThreadPool pool_dynamic(thread_num);
	    std::vector<std::future<int>> results_dynamic;

        for (auto it : edge_changes) {
            graph_hash_of_mixed_weighted_add_edge(instance_graph, it.source, it.terminal, it.new_weight); 
            graph_v_of_v_idealID_add_edge(ideal_g, it.source, it.terminal, it.new_weight);

            auto begin = std::chrono::high_resolution_clock::now();

            /*maintain labels*/
            //WeightIncrease2021(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight, pool_dynamic, results_dynamic);
            WeightDecreaseMaintenance_improv(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight, pool_dynamic, results_dynamic);
            //WeightIncrease2019(ideal_g, mm, it.source, it.terminal, it.old_weight, it.new_weight, pool_dynamic, results_dynamic, 1e1);
            auto end = std::chrono::high_resolution_clock::now();
            avg_maintain_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
        }
}