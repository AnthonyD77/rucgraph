#pragma once
#include <build_in_progress/HL/Hop/graph_v_of_v_idealID_HB_two_hop_labels_v1.h>

/*
    codes for querying distances for hop-bounded
*/

double graph_v_of_v_idealID_two_hop_v1_extract_distance_no_reduc(vector<vector<two_hop_label_v1>>& L, int source, int terminal, int hop_cst, double value_M) {
    /*return std::numeric_limits<double>::max() is not connected*/
    if (hop_cst < 0) {
        return std::numeric_limits<double>::max();
    }
    if (source == terminal) {
        return 0;
    } else if (hop_cst == 0) {
        return std::numeric_limits<double>::max();
    }

    double distance = std::numeric_limits<double>::max();  // if disconnected, return this large value
    auto vector1_check_pointer = L[source].begin();
    auto vector2_check_pointer = L[terminal].begin();
    auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
    while (vector1_check_pointer != pointer_L_s_end) {
        vector2_check_pointer = L[terminal].begin();
        while (vector2_check_pointer != pointer_L_t_end) {
            if (vector2_check_pointer->vertex > vector1_check_pointer->vertex)
                break;
            if (vector1_check_pointer->vertex == vector2_check_pointer->vertex) {
                if (value_M != 0) {
                    int tmp_hop = int(vector1_check_pointer->distance / value_M) + int(vector2_check_pointer->distance / value_M);
                    if (tmp_hop <= hop_cst) {
                        double dis = vector1_check_pointer->distance + vector2_check_pointer->distance - value_M * tmp_hop;
                        if (distance > dis) {
                            distance = dis;
                        }
                    }
                } else {
                    // cout << "hop_cst: " << hop_cst << endl;
                    // cout << vector1_check_pointer->hop << "--" << vector2_check_pointer->hop << endl;
                    if (vector1_check_pointer->hop + vector2_check_pointer->hop <= hop_cst) {
                        double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
                        // cout << "dis: " << dis << endl;
                        if (distance > dis) {
                            distance = dis;
                        }
                    }
                }
            }
            vector2_check_pointer++;
        }
        vector1_check_pointer++;
    }

    return distance;
}

double graph_v_of_v_idealID_two_hop_v1_extract_distance_st_no_R1(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, int source, int terminal, int hop_cst, double value_M = 0) {
    if (hop_cst < 0) {
        return std::numeric_limits<double>::max();
    }
    if (source == terminal) {
        return 0;
    } else if (hop_cst == 0) {
        return std::numeric_limits<double>::max();
    }

    double min_selected_distance = std::numeric_limits<double>::max();

    if (reduction_measures_2019R2[source] == 2) {
        if (reduction_measures_2019R2[terminal] == 2) {
            /*"Both reduced"*/
            // cout << "case 1" << endl;
            auto s_adj_begin = R2_reduced_vertices[source].begin();
            auto s_adj_end = R2_reduced_vertices[source].end();
            auto t_adj_begin = R2_reduced_vertices[terminal].begin();
            auto t_adj_end = R2_reduced_vertices[terminal].end();
            for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++) {
                for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++) {
                    double x = graph_v_of_v_idealID_two_hop_v1_extract_distance_no_reduc(L, it1->first, it2->first, hop_cst - 2, value_M);
                    if (x == std::numeric_limits<double>::max()) {
                        continue;
                    } else {
                        min_selected_distance = min(min_selected_distance, x + double(it1->second) + double(it2->second));
                    }
                }
            }
        } else {
            /*"Only source reduced"*/
            // cout << "case 2" << endl;
            auto s_adj_begin = R2_reduced_vertices[source].begin();
            auto s_adj_end = R2_reduced_vertices[source].end();
            for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++) {
                double x = graph_v_of_v_idealID_two_hop_v1_extract_distance_no_reduc(L, it1->first, terminal, hop_cst - 1, value_M);
                // cout << "dist from " << terminal << " to " << it1->first << ": " << x << " with hop_cst " << hop_cst - 1 << endl;
                if (x == std::numeric_limits<double>::max()) {
                    continue;
                } else {
                    min_selected_distance = min(min_selected_distance, x + double(it1->second));
                }
            }
        }
    } else {
        if (reduction_measures_2019R2[terminal] == 2) {
            /*"Only terminal reduced"*/
            // cout << "case 3" << endl;
            auto t_adj_begin = R2_reduced_vertices[terminal].begin();
            auto t_adj_end = R2_reduced_vertices[terminal].end();
            for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++) {
                double x = graph_v_of_v_idealID_two_hop_v1_extract_distance_no_reduc(L, source, it2->first, hop_cst - 1, value_M);
                // cout << "neighber " << it2->first << " with distance " << x << endl;
                if (x == std::numeric_limits<double>::max()) {
                    continue;
                } else {
                    min_selected_distance = min(min_selected_distance, x + double(it2->second));
                }
            }
        } else {
            /*"Nothing happened"*/
            // cout << "case 4" << endl;
            min_selected_distance = min(min_selected_distance, graph_v_of_v_idealID_two_hop_v1_extract_distance_no_reduc(L, source, terminal, hop_cst, value_M));
        }
    }

    return min_selected_distance;
}

/*
    codes for querying paths for hop-bounded
*/
vector<pair<int, int>> graph_v_of_v_idealID_two_hop_v1_extract_shortest_path_st_no_R1(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, int source, int terminal, int hop_cst, double value_M = 0) {
    vector<pair<int, int>> paths;
    if (source == terminal) {
        return paths;
    }

    double min_dis = std::numeric_limits<double>::max();
    vector<pair<int, int>> partial_edges(2);

    if (reduction_measures_2019R2[source] == 2) {
        auto s_adj_begin = R2_reduced_vertices[source].begin();
        auto s_adj_end = R2_reduced_vertices[source].end();
        /*"Both reduced"*/
        if (reduction_measures_2019R2[terminal] == 2) {
            auto t_adj_begin = R2_reduced_vertices[terminal].begin();
            auto t_adj_end = R2_reduced_vertices[terminal].end();
            for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++) {
                for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++) {
                    double x = graph_v_of_v_idealID_two_hop_v1_extract_distance_no_reduc(L, it1->first, it2->first, hop_cst - 2, value_M) + double(it1->second) + double(it2->second);
                    /*After removing the two edges, it becomes the fourth case: nothing happened*/
                    if (min_dis > x) {
                        min_dis = x;
                        partial_edges[0] = {it1->first, source};
                        partial_edges[1] = {it2->first, terminal};
                    }
                }
            }

            if (min_dis == std::numeric_limits<double>::max()) { /* disconnected */
                return paths;
            }
            paths.push_back(partial_edges[0]);
            paths.push_back(partial_edges[1]);

            /* turn to case 4 */
            vector<pair<int, int>> new_edges;
            if (value_M != 0)
                new_edges = graph_v_of_v_idealID_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, partial_edges[0].first, partial_edges[1].first, hop_cst - 2, value_M);
            else
                new_edges = graph_v_of_v_idealID_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, partial_edges[0].first, partial_edges[1].first, hop_cst - 2);
            if (new_edges.size() > 0) {
                for (int i = new_edges.size() - 1; i >= 0; i--) {
                    paths.push_back(new_edges[i]);
                }
            }
        }
        /*"Only source reduced"*/
        else {
            // cout << "????????????????" << endl;
            for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++) {
                double x = graph_v_of_v_idealID_two_hop_v1_extract_distance_no_reduc(L, it1->first, terminal, hop_cst - 1, value_M) + double(it1->second) - value_M;
                // cout << "neighber " << it1->first << ", dist = " << x << endl;
                // cout << "x: " << x << endl;
                if (min_dis > x) {
                    min_dis = x;
                    partial_edges[0] = {it1->first, source};
                }
            }
            if (min_dis == std::numeric_limits<double>::max()) {
                return paths;
            }
            // cout << "push edge " << partial_edges[0].first << "," << partial_edges[0].second << endl;
            paths.push_back(partial_edges[0]);

            // cout << "hop_cst: " << hop_cst << endl;

            /* turn to case 4 */
            vector<pair<int, int>> new_edges;
            if (value_M != 0)
                new_edges = graph_v_of_v_idealID_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, partial_edges[0].first, terminal, hop_cst - 1, value_M);
            else
                new_edges = graph_v_of_v_idealID_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, partial_edges[0].first, terminal, hop_cst - 1);
            if (new_edges.size() > 0) {
                for (int i = new_edges.size() - 1; i >= 0; i--) {
                    paths.push_back(new_edges[i]);
                }
            }
        }
    } else {
        /*"Only terminal reduced"*/
        if (reduction_measures_2019R2[terminal] == 2) {
            auto t_adj_begin = R2_reduced_vertices[terminal].begin();
            auto t_adj_end = R2_reduced_vertices[terminal].end();
            for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++) {
                double x = graph_v_of_v_idealID_two_hop_v1_extract_distance_no_reduc(L, source, it2->first, hop_cst - 1, value_M) + double(it2->second);
                if (min_dis > x) {
                    min_dis = x;
                    partial_edges[0] = {it2->first, terminal};
                }
            }
            if (min_dis == std::numeric_limits<double>::max()) {
                return paths;
            }
            paths.push_back(partial_edges[0]);

            /* turn to case 4 */
            vector<pair<int, int>> new_edges;
            if (value_M != 0)
                new_edges = graph_v_of_v_idealID_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, source, partial_edges[0].first, hop_cst - 1, value_M);
            else
                new_edges = graph_v_of_v_idealID_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, source, partial_edges[0].first, hop_cst - 1);

            if (new_edges.size() > 0) {
                for (int i = new_edges.size() - 1; i >= 0; i--) {
                    paths.push_back(new_edges[i]);
                }
            }
        }
        /* Nothing happened */
        /* In this case, the problem that the removed vertices appear in the path needs to be solved */
        else {
            int vector1_capped_v_parent = 0, vector2_capped_v_parent = 0;
            double distance = std::numeric_limits<double>::max();
            bool connected = false;
            auto vector1_check_pointer = L[source].begin();
            auto vector2_check_pointer = L[terminal].begin();
            auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
            // while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)

            // 	if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
            // 	{
            // 		connected = true;
            // 		double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
            // 		if (distance > dis)
            // 		{
            // 			distance = dis;
            // 			vector1_capped_v_parent = vector1_check_pointer->parent_vertex;
            // 			vector2_capped_v_parent = vector2_check_pointer->parent_vertex;
            // 		}
            // 		vector1_check_pointer++;
            // 	}
            // 	else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex)
            // 	{
            // 		vector2_check_pointer++;
            // 	}
            // 	else
            // 	{
            // 		vector1_check_pointer++;
            // 	}
            // }
            while (vector1_check_pointer != pointer_L_s_end) {
                vector2_check_pointer = L[terminal].begin();
                while (vector2_check_pointer != pointer_L_t_end) {
                    if (vector2_check_pointer->vertex > vector1_check_pointer->vertex)
                        break;
                    if (vector1_check_pointer->vertex == vector2_check_pointer->vertex) {
                        if (value_M != 0) {
                            int tmp_hop = int(vector1_check_pointer->distance / value_M) + int(vector2_check_pointer->distance / value_M);
                            if (tmp_hop <= hop_cst) {
                                connected = true;
                                double dis = vector1_check_pointer->distance + vector2_check_pointer->distance - value_M * tmp_hop;
                                if (distance > dis) {
                                    distance = dis;
                                    vector1_capped_v_parent = vector1_check_pointer->parent_vertex;
                                    vector2_capped_v_parent = vector2_check_pointer->parent_vertex;
                                }
                            }
                        } else {
                            // cout << hop_cst << endl;
                            if (vector1_check_pointer->hop + vector2_check_pointer->hop <= hop_cst) {
                                connected = true;
                                double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
                                // cout << "dis: " << dis << endl;
                                if (distance > dis) {
                                    distance = dis;
                                    vector1_capped_v_parent = vector1_check_pointer->parent_vertex;
                                    vector2_capped_v_parent = vector2_check_pointer->parent_vertex;
                                }
                            }
                        }
                    }
                    vector2_check_pointer++;
                }
                vector1_check_pointer++;
            }

            // cout << "dis = " << distance << endl;
            // cout << "vector1_capped_v_parent: " << vector1_capped_v_parent << endl;
            // cout << "vector2_capped_v_parent: " << vector2_capped_v_parent << endl;

            if (connected) {
                if (source != vector1_capped_v_parent) {
                    paths.push_back({source, vector1_capped_v_parent});
                    source = vector1_capped_v_parent;
                    hop_cst--;
                }
                if (terminal != vector2_capped_v_parent) {
                    paths.push_back({terminal, vector2_capped_v_parent});
                    terminal = vector2_capped_v_parent;
                    hop_cst--;
                }
            } else {
                return paths;
            }

            /* find new */
            // cout << "turn to " << source << "-" << terminal << endl;
            vector<pair<int, int>> new_edges;
            if (value_M != 0)
                new_edges = graph_v_of_v_idealID_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, source, terminal, hop_cst, value_M);
            else
                new_edges = graph_v_of_v_idealID_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, source, terminal, hop_cst);

            if (new_edges.size() > 0) {
                for (int i = new_edges.size() - 1; i >= 0; i--) {
                    paths.push_back(new_edges[i]);
                }
            }
        }
    }

    return paths;
}