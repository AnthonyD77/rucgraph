#pragma once
#include <iostream>
#include <tool_functions/ThreadPool.h>
#include <shared_mutex>
#include <chrono>
#include <boost/heap/fibonacci_heap.hpp>
#include <dgraph_v_of_v/dgraph_v_of_v.h>

void print_label_set(const std::vector<std::vector<two_hop_label>> &_L)
{
    for (int i = 0; i < _L.size(); i++)
    {
        std::cout << "[" << i << "] ";
        for (auto &cur : _L[i])
        {
            std::cout << "(" << cur.vertex << "," << cur.distance << ") ";
        }
        std::cout << "\n";
    }
}

void propagate(const std::vector<dgraph_v_of_v<two_hop_weight_type>> &G, int k, int u)
{
    mtx.lock();
    int current_tid = thread_id.front();    
    //std::cout << "thread id:" << current_tid << " k:"
    //          <<k << " u:" << u <<endl;
    //cout << endl;
    thread_id.pop();
    mtx.unlock();

    for (auto &cur : L_PSL[k][u])
    {
        int v = cur.vertex;
        two_hop_weight_type d = cur.distance;
        if (dirt[current_tid][v])
        {
            dirt[current_tid][v] = 0;
            dmin[current_tid][v] = d;
        }
        else
        {
            dmin[current_tid][v] = std::min(dmin[current_tid][v], d);
        }
    }

    for (auto label : G[k].OUTs[u])
    {
        int v = label.first;
        two_hop_weight_type dd = label.second;

        for (int i = pos_595[k][v]; i < pos_595[k][v] + increment[k][v]; i++)
        //for (int i = 0; i < pos_595[k][v] + increment[k][v]; i++)
        {
            int x = L_PSL[k][v][i].vertex;
            if (x > u)
            { // higher rank means lower value of `rank`
                continue;
            }
            two_hop_weight_type d_temp = dd + L_PSL[k][v][i].distance;

            bool flag = false;
            for (int j = 0; j < L_PSL[k ^ 1][x].size(); j++)
            {
                int y = L_PSL[k ^ 1][x][j].vertex;
                if (!dirt[current_tid][y] && ((dmin[current_tid][y] + L_PSL[k ^ 1][x][j].distance) <= d_temp))
                {
                    flag = true;
                    break;
                }
            }
            if (flag)
            {
                continue;
            }

            // Line 17
            L_PSL_temp[k][u].push_back({x, d_temp});
        }
    }

    for (auto &cur : L_PSL[k][u])
    {
        int v = cur.vertex;
        dirt[current_tid][v] = 1;
    }

    mtx.lock();
    thread_id.push(current_tid);
    mtx.unlock();
}

void propagate_v1(const std::vector<dgraph_v_of_v<two_hop_weight_type>> &G, int k, int u)
{
    mtx.lock();
    int current_tid = thread_id.front();
    // std::cout << "thread id:" << current_tid << " k:"
    //           <<k << " u:" << u <<endl;
    // cout << endl;
    thread_id.pop();
    mtx.unlock();

    int N = G[0].INs.size();

    vector<two_hop_weight_type> dmin(N, std::numeric_limits<two_hop_weight_type>::max());
    vector<int> dirt(N, std::numeric_limits<two_hop_weight_type>::max());


    for (auto &cur : L_PSL[k][u])
    {
        int v = cur.vertex;
        two_hop_weight_type d = cur.distance;
        if (dirt[v])
        {
            dirt[v] = 0;
            dmin[v] = d;
        }
        else
        {
            dmin[v] = std::min(dmin[v], d);
        }
    }

    for (auto label : G[k].OUTs[u])
    {
        int v = label.first;
        two_hop_weight_type dd = label.second;

        for (int i = pos_595[k][v]; i < pos_595[k][v] + increment[k][v]; i++)
        // for (int i = 0; i < pos_595[k][v] + increment[k][v]; i++)
        {
            int x = L_PSL[k][v][i].vertex;
            if (x > u)
            { // higher rank means lower value of `rank`
                continue;
            }
            two_hop_weight_type d_temp = dd + L_PSL[k][v][i].distance;

            bool flag = false;
            for (int j = 0; j < L_PSL[k ^ 1][x].size(); j++)
            {
                int y = L_PSL[k ^ 1][x][j].vertex;
                if (!dirt[y] && ((dmin[y] + L_PSL[k ^ 1][x][j].distance) <= d_temp))
                {
                    flag = true;
                    break;
                }
            }
            if (flag)
            {
                continue;
            }

            // Line 17
            L_PSL_temp[k][u].push_back({x, d_temp});
        }
    }

    /*
    for (auto &cur : L_PSL[k][u])
    {
        int v = cur.vertex;
        dirt[v] = 1;
    }
    */

    mtx.lock();
    thread_id.push(current_tid);
    mtx.unlock();
}

void append(int k, int u)
{   
    pos_2_595[k][u] = pos_595[k][u];
    pos_595[k][u] += increment[k][u];
    increment[k][u] = L_PSL_temp[k][u].size();

    L_PSL[k][u].insert(L_PSL[k][u].end(), L_PSL_temp[k][u].begin(), L_PSL_temp[k][u].end());
    vector<two_hop_label>(L_PSL[k][u]).swap(L_PSL[k][u]);
}

void shrink(int k, int u, int N)
{
    mtx.lock();
    int current_tid = thread_id.front();
    thread_id.pop();
    mtx.unlock();

    vector<int> in_psl(N,0);
    vector<two_hop_weight_type> min_distance_in_hub(N);

    for (auto &cur : L_PSL[k][u])
    {
        int v = cur.vertex;
        two_hop_weight_type d = cur.distance;
        if (!in_psl[v])
        {
            in_psl[v] = 1;
            min_distance_in_hub[v] = d;
        }
        else if (d < min_distance_in_hub[v])
        {
            min_distance_in_hub[v] = d;
        }
    }

    for (int i = 0; i < N; i++)
    {
        if (!in_psl[i])
        {
            continue;
        }
        L_PSL[k][u].clear();
        L_PSL_temp[k][u].push_back({i, min_distance_in_hub[i]});
    }

    mtx.lock();
    thread_id.push(current_tid);
    mtx.unlock();
}

void dgraph_PSL_v3(dgraph_v_of_v<two_hop_weight_type> &input_graph, int N, int num_of_threads, dgraph_case_info_v1 &case_info)
{
    std::vector<dgraph_v_of_v<two_hop_weight_type>> G(2);
    G[0] = input_graph;
    G[1] = input_graph.reverse_graph();

    // int N = ideal_dgraph.INs.size();

    // prepare for thread pool
    dmin.resize(num_of_threads);
    dirt.resize(num_of_threads);
    for (int i = 0; i < num_of_threads; i++)
    {
        dmin[i].resize(N);
        dirt[i].resize(N);
        for (int j = 0; j < N; j++)
            dirt[i][j] = true;
        thread_id.push(i);
    }

    for (int k = 0; k < 2; k++)
    {
        L_PSL[k].resize(N);
        L_PSL_temp[k].resize(N);
        // endpos1[k].resize(N);
        pos_595[k].resize(N);
        pos_2_595[k].resize(N);
        increment[k].resize(N);
        for (int i = 0; i < N; i++)
        {
            L_PSL[k][i].push_back({i, 0});
            pos_595[k][i] = 0;
            pos_2_595[k][i] = 0;
            increment[k][i] = 1;
        }
    }

    ThreadPool pool(num_of_threads);
    std::vector<std::future<int>> results;

    while (true)
    {
        // Line 4-20
        bool is_empty[2] = {true, true};
        for (int k = 0; k < 2; k++)
        {
            for (int u = 0; u < N; u++)
            {

                results.emplace_back(pool.enqueue([u,k ,&G] {
                    propagate_v1(G, k, u);                   
                    return 1;
                }));
            }
            for (auto &&result : results)
            {
                result.get();
            }
            results.clear();

            for (int u = 0; u < N; u++)
            {
                results.emplace_back(pool.enqueue([k,u] {
                    append(k, u);                   
                    return 1;
                }));
            }
            for (auto &&result : results)
            {
                result.get();
            }
            results.clear();
            
            for (int u = 0; u < N; u++)
            {
                if (!L_PSL_temp[k][u].empty())
                {
                    is_empty[k] = false;
                }
                // endpos1[k][u] = L_PSL[k][u].size() - L_PSL_temp[k][u].size();
                vector<two_hop_label>().swap(L_PSL_temp[k][u]);
            }
        }

        if (is_empty[0] && is_empty[1])
        {
            break;
        }
    } // while loop

    for (int k = 0; k < 2; k++)
    {
        for (int u = 0; u < N; u++)
        {
            /* save the smallest label for each hub of u in T */
            int u_label_size = L_PSL[k][u].size();
            for (int i = 0; i < u_label_size; i++)
            {
                int w = L_PSL[k][u][i].vertex;
                double dis = L_PSL[k][u][i].distance;
                if (dirt[0][w]) // the same hub may have redundancy, record the shortest distance
                {
                    dirt[0][w] = false;
                    dmin[0][w] = dis;
                }
                else if (dis < dmin[0][w])
                {
                    dmin[0][w] = dis;
                }
            }

            for (int i = 0; i < N; i++)
            {
                if (!dirt[0][i])
                {
                    dirt[0][i] = true;
                    two_hop_label xx;
                    xx.vertex = i;
                    xx.distance = dmin[0][i];
                    L_PSL_temp[k][u].push_back(xx);
                }
            }

            vector<two_hop_label>().swap(L_PSL[k][u]); // do not have two L in RAM simultaneously
        }
    }

    L_temp_in = L_PSL_temp[1];
    L_temp_out = L_PSL_temp[0];

    if (case_info.use_canonical_repair)
    {
        case_info.label_size_before_canonical_repair = compute_label_bit_size(L_temp_in, L_temp_out);
        canonical_repair_multi_threads(num_of_threads, &case_info.L_in, &case_info.L_out);  
        case_info.label_size_after_canonical_repair = compute_label_bit_size(case_info.L_in, case_info.L_out);
    }
    else {
        case_info.L_in = L_temp_in;
        case_info.L_out = L_temp_out;
    }

    vector<vector<two_hop_label>>().swap(L_temp_in);
    vector<vector<two_hop_label>>().swap(L_temp_out);

    dgraph_clear_global_values_PSL();
}