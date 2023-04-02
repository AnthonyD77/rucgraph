#pragma once
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <tool_functions/ThreadPool.h>
#include <shared_mutex>
#include <chrono>
#include <boost/heap/fibonacci_heap.hpp>
#include <dgraph_v_of_v/dgraph_v_of_v.h>
#include <build_in_progress/HL/dgraph/dgraph_two_hop_label.h>

/*
    to include this file, please copy
    #include <dgraph_v_of_v/dgraph_dijkstra.h>
*/

void print_L(vector<vector<two_hop_label>> &L_in, vector<vector<two_hop_label>> &L_out)
{
    cout << "in label" << std::endl;

    int size1 = L_in.size();
    for (int i = 0; i < size1; i++)
    {
        cout << i << ":\t";
        int size2 = L_in[i].size();
        for (int j = 0; j < size2; j++)
        {
            cout << "(" << L_in[i][j].vertex << "," << L_in[i][j].distance << ") ";
        }
        cout << std::endl;
    }

    cout << std::endl << "out label" << std::endl;
    
    size1 = L_out.size();
    for (int i = 0; i < size1; i++)
    {
        cout << i << ":\t";
        int size2 = L_out[i].size();
        for (int j = 0; j < size2; j++)
        {
            cout << "(" << L_out[i][j].vertex << "," << L_out[i][j].distance << ") ";
        }
        cout << std::endl;
    }
    cout << std::endl;
    cout << std::endl;
}


/* struct used for dijkstra extra_min */
struct node_for_dij
{
  public:
    int vertex;
    two_hop_weight_type priority_value;
};

bool operator<(node_for_dij const &x, node_for_dij const &y)
{   
    return x.priority_value > y.priority_value; // < is the max-heap; > is the min heap
}
typedef typename boost::heap::fibonacci_heap<node_for_dij>::handle_type dgraph_heap_pointer;

/* single source dijkstra */
void dgraph_dijkstra(int v_k, int N, string in_out)
{
    int debug = 0;

    mtx_595[max_N_ID_for_mtx_595 - 1].lock();
    int thread_id = Qid_595.front();
    Qid_595.pop();
    mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

    queue<int> dij_dist_changed_vertices;
    vector<dgraph_heap_pointer> Q_pointer(N);
    node_for_dij node;
    boost::heap::fibonacci_heap<node_for_dij> Q;

    two_hop_label xx;
    node.vertex = v_k;
    node.priority_value = 0;

    Q_pointer[v_k] = Q.push(node);
    dij_dist[thread_id][v_k] = 0;
    dij_dist_changed_vertices.push(v_k);

    long long int new_label_num = 0;

    while (Q.size() > 0)
    {
        /* find the shortest node */
        node = Q.top();
        Q.pop();
        int u = node.vertex;
        two_hop_weight_type dist_u = node.priority_value;

        if (debug)
        {
            cout << "source: " << u << endl;
        }

        /* add new label - in part */
        if (in_out == "in")
        {
            xx.vertex = v_k;
            xx.distance = dist_u;
            // xx.parent_vertex = u_parent;
            mtx_595[u].lock();
            L_temp_in[u].push_back(xx);
            mtx_595[u].unlock();
        }
        /* add new label - out part */
        else
        {
            xx.vertex = v_k;
            xx.distance = dist_u;
            // xx.parent_vertex = u_parent;
            mtx_595[u].lock();
            L_temp_out[u].push_back(xx);
            mtx_595[u].unlock();
        }
        new_label_num++;

        /* update part of dijstra */
        int u_adj_size = ideal_dgraph.OUTs[u].size();
        for (int i = 0; i < u_adj_size; i++)
        {
            /* adj_v is a neighbor of u */
            int adj_v = ideal_dgraph.OUTs[u][i].first;
            two_hop_weight_type new_dist = ideal_dgraph.OUTs[u][i].second + dist_u;

            if (debug)
            {
                cout << "neighbor: " << adj_v << "\tdistance: " << new_dist << endl;
            }

            if (dij_dist[thread_id][adj_v] == std::numeric_limits<two_hop_weight_type>::max())
            {
                node.vertex = adj_v;
                // node.parent_vertex = u;
                node.priority_value = new_dist;
                Q_pointer[adj_v] = Q.push(node);
                dij_dist[thread_id][adj_v] = node.priority_value;
                dij_dist_changed_vertices.push(adj_v);
            }
            else // v is already in the dij_dist, only need to check if update the dij_dist
            {
                if (dij_dist[thread_id][adj_v] > new_dist)
                {
                    node.vertex = adj_v;
                    // node.parent_vertex = u;
                    node.priority_value = new_dist;
                    Q.update(Q_pointer[adj_v], node);
                    dij_dist[thread_id][adj_v] = node.priority_value;
                }
            }
        }
    }

    while (dij_dist_changed_vertices.size() > 0)
    {
        dij_dist[thread_id][dij_dist_changed_vertices.front()] = std::numeric_limits<two_hop_weight_type>::max();
        dij_dist_changed_vertices.pop();
    }

    mtx_595[v_k].lock();
    vector<two_hop_label>(L_temp_in[v_k]).swap(L_temp_in[v_k]);
    mtx_595[v_k].unlock();

    mtx_595[max_N_ID_for_mtx_595 - 1].lock();
    Qid_595.push(thread_id);
    labal_size_595 = labal_size_595 + new_label_num;
    mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

    // if (labal_size_595 > max_labal_size_595)
    // {
    //     throw reach_limit_error_string_MB;
    //     // after catching error, must call terminate_procedures_595(), otherwise this
    //     // PLL cannot be reused
    // }

    // if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() -
    // begin_time_595)
    //         .count() > max_run_time_nanoseconds_595)
    // {
    //     throw reach_limit_error_string_time;
    //     // after catching error, must call terminate_procedures_595(), otherwise
    //     // this PLL cannot be reused
    // }
}

void dgraph_pruned_in_out_dijstra(int v_k, int N)
{
    int debug = 0;

    mtx_595[max_N_ID_for_mtx_595 - 1].lock();
    int thread_id = Qid_595.front();
    Qid_595.pop();
    // cout << "thread id:" << thread_id << " v_k:"<<v_k << endl;
    mtx_595[max_N_ID_for_mtx_595 - 1].unlock();
    
    //change to vector
    queue<int> dij_dist_changed_vertices, L_vk_changed_vertices;

    /* get the L[vk] first so qeury can be faster, since L[vk] will not change in this loop */
    mtx_595[v_k].lock();
    int L_out_vk_size = L_temp_out[v_k].size();
    for (int i = 0; i < L_out_vk_size; i++)
    {
        int L_vk_i_vertex = L_temp_out[v_k][i].vertex;
        L_vk_tmp[thread_id][L_vk_i_vertex] = L_temp_out[v_k][i].distance;
        L_vk_changed_vertices.push(L_vk_i_vertex);
    }
    mtx_595[v_k].unlock();

    vector<dgraph_heap_pointer> Q_pointer(N);
    node_for_dij node;
    boost::heap::fibonacci_heap<node_for_dij> Q;

    two_hop_label xx;
    node.vertex = v_k;
    node.priority_value = 0;

    Q_pointer[v_k] = Q.push(node);
    dij_dist[thread_id][v_k] = 0;
    dij_dist_changed_vertices.push(v_k);

    long long int new_label_num = 0;

    while (Q.size() > 0)
    {
        /* find the shortest node */
        node = Q.top();
        Q.pop();
        int u = node.vertex;
        two_hop_weight_type dist_u = node.priority_value;
        

        if (debug)
        {
            cout << "source: " << u << endl;
        }

        if (v_k > u)
            continue;

        /* pruned condition -- key part of PLL */
        /* query: vk -> mid_node -> u, so we refer L_out[vk] and L_in[u] */
        two_hop_weight_type query_vk_u = std::numeric_limits<two_hop_weight_type>::max();

        mtx_595[u].lock();
        int L_in_u_size = L_temp_in[u].size();
        mtx_595[u].unlock();
        for (int i = 0; i < L_in_u_size; i++)
        {   
            mtx_595[u].lock();
            two_hop_weight_type dis = L_temp_in[u][i].distance + L_vk_tmp[thread_id][L_temp_in[u][i].vertex];
            mtx_595[u].unlock();
            if (query_vk_u > dis)
            {
                query_vk_u = dis;
            }
        }
        
        if (query_vk_u <= dist_u)
            continue;

        /* add new label - in part */
        xx.vertex = v_k;
        xx.distance = dist_u;

        mtx_595[u].lock();
        L_temp_in[u].push_back(xx);
        mtx_595[u].unlock();

        new_label_num++;

        /* update part of dijstra */
        int u_adj_size;
        u_adj_size = ideal_dgraph.OUTs[u].size();

        for (int i = 0; i < u_adj_size; i++)
        {
            /* adj_v is a neighbor of u */
            int adj_v;
            two_hop_weight_type new_dist;

            adj_v = ideal_dgraph.OUTs[u][i].first;
            new_dist = ideal_dgraph.OUTs[u][i].second + dist_u;

            if (debug)
            {
                cout << "neighbor: " << adj_v << "\tdistance: " << new_dist << endl;
            }

            if (dij_dist[thread_id][adj_v] == std::numeric_limits<two_hop_weight_type>::max())
            {
                node.vertex = adj_v;
                node.priority_value = new_dist;
               
                Q_pointer[adj_v] = Q.push(node);
                dij_dist[thread_id][adj_v] = node.priority_value;
                dij_dist_changed_vertices.push(adj_v);
            }
            else // v is already in the dij_dist, only need to check if update the dij_dist
            {
                if (dij_dist[thread_id][adj_v] > new_dist)
                {
                    node.vertex = adj_v;
                    node.priority_value = new_dist;
                   
                    Q.update(Q_pointer[adj_v], node);
                    dij_dist[thread_id][adj_v] = node.priority_value;
                }
            }
        }
    }

    /* reset dij_dist and L_vk_tmp */
    while (dij_dist_changed_vertices.size() > 0)
    {
        dij_dist[thread_id][dij_dist_changed_vertices.front()] = std::numeric_limits<two_hop_weight_type>::max();
        dij_dist_changed_vertices.pop();
    }
    while (L_vk_changed_vertices.size() > 0)
    {
        L_vk_tmp[thread_id][L_vk_changed_vertices.front()] = std::numeric_limits<two_hop_weight_type>::max();
        L_vk_changed_vertices.pop();
    }

    /* clear excess space */
    mtx_595[v_k].lock();
    vector<two_hop_label>(L_temp_in[v_k]).swap(L_temp_in[v_k]);
    mtx_595[v_k].unlock();

    mtx_595[v_k].lock();
    int L_in_vk_size = L_temp_in[v_k].size();
    for (int i = 0; i < L_in_vk_size; i++)
    {
        int L_vk_i_vertex = L_temp_in[v_k][i].vertex;
        L_vk_tmp[thread_id][L_vk_i_vertex] = L_temp_in[v_k][i].distance;
        L_vk_changed_vertices.push(L_vk_i_vertex);
    }

    mtx_595[v_k].unlock();

    //vector<dgraph_heap_pointer> Q_pointer(N);
    //boost::heap::fibonacci_heap<node_for_dij> Q;

    node.vertex = v_k;
    node.priority_value = 0;
   
    Q_pointer[v_k] = Q.push(node);
    dij_dist[thread_id][v_k] = 0;
    dij_dist_changed_vertices.push(v_k);

    while (Q.size() > 0)
    {
        /* find the shortest node */
        node = Q.top();
        Q.pop();
        int u = node.vertex;
        two_hop_weight_type dist_u = node.priority_value;
        

        if (debug)
        {
            cout << "source: " << u << endl;
        }

        if (v_k > u)
            continue;

        /* pruned condition -- key part of PLL */
        /* query: vk -> mid_node -> u, so we refer L_out[vk] and L_in[u] */
        two_hop_weight_type query_vk_u = std::numeric_limits<two_hop_weight_type>::max();
        
        mtx_595[u].lock();
        int L_out_u_size = L_temp_out[u].size();
        mtx_595[u].unlock();
        for (int i = 0; i < L_out_u_size; i++)
        {   
            mtx_595[u].lock();
            two_hop_weight_type dis = L_temp_out[u][i].distance + L_vk_tmp[thread_id][L_temp_out[u][i].vertex];
            mtx_595[u].unlock();
            if (query_vk_u > dis)
            {
                query_vk_u = dis;
            }
        }

        if (query_vk_u <= dist_u)
            continue;

        xx.vertex = v_k;
        xx.distance = dist_u;
        mtx_595[u].lock();
        L_temp_out[u].push_back(xx);
        mtx_595[u].unlock();

        new_label_num++;

        /* update part of dijstra */
        int u_adj_size;
        u_adj_size = revse_dgraph.OUTs[u].size();

        for (int i = 0; i < u_adj_size; i++)
        {
            /* adj_v is a neighbor of u */
            int adj_v;
            two_hop_weight_type new_dist;

            adj_v = revse_dgraph.OUTs[u][i].first;
            new_dist = revse_dgraph.OUTs[u][i].second + dist_u;

            if (debug)
            {
                cout << "neighbor: " << adj_v << "\tdistance: " << new_dist << endl;
            }

            if (dij_dist[thread_id][adj_v] == std::numeric_limits<two_hop_weight_type>::max())
            {
                node.vertex = adj_v;
                node.priority_value = new_dist;
                
                Q_pointer[adj_v] = Q.push(node);
                dij_dist[thread_id][adj_v] = node.priority_value;
                dij_dist_changed_vertices.push(adj_v);
            }
            else // v is already in the dij_dist, only need to check if update the dij_dist
            {
                if (dij_dist[thread_id][adj_v] > new_dist)
                {
                    node.vertex = adj_v;
                    node.priority_value = new_dist;
                    
                    Q.update(Q_pointer[adj_v], node);
                    dij_dist[thread_id][adj_v] = node.priority_value;
                }
            }
        }
    }

    /* reset dij_dist and L_vk_tmp */
    while (dij_dist_changed_vertices.size() > 0)
    {
        dij_dist[thread_id][dij_dist_changed_vertices.front()] = std::numeric_limits<two_hop_weight_type>::max();
        dij_dist_changed_vertices.pop();
    }
    while (L_vk_changed_vertices.size() > 0)
    {
        L_vk_tmp[thread_id][L_vk_changed_vertices.front()] = std::numeric_limits<two_hop_weight_type>::max();
        L_vk_changed_vertices.pop();
    }

    /* clear excess space */
    mtx_595[v_k].lock();
    vector<two_hop_label>(L_temp_out[v_k]).swap(L_temp_out[v_k]);
    mtx_595[v_k].unlock();

    /* recycle thread */
    mtx_595[max_N_ID_for_mtx_595 - 1].lock();
    Qid_595.push(thread_id);
    labal_size_595 = labal_size_595 + new_label_num;
    mtx_595[max_N_ID_for_mtx_595 - 1].unlock();
}

/* pruned dijkstra used for PLL */
void dgraph_pruned_dijkstra(int v_k, int N, int in_out)
{
    int debug = 0;

    //cout << "enter dgraph_dijkstra" << endl;

    mtx_595[max_N_ID_for_mtx_595 - 1].lock();
    int thread_id = Qid_595.front();
    Qid_595.pop();
    //cout << "this is thread:" << thread_id << " v_k=" << v_k << " in_out=" << in_out<<endl;
    mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

    vector<int> dij_dist_changed_vertices(N,0), L_vk_changed_vertices(N,0);

    /* get the L[vk] first so qeury can be faster, since L[vk] will not change in this loop */
    
    mtx_595[v_k].lock();
    if (in_out == 0)
    {   
        //mtx_595[v_k].lock();
        int L_out_vk_size = L_temp_out[v_k].size();
        for (int i = 0; i < L_out_vk_size; i++)
        {   
            //mtx_595[v_k].lock();
            int L_vk_i_vertex = L_temp_out[v_k][i].vertex;
            L_vk_tmp[thread_id][L_vk_i_vertex] = L_temp_out[v_k][i].distance;
            //mtx_595[v_k].unlock();
            L_vk_changed_vertices[L_vk_i_vertex] = 1;
        }
    }
    else // reverse case, we refer to L_in[vk]
    {   
        //mtx_595[v_k].lock();
        int L_in_vk_size = L_temp_in[v_k].size();        
        for (int i = 0; i < L_in_vk_size; i++)
        {   
            //mtx_595[v_k].lock();
            int L_vk_i_vertex = L_temp_in[v_k][i].vertex;
            L_vk_tmp[thread_id][L_vk_i_vertex] = L_temp_in[v_k][i].distance;
            //mtx_595[v_k].unlock();
            L_vk_changed_vertices[L_vk_i_vertex] = 1;
        }
    }
    mtx_595[v_k].unlock();

    vector<dgraph_heap_pointer> Q_pointer(N);
    node_for_dij node;
    boost::heap::fibonacci_heap<node_for_dij> Q;

    two_hop_label xx;
    node.vertex = v_k;
    node.priority_value = 0;
    

    Q_pointer[v_k] = Q.push(node);
    //mtx_595[v_k].lock();
    dij_dist[thread_id][v_k] = 0;
    //mtx_595[v_k].unlock();
    dij_dist_changed_vertices[v_k] = 1;

    long long int new_label_num = 0;

    while (Q.size() > 0)
    {
        /* find the shortest node */
        node = Q.top();
        Q.pop();
        int u = node.vertex;
        two_hop_weight_type dist_u = node.priority_value;
        

        if (debug)
        {
            cout << "source: " << u << endl;
        }

        if (v_k <= u)
        {
            /* pruned condition -- key part of PLL */
            /* query: vk -> mid_node -> u, so we refer L_out[vk] and L_in[u] */
            two_hop_weight_type query_vk_u = std::numeric_limits<two_hop_weight_type>::max();
            //mtx_595[u].lock();
            if (in_out == 0)
            {   
                mtx_595[u].lock();
                int L_in_u_size = L_temp_in[u].size();
                mtx_595[u].unlock();
                for (int i = 0; i < L_in_u_size; i++)
                {
                    mtx_595[u].lock();
                    two_hop_weight_type dis = L_temp_in[u][i].distance + L_vk_tmp[thread_id][L_temp_in[u][i].vertex];
                    mtx_595[u].unlock();
                    if (query_vk_u > dis)
                    {
                        query_vk_u = dis;
                    }
                }
            }
            else // reverse case
            {   
                mtx_595[u].lock();
                int L_out_u_size = L_temp_out[u].size();
                mtx_595[u].unlock();
                for (int i = 0; i < L_out_u_size; i++)
                {
                    mtx_595[u].lock();
                    two_hop_weight_type dis = L_temp_out[u][i].distance + L_vk_tmp[thread_id][L_temp_out[u][i].vertex];
                    mtx_595[u].unlock();
                    if (query_vk_u > dis)
                    {
                        query_vk_u = dis;
                    }
                }
            }
            //mtx_595[u].unlock();
            if (query_vk_u <= dist_u)
                continue;

            /* add new label - in part */
            if (in_out == 0)
            {
                xx.vertex = v_k;
                xx.distance = dist_u;

                mtx_595[u].lock();
                L_temp_in[u].push_back(xx);
                mtx_595[u].unlock();
            }
            /* add new label - out part */
            else
            {
                xx.vertex = v_k;
                xx.distance = dist_u;

                mtx_595[u].lock();
                L_temp_out[u].push_back(xx);
                mtx_595[u].unlock();
            }
            new_label_num++;

            /* update part of dijstra */
            int u_adj_size;
            if (in_out == 0)
                u_adj_size = ideal_dgraph.OUTs[u].size();
            else
                u_adj_size = revse_dgraph.OUTs[u].size();

            for (int i = 0; i < u_adj_size; i++)
            {
                /* adj_v is a neighbor of u */
                int adj_v;
                two_hop_weight_type new_dist;
                if (in_out == 0)
                {
                    adj_v = ideal_dgraph.OUTs[u][i].first;
                    new_dist = ideal_dgraph.OUTs[u][i].second + dist_u;
                }
                else
                {
                    adj_v = revse_dgraph.OUTs[u][i].first;
                    new_dist = revse_dgraph.OUTs[u][i].second + dist_u;
                }

                if (debug)
                {
                    cout << "neighbor: " << adj_v << "\tdistance: " << new_dist << endl;
                }

                //if (dij_dist[thread_id][adj_v] == std::numeric_limits<two_hop_weight_type>::max())
                if (!dij_dist_changed_vertices[adj_v])
                {
                    node.vertex = adj_v;
                    node.priority_value = new_dist;

                    Q_pointer[adj_v] = Q.push(node);
                    dij_dist[thread_id][adj_v] = new_dist;
                    dij_dist_changed_vertices[adj_v] = 1;
                }        
                
                else // v is already in the dij_dist, only need to check if update the dij_dist
                {
                    
                    if (dij_dist[thread_id][adj_v] > new_dist)
                    {
                        node.vertex = adj_v;
                        node.priority_value = new_dist;
                        
                        Q.update(Q_pointer[adj_v], node);
                        dij_dist[thread_id][adj_v] = new_dist;
                    }
                }
                
            }
        }        
    }

    //cout << "end while" << endl;

    for (int i=0;i<N;i++)
    {
        dij_dist[thread_id][i] = std::numeric_limits<two_hop_weight_type>::max();
        L_vk_tmp[thread_id][i] = std::numeric_limits<two_hop_weight_type>::max();
    }

    /*
    mtx.lock();
    cout << "vk:" << v_k << " in_out:" << in_out<<endl;
    print_L(L_temp_in,L_temp_out);
    
    //if (in_out ==0 )
    //    vector<two_hop_label>(L_temp_in[v_k]).swap(L_temp_in[v_k]);
    //else
    //    vector<two_hop_label>(L_temp_out[v_k]).swap(L_temp_out[v_k]);
    mtx.unlock();
    */

    /* recycle thread */
    mtx_595[max_N_ID_for_mtx_595 - 1].lock();
    Qid_595.push(thread_id);
    labal_size_595 = labal_size_595 + new_label_num;
    mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

    // if (labal_size_595 > max_labal_size_595)
    // {
    //     throw reach_limit_error_string_MB;
    //     // after catching error, must call terminate_procedures_595(), otherwise this
    //     // PLL cannot be reused
    // }

    // if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() -
    // begin_time_595)
    //         .count() > max_run_time_nanoseconds_595)
    // {
    //     throw reach_limit_error_string_time;
    //     // after catching error, must call terminate_procedures_595(), otherwise
    //     // this PLL cannot be reused
    // }
}

void dgraph_pruned_dijkstra_v1(int v_k, int N, int in_out)
{
    int debug = 0;

    // cout << "enter dgraph_dijkstra" << endl;

    mtx_595[max_N_ID_for_mtx_595 - 1].lock();
    int thread_id = Qid_595.front();
    Qid_595.pop();
    // cout << "this is thread:" << thread_id << " v_k=" << v_k << " in_out=" << in_out<<endl;
    mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

    vector<int> dij_dist_changed_vertices(N, 0), L_vk_changed_vertices(N, 0);
    vector<two_hop_weight_type> L_vk_temp_in_function(N, std::numeric_limits<two_hop_weight_type>::max());
    vector<two_hop_weight_type> dist(N, std::numeric_limits<two_hop_weight_type>::max());

    /* get the L[vk] first so qeury can be faster, since L[vk] will not change in this loop */

    mtx_595[v_k].lock();
    if (in_out == 0)
    {
        // mtx_595[v_k].lock();
        int L_out_vk_size = L_temp_out[v_k].size();
        for (int i = 0; i < L_out_vk_size; i++)
        {
            // mtx_595[v_k].lock();
            int L_vk_i_vertex = L_temp_out[v_k][i].vertex;
            L_vk_temp_in_function[L_vk_i_vertex] = L_temp_out[v_k][i].distance;
            // mtx_595[v_k].unlock();
            L_vk_changed_vertices[L_vk_i_vertex] = 1;
        }
    }
    else // reverse case, we refer to L_in[vk]
    {
        // mtx_595[v_k].lock();
        int L_in_vk_size = L_temp_in[v_k].size();
        for (int i = 0; i < L_in_vk_size; i++)
        {
            // mtx_595[v_k].lock();
            int L_vk_i_vertex = L_temp_in[v_k][i].vertex;
            L_vk_temp_in_function[L_vk_i_vertex] = L_temp_in[v_k][i].distance;
            // mtx_595[v_k].unlock();
            L_vk_changed_vertices[L_vk_i_vertex] = 1;
        }
    }
    mtx_595[v_k].unlock();

    vector<dgraph_heap_pointer> Q_pointer(N);
    node_for_dij node;
    boost::heap::fibonacci_heap<node_for_dij> Q;

    two_hop_label xx;
    node.vertex = v_k;
    node.priority_value = 0;

    Q_pointer[v_k] = Q.push(node);
    // mtx_595[v_k].lock();
    dist[v_k] = 0;
    // mtx_595[v_k].unlock();
    dij_dist_changed_vertices[v_k] = 1;

    long long int new_label_num = 0;

    while (Q.size() > 0)
    {
        /* find the shortest node */
        node = Q.top();
        Q.pop();
        int u = node.vertex;
        two_hop_weight_type dist_u = node.priority_value;

        if (debug)
        {
            cout << "source: " << u << endl;
        }

        if (v_k <= u)
        {
            /* pruned condition -- key part of PLL */
            /* query: vk -> mid_node -> u, so we refer L_out[vk] and L_in[u] */
            two_hop_weight_type query_vk_u = std::numeric_limits<two_hop_weight_type>::max();
            mtx_595[u].lock();
            if (in_out == 0)
            {
                //mtx_595[u].lock();
                int L_in_u_size = L_temp_in[u].size();
                //mtx_595[u].unlock();
                for (int i = 0; i < L_in_u_size; i++)
                {
                    //mtx_595[u].lock();
                    two_hop_weight_type dis =
                        L_temp_in[u][i].distance + L_vk_temp_in_function[L_temp_in[u][i].vertex];
                    //mtx_595[u].unlock();
                    if (query_vk_u > dis)
                    {
                        query_vk_u = dis;
                    }
                }
            }
            else // reverse case
            {
                //mtx_595[u].lock();
                int L_out_u_size = L_temp_out[u].size();
                //mtx_595[u].unlock();
                for (int i = 0; i < L_out_u_size; i++)
                {
                    //mtx_595[u].lock();
                    two_hop_weight_type dis =
                        L_temp_out[u][i].distance + L_vk_temp_in_function[L_temp_out[u][i].vertex];
                    //mtx_595[u].unlock();
                    if (query_vk_u > dis)
                    {
                        query_vk_u = dis;
                    }
                }
            }
            mtx_595[u].unlock();
            if (query_vk_u <= dist_u)
                continue;

            /* add new label - in part */
            if (in_out == 0)
            {
                xx.vertex = v_k;
                xx.distance = dist_u;

                mtx_595[u].lock();
                L_temp_in[u].push_back(xx);
                mtx_595[u].unlock();
            }
            /* add new label - out part */
            else
            {
                xx.vertex = v_k;
                xx.distance = dist_u;

                mtx_595[u].lock();
                L_temp_out[u].push_back(xx);
                mtx_595[u].unlock();
            }
            new_label_num++;

            /* update part of dijstra */
            int u_adj_size;
            if (in_out == 0)
                u_adj_size = ideal_dgraph.OUTs[u].size();
            else
                u_adj_size = revse_dgraph.OUTs[u].size();

            for (int i = 0; i < u_adj_size; i++)
            {
                /* adj_v is a neighbor of u */
                int adj_v;
                two_hop_weight_type new_dist;
                if (in_out == 0)
                {
                    adj_v = ideal_dgraph.OUTs[u][i].first;
                    new_dist = ideal_dgraph.OUTs[u][i].second + dist_u;
                }
                else
                {
                    adj_v = revse_dgraph.OUTs[u][i].first;
                    new_dist = revse_dgraph.OUTs[u][i].second + dist_u;
                }

                if (debug)
                {
                    cout << "neighbor: " << adj_v << "\tdistance: " << new_dist << endl;
                }

                // if (dij_dist[thread_id][adj_v] == std::numeric_limits<two_hop_weight_type>::max())
                if (!dij_dist_changed_vertices[adj_v])
                {
                    node.vertex = adj_v;
                    node.priority_value = new_dist;

                    Q_pointer[adj_v] = Q.push(node);
                    dist[adj_v] = new_dist;
                    dij_dist_changed_vertices[adj_v] = 1;
                }               
                else // v is already in the dij_dist, only need to check if update the dij_dist
                {

                    if (dist[adj_v] > new_dist)
                    {
                        node.vertex = adj_v;
                        node.priority_value = new_dist;

                        Q.update(Q_pointer[adj_v], node);
                        dist[adj_v] = new_dist;
                    }
                }
                
            }
        }
    }

    // cout << "end while" << endl;

    for (int i = 0; i < N; i++)
    {
        dij_dist[thread_id][i] = std::numeric_limits<two_hop_weight_type>::max();
        L_vk_tmp[thread_id][i] = std::numeric_limits<two_hop_weight_type>::max();
    }

    /*
    mtx.lock();
    cout << "vk:" << v_k << " in_out:" << in_out<<endl;
    print_L(L_temp_in,L_temp_out);

    //if (in_out ==0 )
    //    vector<two_hop_label>(L_temp_in[v_k]).swap(L_temp_in[v_k]);
    //else
    //    vector<two_hop_label>(L_temp_out[v_k]).swap(L_temp_out[v_k]);
    mtx.unlock();
    */

    /* recycle thread */
    mtx_595[max_N_ID_for_mtx_595 - 1].lock();
    Qid_595.push(thread_id);
    labal_size_595 = labal_size_595 + new_label_num;
    mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

    // if (labal_size_595 > max_labal_size_595)
    // {
    //     throw reach_limit_error_string_MB;
    //     // after catching error, must call terminate_procedures_595(), otherwise this
    //     // PLL cannot be reused
    // }

    // if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() -
    // begin_time_595)
    //         .count() > max_run_time_nanoseconds_595)
    // {
    //     throw reach_limit_error_string_time;
    //     // after catching error, must call terminate_procedures_595(), otherwise
    //     // this PLL cannot be reused
    // }
}
