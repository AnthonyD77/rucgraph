#pragma once

#include <unordered_map>
#include <boost/heap/fibonacci_heap.hpp> 
#include <dgraph_v_of_v/dgraph_v_of_v.h>
/*in a DBLP graph with 150k vertices and 500k edges,
the below code based on pairing_heap uses 0.72s;
the below code based on fibonacci_heap uses 0.72s;
while the boost code uses 0.2s;
the slowness is likely due to the frequent use of unordered_map*/

/*
struct node_for_dij
{
  public:
    int vertex, parent_vertex;
    double priority_value;
};

bool operator<(node_for_dij const &x, node_for_dij const &y)
{
    return x.priority_value > y.priority_value; // < is the max-heap; > is the min heap
}
*/
typedef typename boost::heap::fibonacci_heap<node_for_dij>::handle_type dgraph_heap_pointer;

//use dgraph_dijstar.h but a little change
void dgraph_shortest_paths_source_to_all(dgraph_v_of_v<two_hop_weight_type> &input_graph, int v_k,
                                         unordered_map<int, two_hop_weight_type> &distances)
{   
    int N = input_graph.INs.size();

    for (int i=0;i<N;i++)
    {
        distances[i] = std::numeric_limits<two_hop_weight_type>::max();
    }
    distances[v_k] = 0;

    queue<int> dij_dist_changed_vertices;
    vector<dgraph_heap_pointer> Q_pointer(N);
    node_for_dij node;
    boost::heap::fibonacci_heap<node_for_dij> Q;

    node.vertex = v_k;

    node.priority_value = 0;

    Q_pointer[v_k] = Q.push(node);


    while (Q.size() > 0)
    {
        node = Q.top();
        Q.pop();

        int u = node.vertex;
       
        two_hop_weight_type dist_u = node.priority_value;

        int u_adj_size = input_graph.OUTs[u].size();
        for (int i=0;i<u_adj_size;i++)
        {
            int adj_v = input_graph.OUTs[u][i].first;
            two_hop_weight_type new_dist = input_graph.OUTs[u][i].second + dist_u;

            if (distances[adj_v] == std::numeric_limits<two_hop_weight_type>::max())
            {
                node.vertex = adj_v;
                node.priority_value = new_dist;
                Q_pointer[adj_v] = Q.push(node);
                distances[adj_v] = node.priority_value;           
            }
            else
            {
                if (distances[adj_v] > new_dist)
                {
                    node.vertex = adj_v;
                    node.priority_value = new_dist;
                    Q.update(Q_pointer[adj_v], node);
                    distances[adj_v] = node.priority_value;                
                }
            }
        }
    }
}