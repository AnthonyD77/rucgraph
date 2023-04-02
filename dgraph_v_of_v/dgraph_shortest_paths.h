#pragma once

#include <unordered_map>
#include <boost/heap/fibonacci_heap.hpp> 
#include <dgraph_v_of_v/dgraph_v_of_v.h>

struct node_for_shortest_distances
{
public:
    int vertex;
    double priority_value;
};
bool operator<(node_for_shortest_distances const& x, node_for_shortest_distances const& y)
{
    return x.priority_value > y.priority_value; // < is the max-heap; > is the min heap
}
typedef typename boost::heap::fibonacci_heap<node_for_shortest_distances>::handle_type heap_pointer_dgraph_shortest_distances_source_to_all;

template <typename weight_type>
vector<weight_type> dgraph_shortest_distances_source_to_all(dgraph_v_of_v<weight_type> &input_graph, int v_k) {

    int N = input_graph.INs.size();
    vector<weight_type> distances(N, std::numeric_limits<weight_type>::max());
      
    node_for_shortest_distances node;
    boost::heap::fibonacci_heap<node_for_shortest_distances> Q;
    vector<heap_pointer_dgraph_shortest_distances_source_to_all> Q_pointer(N);

    distances[v_k] = 0;
    node.vertex = v_k;
    node.priority_value = 0;
    Q_pointer[v_k] = Q.push(node);

    while (Q.size() > 0) {
        node = Q.top();
        Q.pop();
        int u = node.vertex;
       
        weight_type dist_u = node.priority_value;

        int u_adj_size = input_graph.OUTs[u].size();
        for (int i=0;i<u_adj_size;i++) {
            int adj_v = input_graph.OUTs[u][i].first;
            weight_type new_dist = input_graph.OUTs[u][i].second + dist_u;

            if (distances[adj_v] == std::numeric_limits<weight_type>::max()) {
                node.vertex = adj_v;
                node.priority_value = new_dist;
                Q_pointer[adj_v] = Q.push(node);
                distances[adj_v] = node.priority_value;           
            }
            else {
                if (distances[adj_v] > new_dist) {
                    node.vertex = adj_v;
                    node.priority_value = new_dist;
                    Q.update(Q_pointer[adj_v], node);
                    distances[adj_v] = node.priority_value;                
                }
            }
        }
    }

    return distances;
}