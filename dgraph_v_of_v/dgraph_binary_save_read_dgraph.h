#pragma once
#include <dgraph_v_of_v/dgraph_v_of_v.h>
#include <text_mining/binary_save_read_vector_of_vectors.h>

template <typename weight_type>
void dgraph_binary_save_dgraph(std::string save_path, dgraph_v_of_v<weight_type>& input_graph) {

    std::vector<std::vector<std::pair<int, weight_type>>> save_vectors = input_graph.INs;
    save_vectors.insert(save_vectors.end(), input_graph.OUTs.begin(), input_graph.OUTs.end()); // save INs+OUTs

    binary_save_vector_of_vectors(save_path, save_vectors);
}

template <typename weight_type>
void dgraph_binary_read_dgraph(std::string read_path, dgraph_v_of_v<weight_type>& read_graph) {

    std::vector<std::vector<std::pair<int, weight_type>>> save_vectors;
    binary_read_vector_of_vectors(read_path, save_vectors); // read INs+OUTs

    int N = save_vectors.size() / 2;
    read_graph.INs = { save_vectors.begin(), save_vectors.begin() + N };
    read_graph.OUTs = { save_vectors.begin() + N, save_vectors.begin() + 2 * N };
}




#include <dgraph_v_of_v/dgraph_compare_two_dgraphs.h>

void test_dgraph_binary_save_read_dgraph() {

    int iteration_graph_times = 100;
    int V = 100, E = 500, precision = 1;
    two_hop_weight_type ec_min = 1, ec_max = 10;

    for (int i = 0; i < iteration_graph_times; i++) {
        cout << i << endl;
        dgraph_v_of_v<two_hop_weight_type> instance_graph = dgraph_generate_random_dgraph(V, E, ec_min, ec_max, precision, boost_random_time_seed);
        dgraph_binary_save_dgraph("x.bin", instance_graph);
        dgraph_v_of_v<two_hop_weight_type> g;
        dgraph_binary_read_dgraph("x.bin", g);

        if (!dgraph_compare_two_dgraphs_not_eaxct_same_weight(g, instance_graph)) {
            cout << "instance_graph != g" << endl;
            getchar();
        }
    }
}