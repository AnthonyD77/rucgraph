#pragma once
// #define DEBUG_MODE
// #define FIX_GRAPH
// #define CHECK_PROCESS

// #ifndef weightTYPE
// 	#define weightTYPE double
// #endif

class graph_edge_change{
    public:
    int source,terminal;
	int terminal_id;
    weightTYPE old_weight, new_weight;
    graph_edge_change(){}
	graph_edge_change(int _source, int _terminal, int _terminal_id, weightTYPE _old_weight, weightTYPE _new_weight) {
		source = _source;
		terminal = _terminal;
		terminal_id = _terminal_id;
		old_weight = _old_weight;
		new_weight = _new_weight;
	}
};

vector<vector<pair<weightTYPE, int>>> Dis_batch;
vector<vector<weightTYPE>> Q_value_batch;

void initialize_global_values_dynamic_batch(int N, int thread_num) {
	Dis_batch.resize(thread_num);
	Q_value_batch.resize(thread_num);
	for (int i = 0; i < thread_num; i++) {
		Dis_batch[i].resize(N, {-1, -1});
		Q_value_batch[i].resize(N, MAX_VALUE);
	}
}