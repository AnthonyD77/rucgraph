#pragma once
#include <queue>
#include <vector>
#include <climits>
#include <cmath>
#include <fstream>
#include <graph_v_of_v_idealID/graph_v_of_v_idealID.h>
#include <graph_hash_of_mixed_weighted/HL/two_hop_v1/graph_hash_of_mixed_weighted_two_hop_labels_v1.h>
#include <graph_hash_of_mixed_weighted/HL/two_hop_v1/graph_hash_of_mixed_weighted_PLL_v1.h>
#include <graph_hash_of_mixed_weighted/HL/two_hop_v1/graph_hash_of_mixed_weighted_PSL_v1.h>
#include <graph_hash_of_mixed_weighted/HL/two_hop_v1/graph_hash_of_mixed_weighted_PLL_dummy_v1.h>







struct midnode
{
	int left;
	int right;
};

struct node_degree
{
	int vertex, degree;
	bool operator < (const node_degree& nd) const
	{
		return degree > nd.degree;
	}
};

class graph_hash_of_mixed_weighted_CT_v1_case_info {
public:

	/*parameters*/
	int thread_num = 1;
	bool use_PLL = true; // 0 use PSL
	int max_non_dummy_ID = INT_MAX;
	int d = 3;
	// set d as an large number to test the correctness of CT-index
	// set d as 0 to test the correctness of PLL on core
	// set d as a moderate number to test the correctness of the mixed labels

	/*labels*/
	graph_hash_of_mixed_weighted_two_hop_case_info_v1 two_hop_case_info;
	std::vector<std::vector<pair<int, midnode>>> merging_predecessors; // used in case 1 of query path
	graph_hash_of_mixed_weighted core_graph; // in CT, it is required to save the merged graph for querying using R1 and R2
	std::vector<std::vector<pair<int, float>>> Bags; // bag nodes of decomposied tree
	vector<bool> isIntree; // isIntree[v]=1 means vertex v is in the tree_index
	vector<int> root;
	vector<vector<int>> tree_st; // for lca
	vector<vector<int>> tree_st_r; // for lca
	vector<int> first_pos; // for lca
	vector<int> lg; // for lca
	vector<int> dep; // for lca

	/*running limits*/
	long long int max_bit_size = 1e12;
	double max_run_time_seconds = 1e12; // s


	/*compute label size*/
	long long int compute_label_bit_size() {
		long long int size = 0;
		size = size + two_hop_case_info.compute_label_bit_size();  // L includes both tree and core indexes
		for (auto it = merging_predecessors.begin(); it != merging_predecessors.end(); it++) {
			size = size + (*it).size() * sizeof(pair<int, midnode>);
		}
		size = size + graph_hash_of_mixed_weighted_num_vertices(core_graph) * sizeof(std::vector <pair<int, double>>) + graph_hash_of_mixed_weighted_num_edges(core_graph) * sizeof(pair<int, double>);
		for (auto it = Bags.begin(); it != Bags.end(); it++) {
			size = size + (*it).size() * sizeof(pair<int, float>);
		}
		size = size + isIntree.size() * sizeof(bool);
		size = size + root.size() * 4;
		for (auto it = tree_st.begin(); it != tree_st.end(); it++) {
			size = size + (*it).size() * 4;
		}
		for (auto it = tree_st_r.begin(); it != tree_st_r.end(); it++) {
			size = size + (*it).size() * 4;
		}
		size = size + first_pos.size() * 4;
		size = size + lg.size() * 4;
		size = size + dep.size() * 4;
		return size;
	}

	/*clear labels*/
	void clear_labels() {
		two_hop_case_info.clear_labels();
		std::vector<std::vector<pair<int, midnode>>>().swap(merging_predecessors);
		core_graph.clear();
		std::vector<std::vector<pair<int, float>>>().swap(Bags);
		vector<bool>().swap(isIntree);
		vector<int>().swap(root);
		vector<vector<int>>().swap(tree_st);
		vector<vector<int>>().swap(tree_st_r);
		vector<int>().swap(first_pos);
		vector<int>().swap(lg);
		vector<int>().swap(dep);
	}


	/*indexing times*/
	double time1_initialization = 0;
	double time2_tree_decomposition = 0;
	double time3_tree_indexs = 0;
	double time4_lca = 0;
	double time5_core_indexs = 0;
	double time6_post = 0;
	double time_total = 0;


	/*printing*/
	void print_root() {
		cout << "print_root:" << endl;
		for (int i = 0; i < root.size(); i++) {
			cout << "root[" << i << "]=" << root[i] << endl;
		}
	}
	void print_isIntree() {
		cout << "print_isIntree:" << endl;
		for (int i = 0; i < isIntree.size(); i++) {
			cout << "isIntree[" << i << "]=" << isIntree[i] << endl;
		}
	}


	/*record_all_details*/
	void record_all_details(string save_name) {

		ofstream outputFile;
		outputFile.precision(6);
		outputFile.setf(ios::fixed);
		outputFile.setf(ios::showpoint);
		outputFile.open(save_name + ".txt");


		outputFile << "CT info:" << endl;
		outputFile << "thread_num=" << thread_num << endl;
		outputFile << "use_PLL=" << use_PLL << endl;
		outputFile << "max_non_dummy_ID=" << max_non_dummy_ID << endl;
		outputFile << "d=" << d << endl;

		outputFile << "max_bit_size=" << max_bit_size << endl;
		outputFile << "max_run_time_seconds=" << max_run_time_seconds << endl;

		outputFile << "compute_label_bit_size()=" << compute_label_bit_size() << endl;

		outputFile << "time1_initialization=" << time1_initialization << endl;
		outputFile << "time2_tree_decomposition=" << time2_tree_decomposition << endl;
		outputFile << "time3_tree_indexs=" << time3_tree_indexs << endl;
		outputFile << "time4_lca=" << time4_lca << endl;
		outputFile << "time5_core_indexs=" << time5_core_indexs << endl;
		outputFile << "time6_post=" << time6_post << endl;
		outputFile << "time_total=" << time_total << endl;




		outputFile << endl << endl;
		outputFile << "Core info:" << endl;
		outputFile << "use_2019R1=" << two_hop_case_info.use_2019R1 << endl;
		outputFile << "use_2019R2=" << two_hop_case_info.use_2019R2 << endl;
		outputFile << "use_enhanced2019R2=" << two_hop_case_info.use_enhanced2019R2 << endl;
		outputFile << "use_non_adj_reduc_degree=" << two_hop_case_info.use_non_adj_reduc_degree << endl;
		outputFile << "max_degree_MG_enhanced2019R2=" << two_hop_case_info.max_degree_MG_enhanced2019R2 << endl;
		outputFile << "use_dummy_dij_search_in_PLL=" << two_hop_case_info.use_dummy_dij_search_in_PLL << endl;
		outputFile << "reduce_V_num_2019R1=" << two_hop_case_info.reduce_V_num_2019R1 << endl;
		outputFile << "MG_num=" << two_hop_case_info.MG_num << endl;

		outputFile << "time_2019R1=" << two_hop_case_info.time_2019R1 << endl;
		outputFile << "time_2019R2_or_enhanced_pre=" << two_hop_case_info.time_2019R2_or_enhanced_pre << endl;
		outputFile << "time_2019R2_or_enhanced_fixlabels=" << two_hop_case_info.time_2019R2_or_enhanced_fixlabels << endl;
		outputFile << "time_initialization=" << two_hop_case_info.time_initialization << endl;
		outputFile << "time_generate_labels=" << two_hop_case_info.time_generate_labels << endl;
		outputFile << "time_canonical_repair1=" << two_hop_case_info.time_canonical_repair1 << endl;
		outputFile << "time_canonical_repair2=" << two_hop_case_info.time_canonical_repair2 << endl;
		outputFile << "time_update_old_IDs_in_labels=" << two_hop_case_info.time_update_old_IDs_in_labels << endl;

		outputFile << "max_labal_size=" << two_hop_case_info.max_labal_size << endl;
		outputFile << "max_run_time_seconds=" << two_hop_case_info.max_run_time_seconds << endl;

		outputFile << "use_canonical_repair=" << two_hop_case_info.use_canonical_repair << endl;
		outputFile << "label_size_before_canonical_repair=" << two_hop_case_info.label_size_before_canonical_repair << endl;
		outputFile << "label_size_after_canonical_repair=" << two_hop_case_info.label_size_after_canonical_repair << endl;
		outputFile << "canonical_repair_remove_label_ratio=" << two_hop_case_info.canonical_repair_remove_label_ratio << endl;

		outputFile << "compute_label_bit_size()=" << two_hop_case_info.compute_label_bit_size() << endl;



	}


};


/*global values*/
graph_v_of_v_idealID global_ideal_graph_CT;
vector<vector<midnode>> global_midnode_graph;  // midnode_graph[i][j].left 离i近，midnode_graph[i][j].right离j近

void clear_gloval_values_CT() {
	graph_v_of_v_idealID().swap(global_ideal_graph_CT);
	vector<vector<midnode>>().swap(global_midnode_graph);
}



void graph_hash_of_mixed_weighted_to_ideal_graph_of_CT(graph_hash_of_mixed_weighted& input_graph, int max_N_ID)
{
	global_ideal_graph_CT.resize(max_N_ID);
	global_midnode_graph.resize(max_N_ID);

	for (auto it1 = input_graph.hash_of_vectors.begin(); it1 != input_graph.hash_of_vectors.end(); it1++) {
		int i = it1->first;
		auto search = input_graph.hash_of_hashs.find(i);
		if (search != input_graph.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int j = it2->first;
				if (i < j) {
					double ec = it2->second;
					graph_v_of_v_idealID_add_edge(global_ideal_graph_CT, i, j, ec);
				}
			}
		}
		else {
			auto search2 = input_graph.hash_of_vectors.find(i);
			for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
				int j = it2->first;
				if (i < j) {
					double ec = it2->second;
					graph_v_of_v_idealID_add_edge(global_ideal_graph_CT, i, j, ec);
				}
			}
		}
	}
	for (int i = 0; i < max_N_ID; i++) {
		int size = global_ideal_graph_CT[i].size();
		for (int x = 0; x < size; x++) {
			midnode md1 = { global_ideal_graph_CT[i][x].first, i };
			global_midnode_graph[i].push_back(md1);
		}
	}
}

void substitute(int v_x, int v_y, double ec, int l, int r)
{
	//int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(ideal_graph[v_x], v_y);
	//if (pos == -1) { // add an new edge
	//	graph_v_of_v_idealID_add_edge(ideal_graph, v_x, v_y, ec);
	//	pos = graph_hash_of_mixed_weighted_binary_operations_search_position(ideal_graph[v_x], v_y);
	//	midnode_graph[v_x].insert(midnode_graph[v_x].begin() + pos, { l,r });
	//	pos = graph_hash_of_mixed_weighted_binary_operations_search_position(ideal_graph[v_y], v_x);
	//	midnode_graph[v_y].insert(midnode_graph[v_y].begin() + pos, { r,l });
	//}
	//else {
	//	double old_ec = ideal_graph[v_x][pos].second;
	//	if (ec < old_ec) {
	//		graph_v_of_v_idealID_add_edge(ideal_graph, v_x, v_y, ec);
	//		midnode_graph[v_x].insert(midnode_graph[v_x].begin() + pos, { l,r });
	//		pos = graph_hash_of_mixed_weighted_binary_operations_search_position(ideal_graph[v_y], v_x);
	//		midnode_graph[v_y].insert(midnode_graph[v_y].begin() + pos, { r,l });
	//	}
	//}

	//int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(ideal_graph[v_x], v_y);
	//if (pos == -1) { // add an new edge
	//	graph_hash_of_mixed_weighted_binary_operations_insert(ideal_graph[v_x], v_y, ec);
	//	pos = graph_hash_of_mixed_weighted_binary_operations_search_position(ideal_graph[v_x], v_y);
	//	midnode md;
	//	md.left = l;
	//	md.right = r;
	//	midnode_graph[v_x].insert(midnode_graph[v_x].begin() + pos, md);
	//}
	//else {
	//	double old_ec = ideal_graph[v_x][pos].second;
	//	if (ec < old_ec) {
	//		ideal_graph[v_x][pos].second = ec;
	//		midnode md;
	//		md.left = l;
	//		md.right = r;
	//		midnode_graph[v_x].insert(midnode_graph[v_x].begin() + pos, md);
	//	}
	//}

	bool flag = 0; // check if the edge already exists
	int v_adj_size = global_ideal_graph_CT[v_x].size();
	for (int i = 0; i < v_adj_size; i++)
		if (global_ideal_graph_CT[v_x][i].first == v_y) // the edge already exists； since ideal_graph is not sorted here, searching edge is very slow here
		{
			flag = 1;
			if (ec < global_ideal_graph_CT[v_x][i].second) // need to substitute the length
			{
				global_ideal_graph_CT[v_x][i].second = ec;
				global_midnode_graph[v_x][i].left = l;
				global_midnode_graph[v_x][i].right = r;
			}
			break;
		}
	if (!flag) // otherwise add an new edge
	{
		global_ideal_graph_CT[v_x].push_back({ v_y, ec });
		global_midnode_graph[v_x].push_back({ l,r });
	}
}

void remove(int v_x, int v_y, double ec)
{
	pair<int, double> p1(v_y, ec);
	vector<pair<int, double>> ::iterator  it = find(global_ideal_graph_CT[v_x].begin(), global_ideal_graph_CT[v_x].end(), p1);
	if (it != global_ideal_graph_CT[v_x].end())
	{
		global_midnode_graph[v_x].erase(it - global_ideal_graph_CT[v_x].begin() + global_midnode_graph[v_x].begin());
		global_ideal_graph_CT[v_x].erase(it);
	}
}

void two_hop_label_v1_binary_insert(std::vector<two_hop_label_v1>& input_vector, int vertex, int parent_vertex, double distance) {

	/*insert <key, load> into vector, if key is already inside, then load is updated; time complexity O(log n + size()-position ), which is O(n) in the worst case, as
	the time complexity of inserting an element into a vector is the number of elements behind this element*/

	int left = 0, right = input_vector.size() - 1;

	while (left <= right) // it will be skept when input_vector.size() == 0
	{
		int mid = left + ((right - left) / 2); // mid is between left and right (may be equal); 
		if (input_vector[mid].vertex == vertex) { // vertex is already in input_vector, this cannot happen in PLL
			return;
		}
		else if (input_vector[mid].vertex > vertex) {
			right = mid - 1; // the elements after right are always either empty, or have larger keys than input key
		}
		else {
			left = mid + 1; // the elements before left are always either empty, or have smaller keys than input key
		}
	}

	/*the following code is used when key is not in vector, i.e., left > right, specifically, left = right + 1;
	the elements before left are always either empty, or have smaller keys than input key;
	the elements after right are always either empty, or have larger keys than input key;
	so, the input key should be insert between right and left at this moment*/
	two_hop_label_v1 x;
	x.vertex = vertex;
	x.parent_vertex = parent_vertex;
	x.distance = distance;
	input_vector.insert(input_vector.begin() + left, x);

}







/*indexing function*/
void dfs(int& total, vector<int>& first_pos, int x, vector<vector<int>>& son, vector<int>& dfn)
{
	total++;
	dfn[total] = x;
	first_pos[x] = total;
	int s_size = son[x].size();
	for (int i = 0; i < s_size; i++)
	{
		dfs(total, first_pos, son[x][i], son, dfn); // this is a recursive function
		total++; dfn[total] = x;
	}
}

void CT_v1(graph_hash_of_mixed_weighted& input_graph, int max_N_ID, graph_hash_of_mixed_weighted_CT_v1_case_info& case_info) {


	//--------------------------------- step 1: initialization ---------------------------
	//cout << "step 1: initialization" << endl;
	auto begin1 = std::chrono::high_resolution_clock::now();

	auto& Bags = case_info.Bags;
	auto& isIntree = case_info.isIntree;
	auto& root = case_info.root;
	auto& tree_st = case_info.tree_st;
	auto& tree_st_r = case_info.tree_st_r;
	auto& first_pos = case_info.first_pos;
	auto& lg = case_info.lg;
	auto& dep = case_info.dep;

	int N = input_graph.hash_of_vectors.size();

	// used to extract the path
	vector<vector<int>> midnode_bag(max_N_ID);
	graph_hash_of_mixed_weighted_to_ideal_graph_of_CT(input_graph, max_N_ID);
	isIntree.resize(max_N_ID, 0); // whether it is in the CT-tree

	/*
	priority_queue for maintaining the degrees of vertices  (we do not update degrees in q, so everytime you pop out a degree in q, you check whether it is the right one, ignore it if wrong)

	We need to maintain the degree of vertexes in graph dynamically, otherwise it would increase the cost in time
	*/
	priority_queue<node_degree> q;
	for (int i = 0; i < N; i++) // change N to max_N_ID cause bugs, why?
	{
		node_degree nd;
		nd.degree = global_ideal_graph_CT[i].size();
		nd.vertex = i;
		q.push(nd);
	}

	auto end1 = std::chrono::high_resolution_clock::now();
	case_info.time1_initialization = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9;
	//---------------------------------------------------------------------------------------------------






	//-------------------------------------------------- step 2: MDE-based tree decomposition ------------------------------------------------------------
	//cout << "step 2: MDE-based tree decomposition" << endl;
	auto begin2 = std::chrono::high_resolution_clock::now();


	/*MDE-based tree decomposition; generating bags*/
	int bound_lambda = N;
	Bags.resize(N);
	vector<int> right_temp(N);
	vector<int> node_order(N + 1); // merging ID to original ID
	//ThreadPool pool(case_info.thread_num);
	std::vector< std::future<int> > results; // return typename: xxx
	for (int i = 1; i <= N; i++)  // repeat
	{
		node_degree nd;
		while (1)
		{
			nd = q.top();
			q.pop();
			if (!isIntree[nd.vertex] && global_ideal_graph_CT[nd.vertex].size() == nd.degree) break; // nd.vertex is the lowest degree vertex not in tree
		}
		int v_x = nd.vertex;  // the node with the minimum degree in G 
		if (nd.degree >= case_info.d) // reach the boudary
		{
			bound_lambda = i - 1;
			q.push(nd);
			break;  // until |Ni| >= d
		}

		isIntree[v_x] = 1; // add to CT-tree
		node_order[i] = v_x;

		vector<pair<int, double>> adj_temp = global_ideal_graph_CT[v_x];
		int v_adj_size = global_ideal_graph_CT[v_x].size();
		for (int j = 0; j < v_adj_size; j++)
		{
			auto it = &global_ideal_graph_CT[v_x][j];
			Bags[v_x].push_back({ it->first,it->second });
			auto it2 = &global_midnode_graph[v_x][j];
			right_temp[j] = it2->right;
			midnode_bag[v_x].push_back(it2->left); // midnode_bag对应Bags？
		}

		/*add new edge; this part has large time complexity and can be further optimized*/
		for (int j = 0; j < v_adj_size; j++) {
			for (int k = j + 1; k < v_adj_size; k++) {
				int adj_j = adj_temp[j].first;
				int adj_k = adj_temp[k].first;
				double new_ec = adj_temp[j].second + adj_temp[k].second;
				int right_j = right_temp[j];
				int right_k = right_temp[k];
				substitute(adj_j, adj_k, new_ec, right_j, right_k);
				substitute(adj_k, adj_j, new_ec, right_k, right_j);
			}
		}

		/*下面的并行在d=20时更慢*/
		//for (int j = 0; j < v_adj_size; j++) {
		//	for (int k = j + 1; k < v_adj_size; k++) {
		//		int adj_j = adj_temp[j].first;
		//		int adj_k = adj_temp[k].first;
		//		double new_ec = adj_temp[j].second + adj_temp[k].second;
		//		int right_j = right_temp[j];
		//		int right_k = right_temp[k];
		//		results.emplace_back(
		//			pool.enqueue([adj_j, adj_k, new_ec, right_j, right_k] { // pass const type value j to thread; [] can be empty
		//				substitute(adj_j, adj_k, new_ec, right_j, right_k);
		//				substitute(adj_k, adj_j, new_ec, right_k, right_j);
		//				return 1;
		//				})
		//		);
		//	}
		//}
		//for (auto&& result : results) {
		//	result.get();
		//}
		//results.clear();

		// delete edge related to v_x;
		for (int j = 0; j < v_adj_size; j++)
			remove(global_ideal_graph_CT[v_x][j].first, v_x, global_ideal_graph_CT[v_x][j].second);
		// delete v_x from ideal graph directly
		vector<pair<int, double>>().swap(global_ideal_graph_CT[v_x]);
		vector<midnode>().swap(global_midnode_graph[v_x]);


		///*add new edge; this part has large time complexity and can be further optimized*/
		//for (int j = 0; j < v_adj_size; j++)
		//	for (int k = j + 1; k < v_adj_size; k++)
		//		substitute(global_ideal_graph_CT, adj_temp[j].first, adj_temp[k].first, adj_temp[j].second + adj_temp[k].second, midnode_graph, right_temp[j], right_temp[k]);
		//// delete v_x from ideal graph
		//for (auto it = global_ideal_graph_CT[v_x].begin(); it != global_ideal_graph_CT[v_x].end(); it++) {
		//	int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(global_ideal_graph_CT[it->first], v_x);
		//	global_ideal_graph_CT[it->first].erase(global_ideal_graph_CT[it->first].begin() + pos);
		//	midnode_graph[it->first].erase(midnode_graph[it->first].begin() + pos);
		//}
		//vector<pair<int, double>>().swap(global_ideal_graph_CT[v_x]);
		//vector<midnode>().swap(midnode_graph[v_x]);


		for (int j = 0; j < v_adj_size; j++)
		{
			struct node_degree nd;
			nd.vertex = adj_temp[j].first;
			nd.degree = global_ideal_graph_CT[nd.vertex].size();
			q.push(nd);
		}
	}


	auto end2 = std::chrono::high_resolution_clock::now();
	case_info.time2_tree_decomposition = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2).count() / 1e9; // s
	//--------------------------------------------------------------------------------------------------------


	//---------------------------------------------- step 3: generate CT-tree indexs ------------------------------------------------
	//cout << "step 3: generate CT-tree indexs" << endl;
	auto begin3 = std::chrono::high_resolution_clock::now();

	/* generate CT-tree indexs */
	vector<vector<two_hop_label_v1>> L1(N); // Labels for CT-tree index merge later, otherwise may increase query time on PLL
	vector<int> fa(N);
	root.resize(N);
	vector<double> temp_dis(N);
	vector<int> temp_parent(N);
	vector<int> order_mapping(N + 1); // original ID to merging ID
	for (int i = 1; i <= bound_lambda; i++)  order_mapping[node_order[i]] = i;
	vector<bool> popped_isIncore_node(N, 0);
	for (int i = bound_lambda + 1; i <= N; i++)  // take advantage of the priority queue
	{
		struct node_degree nd;
		while (1)
		{
			nd = q.top();
			q.pop();
			if (!isIntree[nd.vertex] && !popped_isIncore_node[nd.vertex] && global_ideal_graph_CT[nd.vertex].size() == nd.degree) break;
		}
		node_order[i] = nd.vertex; // merging ID to original ID
		popped_isIncore_node[nd.vertex] = 1;
		order_mapping[nd.vertex] = i; // original ID to merging ID
	}

	dep.resize(N);
	vector<int> islabel(N, 0);
	first_pos.resize(N);
	vector<vector<int>> son(N);
	vector<int> index_node(N);

	vector <int> isneighour(N);
	int neighbournum = 0;
	vector<double> T_temp_n(N);
	vector<int> T_temp_p(N);


	int labelnum = 0;

	for (int i = bound_lambda; i >= 1; i--)
	{
		int v_x = node_order[i]; // reflection

		fa[v_x] = INT_MAX;
		int v_adj_size = Bags[v_x].size();
		for (int j = 0; j < v_adj_size; j++)
			if (order_mapping[Bags[v_x][j].first] < fa[v_x]) fa[v_x] = order_mapping[Bags[v_x][j].first]; // renew fa[v_x]

		if (fa[v_x] > bound_lambda || v_adj_size == 0)  //a root in the forest
		{
			root[v_x] = v_x;  fa[v_x] = -1; dep[v_x] = 0;// multi_fa[v_x][0] = -1;

			two_hop_label_v1 xx;
			for (int j = 0; j < v_adj_size; j++) // (bound_lambda - 1) - local distance to the interface
			{
				xx.vertex = Bags[v_x][j].first;
				xx.distance = Bags[v_x][j].second;
				xx.parent_vertex = midnode_bag[v_x][j];
				L1[v_x].push_back(xx);
			}
		}
		else //a non-root in the forest
		{
			fa[v_x] = node_order[fa[v_x]];
			root[v_x] = root[fa[v_x]]; dep[v_x] = dep[fa[v_x]] + 1;// multi_fa[v_x][0] = fa[v_x];
			son[fa[v_x]].push_back(v_x);

			int index_node_num = 0; labelnum++;  // to ensure the order, we can not use "push_back"

			int root_adj_size = Bags[root[v_x]].size();

			for (int j = 0; j < root_adj_size; j++)   // put the labels to interface in the beginnig position
			{  // add interface
				islabel[Bags[root[v_x]][j].first] = labelnum;
				index_node_num++;
				index_node[index_node_num] = Bags[root[v_x]][j].first;
				temp_dis[Bags[root[v_x]][j].first] = std::numeric_limits<double>::max();
			}


			int v_y = v_x;

			while (fa[v_y] != -1)
			{
				// add ancestor
				if (islabel[fa[v_y]] != labelnum)
				{
					index_node_num++;
					index_node[index_node_num] = fa[v_y];
					islabel[fa[v_y]] = labelnum;
					temp_dis[fa[v_y]] = std::numeric_limits<double>::max();
				}
				v_y = fa[v_y];
			}


			for (int j = 0; j < v_adj_size; j++)
			{
				// add neighbours

				if (islabel[Bags[v_x][j].first] != labelnum)
				{
					islabel[Bags[v_x][j].first] = labelnum;
					index_node.push_back(Bags[v_x][j].first);
					temp_dis[Bags[v_x][j].first] = Bags[v_x][j].second;
					temp_parent[Bags[v_x][j].first] = midnode_bag[v_x][j];
				}
				else
				{
					temp_dis[Bags[v_x][j].first] = Bags[v_x][j].second;
					temp_parent[Bags[v_x][j].first] = midnode_bag[v_x][j];
				}
			}
			// query (bound_lambda - 1)_local_distance to ancestor or the interface through neighbours

			for (int j = 0; j < v_adj_size; j++)
				if (isIntree[Bags[v_x][j].first]) // isneighour and isintree --> can be used as an intermediate node to update labels
				{
					double dis_vj = Bags[v_x][j].second;
					int vj = Bags[v_x][j].first;
					int Lj_size = L1[vj].size();
					for (int k = 0; k < Lj_size; k++)  // update the (bound_lambda-1)_local_distance
					{
						if (islabel[L1[vj][k].vertex] == labelnum && dis_vj + L1[vj][k].distance < temp_dis[L1[vj][k].vertex])
						{
							temp_dis[L1[vj][k].vertex] = dis_vj + L1[vj][k].distance;
							temp_parent[L1[vj][k].vertex] = midnode_bag[v_x][j];
						}
					}
				}


			// add labels to L1
			L1[v_x].resize(index_node_num);
			for (int j = 1; j <= index_node_num; j++)
			{
				two_hop_label_v1 xx;
				xx.vertex = index_node[j];
				xx.distance = temp_dis[index_node[j]];
				xx.parent_vertex = temp_parent[index_node[j]];
				L1[v_x][j - 1] = xx;

			}

			// update conversely
			neighbournum++;
			for (int j = 0; j < v_adj_size; j++)
			{
				isneighour[Bags[v_x][j].first] = neighbournum;
				T_temp_n[Bags[v_x][j].first] = Bags[v_x][j].second;
				T_temp_p[Bags[v_x][j].first] = midnode_bag[v_x][j];
			}
			for (int j = 1; j <= index_node_num; j++)
			{
				int vj = index_node[j];
				int Lj_size = L1[vj].size();
				for (int k = 0; k < Lj_size; k++)
				{
					int vk = L1[vj][k].vertex;
					if ((isneighour[vk] == neighbournum) && (T_temp_n[vk] + L1[vj][k].distance < temp_dis[vj]))
					{
						temp_dis[vj] = T_temp_n[vk] + L1[vj][k].distance;
						temp_parent[vj] = T_temp_p[vk];
					}
				}
			}

			for (int j = 1; j <= index_node_num; j++)
				if (temp_dis[index_node[j]] < L1[v_x][j - 1].distance)
				{
					L1[v_x][j - 1].distance = temp_dis[index_node[j]];
					L1[v_x][j - 1].parent_vertex = temp_parent[index_node[j]];
				}

		}
	}



	/*add distance-0 labels to tree nodes; this is needed in querying functions*/
	two_hop_label_v1 node;
	for (int i = 0; i < N; i++) {
		if (isIntree[i]) {
			node.vertex = i;
			node.distance = 0;
			node.parent_vertex = i;
			L1[i].push_back(node);
		}
	}



	auto end3 = std::chrono::high_resolution_clock::now();
	case_info.time3_tree_indexs = std::chrono::duration_cast<std::chrono::nanoseconds>(end3 - begin3).count() / 1e9; // s
	//-------------------------------------------------------------------------------------------------------






	//------------------------------------------------ step 4: LCA --------------------------------------------------------------
	//cout << "step 4: LCA" << endl;
	auto begin4 = std::chrono::high_resolution_clock::now();

	/* LCA code; already get the root, the father and the depth, here is the preprocessing of querying LCA */
	int total = 0;
	vector<int> dfn(2 * N + 5);
	for (int i = 1; i <= bound_lambda; i++)
	{
		int v_x = node_order[i];
		if (root[v_x] == v_x) dfs(total, first_pos, v_x, son, dfn);
	}

	if (total > 0)
	{
		int multi_step = ceil(log(total) / log(2)) + 2;

		tree_st.resize(total + 5);
		tree_st_r.resize(total + 5);
		for (int i = 1; i <= total; i++)
		{
			tree_st[i].resize(multi_step + 2);
			tree_st_r[i].resize(multi_step + 2);
		}

		vector<int> pow_2(multi_step);

		pow_2[0] = 1;
		for (int i = 1; i < multi_step; i++)
			pow_2[i] = pow_2[i - 1] << 1;


		for (int i = 1; i <= total; i++)
		{
			tree_st[i][0] = dfn[i];
			tree_st_r[i][0] = dfn[i];
		}

		for (int j = 1; j < multi_step; j++)
			for (int i = 1; i <= total; i++)
			{
				int k = i + pow_2[j - 1];
				if (k > total || dep[tree_st[i][j - 1]] <= dep[tree_st[k][j - 1]])
					tree_st[i][j] = tree_st[i][j - 1];
				else tree_st[i][j] = tree_st[k][j - 1];
				k = i - pow_2[j - 1];
				if (k <= 0 || dep[tree_st_r[i][j - 1]] <= dep[tree_st_r[k][j - 1]])
					tree_st_r[i][j] = tree_st_r[i][j - 1];
				else  tree_st_r[i][j] = tree_st_r[k][j - 1];
			}

	}

	lg.resize(total + 1);
	for (int i = 1; i <= total; i++)
		lg[i] = floor(log(i) / log(2));



	/*clear variables not used below*/
	vector<vector<int>>().swap(midnode_bag);
	priority_queue<node_degree>().swap(q);
	vector<int>().swap(node_order);
	vector<int>().swap(fa);
	vector<double>().swap(temp_dis);
	vector<int>().swap(temp_parent);
	vector<int>().swap(order_mapping);
	vector<int>().swap(islabel);
	vector<vector<int>>().swap(son);
	vector<int>().swap(index_node);
	vector <int>().swap(isneighour);
	vector<double>().swap(T_temp_n);
	vector<int>().swap(T_temp_p);
	vector<int>().swap(dfn);
	vector<int>().swap(right_temp);

	auto end4 = std::chrono::high_resolution_clock::now();
	case_info.time4_lca = std::chrono::duration_cast<std::chrono::nanoseconds>(end4 - begin4).count() / 1e9; // s
	//--------------------------------------------------------------------------------------------------------------------------




	//----------------------------------------------- step 5: 2-hop labeling -------------------------------------------
	//cout << "step 5: 2-hop labeling" << endl;
	auto begin5 = std::chrono::high_resolution_clock::now();

	/*update limits*/
	double to_date_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end4 - begin1).count() / 1e9;
	case_info.two_hop_case_info.max_run_time_seconds = case_info.max_run_time_seconds - to_date_time;
	if (case_info.two_hop_case_info.max_run_time_seconds < 0) {
		throw reach_limit_error_string;
	}
	long long int to_date_bit_size = case_info.compute_label_bit_size();
	case_info.two_hop_case_info.max_labal_size = (case_info.max_bit_size - to_date_bit_size) / 12; // this is slightly inaccurate, since reduction measures of R1 R2 are not counted
	if (case_info.two_hop_case_info.max_labal_size < 0) {
		throw reach_limit_error_string;
	}


	/* construct 2-hop labels on core */
	/* graph_v_of_v_idealID to graph_hash_of_mixed_weighted */
	graph_hash_of_mixed_weighted hash_g;
	int size = global_ideal_graph_CT.size();
	for (int i = 0; i < size; i++) {
		int v_size = global_ideal_graph_CT[i].size();
		for (int x = 0; x < v_size; x++) {
			int j = global_ideal_graph_CT[i][x].first;
			if (i < j) {
				double ec = global_ideal_graph_CT[i][x].second;
				graph_hash_of_mixed_weighted_add_edge(hash_g, i, j, ec);
			}
		}
	}
	case_info.core_graph = hash_g;
	if (case_info.use_PLL) {
		graph_hash_of_mixed_weighted_PLL_v1(hash_g, max_N_ID, 1, case_info.thread_num, case_info.two_hop_case_info);
	}
	else {
		graph_hash_of_mixed_weighted_PSL_v1(hash_g, max_N_ID, case_info.thread_num, case_info.two_hop_case_info);
	}
	//case_info.two_hop_case_info.print_L();
	auto end5 = std::chrono::high_resolution_clock::now();
	case_info.time5_core_indexs = std::chrono::duration_cast<std::chrono::nanoseconds>(end5 - begin5).count() / 1e9;
	//--------------------------------------------------------------------------------------------------------------------




	//-------------------------------------------------- step 6: postprocessing -------------------------------------------------------------------
	//cout << "step 6: postprocessing" << endl;
	auto begin6 = std::chrono::high_resolution_clock::now();

	/* merge tree_index: L1 into case_info.two_hop_case_info.L */
	for (int v_k = 0; v_k < N; v_k++)
	{
		if (L1[v_k].size() > 0) {
			vector<two_hop_label_v1>(L1[v_k]).swap(L1[v_k]);
			case_info.two_hop_case_info.L[v_k] = L1[v_k];
			vector<two_hop_label_v1>().swap(L1[v_k]);
		}
	}

	/*build merging_predecessors*/
	std::vector<std::vector<pair<int, midnode>>>& merging_predecessors = case_info.merging_predecessors;
	merging_predecessors.resize(max_N_ID);
	for (int i = 0; i < max_N_ID; i++) {
		int size = global_ideal_graph_CT[i].size();
		auto it1 = global_ideal_graph_CT[i].begin();
		auto it2 = global_midnode_graph[i].begin();
		for (int x = 0; x < size; x++, it1++, it2++) {
			int j = it1->first;
			midnode node = *it2;
			graph_hash_of_mixed_weighted_binary_operations_insert(merging_predecessors[i], j, node);
		}
	}

	auto end6 = std::chrono::high_resolution_clock::now();
	case_info.time6_post = std::chrono::duration_cast<std::chrono::nanoseconds>(end6 - begin6).count() / 1e9;
	//---------------------------------------------------------------------------------------------------------------------------------

	clear_gloval_values_CT();

	case_info.time_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end6 - begin1).count() / 1e9;
}








/*
query function

note that, the predecessors in the tree part of L are true predecessors in the original graph,
while the predecessors in the core part of L are predecessors in the merged core-graph (true predecessors in the original graph should be retrived using search_midnode)
*/

int lca(graph_hash_of_mixed_weighted_CT_v1_case_info& case_info, int x, int y)
{
	auto& first_pos = case_info.first_pos;

	if (first_pos[x] > first_pos[y])
	{
		int t = x; x = y; y = t;
	}

	int len = first_pos[y] - first_pos[x] + 1;
	int j = case_info.lg[len];
	x = case_info.tree_st[first_pos[x]][j]; y = case_info.tree_st_r[first_pos[y]][j];
	if (case_info.dep[x] < case_info.dep[y]) return x; else return y; // return the vertex with minimum depth between x and y in the dfs sequence.
}

double CT_extract_distance(graph_hash_of_mixed_weighted_CT_v1_case_info& case_info, int source, int terminal)
{
	auto& L = case_info.two_hop_case_info.L;
	auto& Bags = case_info.Bags;
	auto& isIntree = case_info.isIntree;
	auto& root = case_info.root;

	if (source == terminal) return 0.0;

	double distance = std::numeric_limits<double>::max(); // if disconnected, return this large value

	if (!isIntree[source] && !isIntree[terminal]) // both on core, use PLL indexs only
	{
		//cout << 1 << endl;

		return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance
		(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, source, terminal);

	}
	else if ((!isIntree[source] && isIntree[terminal]) || (isIntree[source] && !isIntree[terminal])) // one on core and one on tree
	{
		//cout << 2 << endl; 
		// 3-hop, use the interface
		if (isIntree[source])
		{
			int t = source; source = terminal; terminal = t;
		}
		/*the following: source is in core, terminal is in tree*/
		int r = root[terminal];
		//cout << "source=" << source << endl;
		//cout << "terminal=" << terminal << endl;
		//cout << "r=" << r << endl;
		int r_size = Bags[r].size();
		for (int i = 0; i < r_size; i++)
		{
			if (L[terminal][i].distance > distance)  continue; // already exceed the present minimum distance

			int x = L[terminal][i].vertex;
			float d_dis = L[terminal][i].distance;

			float dis = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance
			(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, source, x);
			//cout << "x=" << x << endl;
			//cout << "d_dis=" << d_dis << endl;
			//cout << "dis=" << dis << endl;
			if (distance > dis + d_dis) distance = dis + d_dis;
		}

		return distance;
	}
	else if (root[source] != root[terminal])  // locate in different trees, have not used lamma 9 here yet
	{
		//cout << 3 << endl;

		int r_s = root[source];
		int r_s_size = Bags[r_s].size();
		int r_t = root[terminal];
		int r_t_size = Bags[r_t].size();

		for (int i = 0; i < r_s_size; i++) {
			int u = L[source][i].vertex;
			double d_s_u = L[source][i].distance;
			for (int j = 0; j < r_t_size; j++) {
				int w = L[terminal][j].vertex;
				double d_t_w = L[terminal][j].distance;
				double d_u_w = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance
				(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, u, w);
				double dis = d_s_u + d_t_w + d_u_w;
				if (dis < distance) {
					distance = dis;
				}
			}
		}

		return distance;
	}
	else
	{
		int grand = lca(case_info, source, terminal);

		/*compute d2*/
		unordered_set<int> Bk = { grand };
		for (int i = Bags[grand].size() - 1; i >= 0; i--) {
			Bk.insert(Bags[grand][i].first);
		}
		unordered_map<int, double> Bk_source_dis;
		for (int i = L[source].size() - 1; i >= 0; i--) {
			int v = L[source][i].vertex;
			if (Bk.count(v) > 0) {
				Bk_source_dis[v] = L[source][i].distance;
			}
		}
		for (int i = L[terminal].size() - 1; i >= 0; i--) {
			int v = L[terminal][i].vertex;
			if (Bk_source_dis.count(v) > 0) {
				double d2 = Bk_source_dis[v] + L[terminal][i].distance;
				if (distance > d2) {
					distance = d2;
				}
			}
		}

		/*compute d4*/
		int r = root[source]; // source and terminal have the same root
		int r_size = Bags[r].size();
		for (int i = 0; i < r_size; i++)
		{
			int u = L[source][i].vertex;
			double d_s_u = L[source][i].distance;
			for (int j = 0; j < r_size; j++)
			{
				int w = L[terminal][j].vertex;
				double d_t_w = L[terminal][j].distance;
				double d_u_w = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance
				(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, u, w);
				double d4 = d_s_u + d_u_w + d_t_w;
				if (distance > d4) {
					distance = d4;
				}
			}
		}

		return distance;
	}
}

/*
the predecessor of vi in the original graph for edge (vi, vj) in the merged graph is the left of the return value (only when (vi, vj) is a merged edge, i.e., the left of the return value is not vi);
the predecessor of vj in the original graph for edge (vi, vj) in the merged graph is the right of the return value (only when (vi, vj) is a merged edge, i.e., the right of the return value is not vj);
*/
midnode search_midnode(std::vector<std::vector<pair<int, midnode>>>& merging_predecessors, int vi, int vj) {

	/*we assume that edge (vi, vj) exists; time complexity O(log n)*/

	std::vector<pair<int, midnode>>& input_vector = merging_predecessors[vi];

	int left = 0, right = input_vector.size() - 1;

	while (left <= right) {
		int mid = left + ((right - left) / 2); // mid is between left and right (may be equal); 
		if (input_vector[mid].first == vj) {
			return input_vector[mid].second;
		}
		else if (input_vector[mid].first > vj) {
			right = mid - 1;
		}
		else {
			left = mid + 1;
		}
	}
}

void CT_extract_path(graph_hash_of_mixed_weighted_CT_v1_case_info& case_info, int source, int terminal, vector<pair<int, int>>& path) {

	/*may return INT_MAX, INT_MAX*/

	//cout << "CT_extract_path " << source << " " << terminal << endl;
	//getchar();

	auto& root = case_info.root;
	auto& isIntree = case_info.isIntree;
	auto& Bags = case_info.Bags;
	auto& L = case_info.two_hop_case_info.L;
	auto& merging_predecessors = case_info.merging_predecessors;

	/*if source and terminal are disconnected, then return empty vector; note that, if source==terminal, then also return empty vector*/
	if (source == terminal) return;

	int vector1_capped_v_parent = -1, vector2_capped_v_parent = -1;

	double distance = std::numeric_limits<double>::max();  // record the shortest distance

	if (!isIntree[source] && !isIntree[terminal]) // both on core, use PLL indexs only
	{
		//cout << "a" << endl;

		pair<int, int> added_edge = { INT_MAX, INT_MAX };
		pair<int, int> two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors
		(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, source, terminal);
		if (two_predecessors.first == source && two_predecessors.second == terminal) { // disconnected
			path.push_back({ INT_MAX , INT_MAX });
			return;
		}

		//cout << "two_predecessors " << two_predecessors.first << " " << two_predecessors.second << endl;

		if (source != two_predecessors.first)
		{
			midnode node = search_midnode(case_info.merging_predecessors, source, two_predecessors.first);
			if (two_predecessors.first != node.left) {
				two_predecessors.first = node.left; // node.left should be the neighbour of source in the SP from source to two_predecessors.first
			}
			path.push_back({ source , two_predecessors.first });
			added_edge = { source , two_predecessors.first };
			source = two_predecessors.first;
		}
		if (terminal != two_predecessors.second)
		{
			midnode node = search_midnode(case_info.merging_predecessors, terminal, two_predecessors.second);
			if (two_predecessors.second != node.left) {
				two_predecessors.second = node.left; // node.left should be the neighbour of terminal in the SP from terminal to two_predecessors.second
			}
			if (!(added_edge.first == two_predecessors.second && added_edge.second == terminal)) {
				path.push_back({ two_predecessors.second , terminal });
				terminal = two_predecessors.second; // else, two_predecessors.second == source, terminal should not change to source, since source has changed to terminal above
			}
			else {
				return;
			}
		}

		//cout<< "resursive CT_extract_path " << source << " " << terminal << endl;

		CT_extract_path(case_info, source, terminal, path);
		return;
	}
	else if ((!isIntree[source] && isIntree[terminal]) || (isIntree[source] && !isIntree[terminal])) // source on core and terminal on tree
	{
		//cout << "b" << endl;
		// 3-hop, use the interface
		if (isIntree[source])
		{
			int t = source; source = terminal; terminal = t;
		}
		/*the following: source is in core, terminal is in tree*/
		int r = root[terminal];
		int r_size = Bags[r].size();
		int terminal_predecessor = INT_MAX;
		//cout << "source=" << source << endl;
		//cout << "terminal=" << terminal << endl;
		//cout << "r=" << r << endl;
		for (int i = 0; i < r_size; i++)
		{
			if (L[terminal][i].distance > distance)  continue; // already exceed the present minimum distance

			int x = L[terminal][i].vertex;
			float d_dis = L[terminal][i].distance;

			float dis = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance
			(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, source, x);
			//cout << "x=" << x << endl;
			//cout << "d_dis=" << d_dis << endl;
			//cout << "dis=" << dis << endl;
			if (distance > dis + d_dis) {
				distance = dis + d_dis;
				terminal_predecessor = L[terminal][i].parent_vertex; // terminal is in the tree, so this is the predecessor in the original graph
			}
		}
		if (terminal_predecessor == INT_MAX) { // disconnected
			path.push_back({ INT_MAX , INT_MAX });
			return;
		}
		path.push_back({ terminal_predecessor, terminal });
		terminal = terminal_predecessor;
		CT_extract_path(case_info, source, terminal, path);

		return;
	}
	else  if (root[source] != root[terminal])  // locate in different trees
	{
		//cout << "c" << endl;
		int r_s = root[source];
		int r_s_size = Bags[r_s].size();
		int r_t = root[terminal];
		int r_t_size = Bags[r_t].size();
		int source_predecessor = INT_MAX, terminal_predecessor = INT_MAX;
		for (int i = 0; i < r_s_size; i++) {
			int u = L[source][i].vertex;
			double d_s_u = L[source][i].distance;
			for (int j = 0; j < r_t_size; j++) {
				int w = L[terminal][j].vertex;
				double d_t_w = L[terminal][j].distance;
				double d_u_w = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance
				(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, u, w);
				double dis = d_s_u + d_t_w + d_u_w;
				if (dis < distance) {
					distance = dis;
					source_predecessor = L[source][i].parent_vertex;
					terminal_predecessor = L[terminal][j].parent_vertex;
				}
			}
		}
		if (source_predecessor == INT_MAX || terminal_predecessor == INT_MAX) { // disconnected
			path.push_back({ INT_MAX , INT_MAX });
			return;
		}
		pair<int, int> added_edge = { INT_MAX, INT_MAX };
		if (source != source_predecessor)
		{
			path.push_back({ source , source_predecessor });
			added_edge = { source , source_predecessor };
			source = source_predecessor;
		}
		if (terminal != terminal_predecessor)
		{
			if (!(added_edge.first == terminal_predecessor && added_edge.second == terminal)) {
				path.push_back({ terminal_predecessor , terminal });
				terminal = terminal_predecessor;
			}
			else {
				return;
			}
		}
		CT_extract_path(case_info, source, terminal, path);

		return;
	}
	else
	{
		//cout << "d" << endl;
		int source_predecessor = INT_MAX, terminal_predecessor = INT_MAX;

		int grand = lca(case_info, source, terminal);

		/*compute d2*/
		unordered_set<int> Bk = { grand };
		for (int i = Bags[grand].size() - 1; i >= 0; i--) {
			Bk.insert(Bags[grand][i].first);
		}
		unordered_map<int, pair<double, int>> Bk_source_dis;
		for (int i = L[source].size() - 1; i >= 0; i--) {
			int v = L[source][i].vertex;
			if (Bk.count(v) > 0) {
				Bk_source_dis[v] = { L[source][i].distance , L[source][i].parent_vertex };
			}
		}
		for (int i = L[terminal].size() - 1; i >= 0; i--) {
			int v = L[terminal][i].vertex;
			if (Bk_source_dis.count(v) > 0) {
				double d2 = Bk_source_dis[v].first + L[terminal][i].distance;
				if (distance > d2) {
					distance = d2;
					source_predecessor = Bk_source_dis[v].second;
					terminal_predecessor = L[terminal][i].parent_vertex;
				}
			}
		}

		/*compute d4*/
		int r = root[source]; // source and terminal have the same root
		int r_size = Bags[r].size();
		for (int i = 0; i < r_size; i++)
		{
			int u = L[source][i].vertex;
			double d_s_u = L[source][i].distance;
			for (int j = 0; j < r_size; j++)
			{
				int w = L[terminal][j].vertex;
				double d_t_w = L[terminal][j].distance;
				double d_u_w = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance
				(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, u, w);
				double d4 = d_s_u + d_u_w + d_t_w;
				if (distance > d4) {
					distance = d4;
					source_predecessor = L[source][i].parent_vertex;
					terminal_predecessor = L[terminal][j].parent_vertex;
				}
			}
		}

		if (source_predecessor == INT_MAX || terminal_predecessor == INT_MAX) { // disconnected
			path.push_back({ INT_MAX , INT_MAX });
			return;
		}
		pair<int, int> added_edge = { INT_MAX, INT_MAX };
		if (source != source_predecessor)
		{
			path.push_back({ source , source_predecessor });
			added_edge = { source , source_predecessor };
			source = source_predecessor;
		}
		if (terminal != terminal_predecessor)
		{
			if (!(added_edge.first == terminal_predecessor && added_edge.second == terminal)) {
				path.push_back({ terminal_predecessor , terminal });
				terminal = terminal_predecessor;
			}
			else {
				return;
			}
		}
		CT_extract_path(case_info, source, terminal, path);


		return;

	}

}










