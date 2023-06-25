#pragma once
#include <graph_hash_of_mixed_weighted/HL/two_hop_v1/graph_hash_of_mixed_weighted_CT_v1.h>


/*
suppose that 0,...,max_non_dummy_ID are non-dummy vertices, max_non_dummy_ID+1,...,max_N_ID-1 are dummy vertices


CT_dummy_v1:
all dummy vertices are in the core;  
(see code in MDE:  

if (nd.vertex > max_non_dummy_ID) {
			nd.degree = INT_MAX;
		}

see code in step 3: generate CT-tree indexs:

if (!isIntree[nd.vertex] && !popped_isIncore_node[nd.vertex] && (ideal_graph[nd.vertex].size() == nd.degree || nd.degree == INT_MAX) ) break;

) 
no LCA info;
only use PLL_dummy_v1 to generate core-index;

*/





/*CT_dummy_v1 seems wrong, since LCA may be needed even when all dummy vertices are in core*/

void CT_dummy_v1(graph_hash_of_mixed_weighted& input_graph, int max_non_dummy_ID, int max_N_ID, graph_hash_of_mixed_weighted_CT_v1_case_info& case_info) {


	//--------------------------------- step 1: initialization ---------------------------
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
	vector<vector<midnode>> midnode_graph(max_N_ID); // midnode_graph[i][j].left 离i近，midnode_graph[i][j].right离j近
	vector<vector<int>> midnode_bag(max_N_ID);
	graph_v_of_v_idealID ideal_graph = graph_hash_of_mixed_weighted_to_ideal_graph_of_CT(input_graph, midnode_graph, max_N_ID);
	isIntree.resize(max_N_ID, 0); // whether it is in the CT-tree

	/*
	priority_queue for maintaining the degrees of vertices  (we do not update degrees in q, so everytime you pop out a degree in q, you check whether it is the right one, ignore it if wrong)

	We need to maintain the degree of vertexes in graph dynamically, otherwise it would increase the cost in time
	*/
	priority_queue<node_degree> q;
	for (int i = 0; i < N; i++) // change N to max_N_ID cause bugs, why?
	{
		node_degree nd;
		nd.degree = ideal_graph[i].size();
		nd.vertex = i;

		if (nd.vertex > max_non_dummy_ID) {
			nd.degree = INT_MAX;
		}

		q.push(nd);
	}

	auto end1 = std::chrono::high_resolution_clock::now();
	case_info.time1_initialization = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9; // s
	//---------------------------------------------------------------------------------------------------




	
	//-------------------------------------------------- step 2: MDE-based tree decomposition ------------------------------------------------------------
	auto begin2 = std::chrono::high_resolution_clock::now();

	/*MDE-based tree decomposition; generating bags*/
	int bound_lambda = N;
	Bags.resize(N);
	vector<int> node_order(N + 1); // merging ID to original ID
	for (int i = 1; i <= N; i++)  // repeat
	{
		node_degree nd;
		while (1)
		{
			nd = q.top();
			q.pop();
			if (!isIntree[nd.vertex] && ideal_graph[nd.vertex].size() == nd.degree) break; // nd.vertex is the lowest degree non_dummy vertex not in tree
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

		vector<pair<int, double>> adj_temp = ideal_graph[v_x];
		vector<int> right_temp;
		int v_adj_size = ideal_graph[v_x].size();
		for (int j = 0; j < v_adj_size; j++)
		{
			auto it = &ideal_graph[v_x][j];
			Bags[v_x].push_back({ it->first,it->second });
			auto it2 = &midnode_graph[v_x][j];
			right_temp.push_back(it2->right);
			midnode_bag[v_x].push_back(it2->left); // midnode_bag对应Bags？
		}


		/*add new edge; this part has large time complexity and can be further optimized*/
		for (int j = 0; j < v_adj_size; j++)
			for (int k = 0; k < v_adj_size; k++)
				if (j != k)
				{
					substitute(ideal_graph, adj_temp[j].first, adj_temp[k].first, adj_temp[j].second + adj_temp[k].second, midnode_graph, right_temp[j], right_temp[k]);
					substitute(ideal_graph, adj_temp[k].first, adj_temp[j].first, adj_temp[j].second + adj_temp[k].second, midnode_graph, right_temp[k], right_temp[j]);
				}

		// delete edge related to v_x;
		for (int j = 0; j < v_adj_size; j++)
			remove(ideal_graph, ideal_graph[v_x][j].first, v_x, ideal_graph[v_x][j].second, midnode_graph);


		// delete v_x from ideal graph directly
		vector<pair<int, double>>().swap(ideal_graph[v_x]);
		vector<midnode>().swap(midnode_graph[v_x]);

		for (int j = 0; j < v_adj_size; j++)
		{
			struct node_degree nd;
			nd.vertex = adj_temp[j].first;
			nd.degree = ideal_graph[nd.vertex].size();
			//   q.update(q_handle[v_temp[j]], nd);

			if (nd.vertex > max_non_dummy_ID) {
				nd.degree = INT_MAX;
			}

			q.push(nd);
		}
	}


	auto end2 = std::chrono::high_resolution_clock::now();
	case_info.time2_tree_decomposition = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2).count() / 1e9; // s
	//--------------------------------------------------------------------------------------------------------


	//---------------------------------------------- step 3: generate CT-tree indexs ------------------------------------------------
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
			if (!isIntree[nd.vertex] && !popped_isIncore_node[nd.vertex] && (ideal_graph[nd.vertex].size() == nd.degree || nd.degree == INT_MAX) ) break;
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


	auto end3 = std::chrono::high_resolution_clock::now();
	case_info.time3_tree_indexs = std::chrono::duration_cast<std::chrono::nanoseconds>(end3 - begin3).count() / 1e9; // s
	//-------------------------------------------------------------------------------------------------------



	


	//------------------------------------------------ step 4: LCA --------------------------------------------------------------
	case_info.time4_lca = 0; // s
	//--------------------------------------------------------------------------------------------------------------------------







	//----------------------------------------------- step 5: 2-hop labeling -------------------------------------------
	auto begin5 = std::chrono::high_resolution_clock::now();

	/* construct 2-hop labels on core */
	/* graph_v_of_v_idealID to graph_hash_of_mixed_weighted */
	graph_hash_of_mixed_weighted hash_g;
	int size = ideal_graph.size();
	for (int i = 0; i < size; i++) {
		int v_size = ideal_graph[i].size();
		for (int x = 0; x < v_size; x++) {
			int j = ideal_graph[i][x].first;
			if (i < j) {
				double ec = ideal_graph[i][x].second;
				graph_hash_of_mixed_weighted_add_edge(hash_g, i, j, ec);
			}
		}
	}
	case_info.core_graph = hash_g;
	PLL_dummy_v1(hash_g, max_non_dummy_ID, max_N_ID, 1, case_info.thread_num, case_info.two_hop_case_info);

	auto end5 = std::chrono::high_resolution_clock::now();
	case_info.time5_core_indexs = std::chrono::duration_cast<std::chrono::nanoseconds>(end5 - begin5).count() / 1e9; // s
	//--------------------------------------------------------------------------------------------------------------------




	//-------------------------------------------------- step 6: postprocessing -------------------------------------------------------------------
	auto begin6 = std::chrono::high_resolution_clock::now();

	/* merge tree_index: L1 into case_info.two_hop_case_info.L */
	for (int v_k = 0; v_k < N; v_k++)
	{
		if (L1[v_k].size() > 0) {
			case_info.two_hop_case_info.L[v_k] = L1[v_k];
		}
	}

	/*build merging_predecessors*/
	std::vector<std::vector<pair<int, midnode>>>& merging_predecessors = case_info.merging_predecessors;
	merging_predecessors.resize(max_N_ID);
	for (int i = 0; i < max_N_ID; i++) {
		int size = ideal_graph[i].size();
		auto it1 = ideal_graph[i].begin();
		auto it2 = midnode_graph[i].begin();
		for (int x = 0; x < size; x++, it1++, it2++) {
			int j = it1->first;
			midnode node = *it2;
			graph_hash_of_mixed_weighted_binary_operations_insert(merging_predecessors[i], j, node);
		}
	}

	auto end6 = std::chrono::high_resolution_clock::now();
	case_info.time6_post = std::chrono::duration_cast<std::chrono::nanoseconds>(end6 - begin6).count() / 1e9; // s
	//---------------------------------------------------------------------------------------------------------------------------------
}








