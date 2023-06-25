#pragma once
#include <vector>
#include <ThreadPool.h>
#include <shared_mutex>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>


/*PLL label format*/
class two_hop_label_v1 {
public:
	int vertex, parent_vertex; // std::vector<PLL_with_non_adj_reduction_sorted_label> should be sorted by vertex
	float distance;
};


/*
global values that should be cleared after usig PLL or PSL
unique code for this file: 595
*/
string reach_limit_error_string = "reach limit error";
long long int max_labal_size_595;
long long int labal_size_595;
auto begin_time_595 = std::chrono::high_resolution_clock::now();
double max_run_time_nanoseconds_595;
bool this_parallel_PLL_PSL_is_running_595 = false;
bool if_continue_595;
graph_v_of_v_idealID ideal_graph_595;
vector<vector<two_hop_label_v1>> L_595;
vector<vector<two_hop_label_v1>> L_temp_595; // d轮中d-1轮的label
std::shared_mutex mtx_595_1, mtx_595_2;
int max_N_ID_for_mtx_595 = 1e7;  // this is the max N to run
vector<std::shared_mutex> mtx_595(max_N_ID_for_mtx_595);  // std::mutex has no copy or move constructor, while std::vector::resize() requires that; you cannot resize mtx;    
// moreover, do not change mtx to a pointer and then points to local values, it is very slow!!
queue<int> Qid_595; // IDs of available elements of P T
vector<vector<float>> P_dij_595;
vector<vector<float>> T_dij_595;
vector<vector<int>> P_bfs_595;
vector<vector<int>> T_bfs_595;
vector<int> vertexID_new_to_old_595;
vector<int> pos_595;
vector<int> pos_2_595;
vector<int> increment_595;
vector<vector<bool>> dirty_tag_595;
vector<int> reduction_measures_2019R1_new_ID; 
vector<int> reduction_measures_2019R2_new_ID;
vector<int> f_2019R1_new_ID;
vector<vector<two_hop_label_v1>> incremental_label_vectors;
vector<vector<int>> labels_to_be_removed;
vector <vector<pair<int, double>>> adjs_new_IDs;
vector<pair<int, double>> min_adjs_new_IDs;

void graph_hash_of_mixed_weighted_two_hop_clear_global_values() {

	this_parallel_PLL_PSL_is_running_595 = false;
	graph_v_of_v_idealID().swap(ideal_graph_595);
	vector<vector<two_hop_label_v1>>().swap(L_595);
	vector<vector<two_hop_label_v1>>().swap(L_temp_595);
	queue<int>().swap(Qid_595);
	vector<vector<float>>().swap(P_dij_595);
	vector<vector<float>>().swap(T_dij_595);
	vector<vector<int>>().swap(P_bfs_595);
	vector<vector<int>>().swap(T_bfs_595);
	vector<int>().swap(vertexID_new_to_old_595);
	vector<int>().swap(pos_595);
	vector<int>().swap(pos_2_595);
	vector<int>().swap(increment_595);
	vector<vector<bool>>().swap(dirty_tag_595);
	vector<int>().swap(reduction_measures_2019R1_new_ID);
	vector<int>().swap(reduction_measures_2019R2_new_ID);
	vector<int>().swap(f_2019R1_new_ID);
	vector<vector<two_hop_label_v1>>().swap(incremental_label_vectors);
	vector<vector<int>>().swap(labels_to_be_removed);
	vector <vector<pair<int, double>>>().swap(adjs_new_IDs);
	vector<pair<int, double>>().swap(min_adjs_new_IDs);
}

/*
global querying values; 

every time PLL or PSL is used, the following global querying values are updated, and then the previous PLL or PSL indexes cannot be reused again
*/
vector <vector<pair<int, double>>> adjs;
vector<pair<int, double>> min_adjs;
bool full_two_hop_labels = true;

void graph_hash_of_mixed_weighted_two_hop_clear_global_values2() {

	vector <vector<pair<int, double>>>().swap(adjs);
	vector<pair<int, double>>().swap(min_adjs);
}







class graph_hash_of_mixed_weighted_two_hop_case_info_v1 {
public:

	/*use reduction info*/
	bool use_2019R1 = false;
	bool use_2019R2 = false;
	bool use_enhanced2019R2 = false;
	bool use_non_adj_reduc_degree = false; // use_enhanced2019R2 和 use_non_adj_reduc_degree 貌似效果基本一样
	bool use_dummy_dij_search_in_PLL = false;
	int max_degree_MG_enhanced2019R2 = 100;
	int reduce_V_num_2019R1 = 0, MG_num = 0;


	/*running time records*/
	double time_2019R1 = 0;
	double time_2019R2_or_enhanced_pre = 0;
	double time_2019R2_or_enhanced_fixlabels = 0; // s
	double time_initialization = 0;
	double time_generate_labels = 0;
	double time_canonical_repair1 = 0;
	double time_canonical_repair2 = 0;
	double time_update_old_IDs_in_labels = 0;

	/*running limits*/
	long long int max_labal_size = 1e12; // 2-hop-label-num
	double max_run_time_seconds = 1e12; // s

	/*labels*/
	vector<int> reduction_measures_2019R2; // for 2019 R2
	vector<int> reduction_measures_2019R1; // for 2019 R1;  11 means equivalent_1 relation (no edge between), 12 means equivalent_2 relation (edge between)
	vector<int> f_2019R1; // for 2019 R1
	vector<vector<two_hop_label_v1>> L;
	

	/*canonical_repair info*/
	bool use_canonical_repair = false;
	long long int label_size_before_canonical_repair = 0;
	long long int label_size_after_canonical_repair = 0;
	double canonical_repair_remove_label_ratio = 0; 


	/*compute label size; this should equal label_size_after_canonical_repair when use_canonical_repair==true*/
	long long int compute_label_bit_size() {
		long long int size = 0;
		size = size + reduction_measures_2019R2.size() * 4;
		size = size + reduction_measures_2019R1.size() * 4;
		size = size + f_2019R1.size() * 4;
		for (auto it = L.begin(); it != L.end(); it++) {
			size = size + (*it).size() * sizeof(two_hop_label_v1); // 12 bit per two_hop_label_v1
		}
		return size;
	}


	/*clear labels*/
	void clear_labels() {
		vector<int>().swap(reduction_measures_2019R2);
		vector<int>().swap(reduction_measures_2019R1);
		vector<int>().swap(f_2019R1);
		vector<vector<two_hop_label_v1>>().swap(L);	
	}



	/*printing*/
	void print_L() {
		cout << "print_L:" << endl;
		for (int i = 0; i < L.size(); i++) {
			cout << "L[" << i << "]=";
			for (int j = 0; j < L[i].size(); j++) {
				cout << "{" << L[i][j].vertex << "," << L[i][j].distance << "," << L[i][j].parent_vertex << "}";
			}
			cout << endl;
		}
	}
	void print_reduction_measures_2019R1() {
		cout << "print_reduction_measures_2019R1:" << endl;
		for (int i = 0; i < reduction_measures_2019R1.size(); i++) {
			cout << "reduction_measures_2019R1[" << i << "]=" << reduction_measures_2019R1[i] << endl;
		}
	}
	void print_reduction_measures_2019R2() {
		cout << "print_reduction_measures_2019R2:" << endl;
		for (int i = 0; i < reduction_measures_2019R2.size(); i++) {
			cout << "reduction_measures_2019R2[" << i << "]=" << reduction_measures_2019R2[i] << endl;
		}
	}
	void print_f_2019R1() {
		cout << "print_f_2019R1:" << endl;
		for (int i = 0; i < f_2019R1.size(); i++) {
			cout << "f_2019R1[" << i << "]=" << f_2019R1[i] << endl;
		}
	}

};







/*common functions shared by PLL and PSL*/
bool compare_pair_second_large_to_small(const pair<int, int>& i, pair<int, int>& j)
{
	/*< is nearly 10 times slower than >*/
	return i.second > j.second;  // < is from small to big; > is from big to small.  sort by the second item of pair<int, int>
}

bool compare_two_hop_label_small_to_large(two_hop_label_v1& i, two_hop_label_v1& j)
{
	return i.vertex < j.vertex;  // < is from small to big; > is from big to small
}

/*the following locks are used in PLL search process and canonical_repair*/

/*
canonical_repair

graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc_for_canonical_repair is modified to only use target_label_vector_before_checking_label;

graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1_for_canonical_repair and 
graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_for_canonical_repair are generally the same with the non "for_canonical_repair" versions
*/

float graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc_for_canonical_repair(int source, int terminal, int target_v) {

	/*return std::numeric_limits<float>::max() is not connected*/

	if (source == terminal) {
		return 0;
	}

	float distance = std::numeric_limits<float>::max(); // if disconnected, return this large value

	vector<two_hop_label_v1>::iterator vector1_check_pointer, vector2_check_pointer, pointer_L_s_end, pointer_L_t_end;
	if (source == target_v) {
		vector1_check_pointer = incremental_label_vectors[target_v].begin();
		pointer_L_s_end = incremental_label_vectors[target_v].end();
	}
	else {
		vector1_check_pointer = L_temp_595[source].begin();
		pointer_L_s_end = L_temp_595[source].end();
	}
	if (terminal == target_v) {
		vector2_check_pointer = incremental_label_vectors[target_v].begin();
		pointer_L_t_end = incremental_label_vectors[target_v].end();
	}
	else {
		vector2_check_pointer = L_temp_595[terminal].begin();
		pointer_L_t_end = L_temp_595[terminal].end();
	}

	while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end) {
		if (vector1_check_pointer->vertex == vector2_check_pointer->vertex) {
			float dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
			if (distance > dis) {
				distance = dis;
			}
			vector1_check_pointer++;
		}
		else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex) {
			vector2_check_pointer++;
		}
		else {
			vector1_check_pointer++;
		}
	}

	return distance;

}

float graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1_for_canonical_repair /* we assume that source and terminal are not reduced by 2019R1 here*/
(int source, int terminal, int target_v)
{

	/* we assume that source and terminal are not reduced by 2019R1*/

	if (source == terminal) {
		return 0;
	}

	vector<float> selected_distance = { std::numeric_limits<float>::max() }; // Store all path lengths to be selected; disconnected = std::numeric_limits<float>::max()
	auto s_adj_begin = adjs_new_IDs[source].begin();
	auto s_adj_end = adjs_new_IDs[source].end();
	auto t_adj_begin = adjs_new_IDs[terminal].begin();
	auto t_adj_end = adjs_new_IDs[terminal].end();
	if (reduction_measures_2019R2_new_ID[source] == 2)
	{
		if (reduction_measures_2019R2_new_ID[terminal] == 2)
		{
			/*"Both reduced"*/
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
				{
					if (f_2019R1_new_ID[it1->first] == it1->first && f_2019R1_new_ID[it2->first] == it2->first) {
						float x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc_for_canonical_repair(it1->first, it2->first, target_v);
						if (x == std::numeric_limits<float>::max()) {
							return x;
						}
						else {
							selected_distance.push_back(x + float(it1->second) + float(it2->second));
						}
					}
				}
			}
		}
		else
		{
			/*"Only source reduced"*/
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				if (f_2019R1_new_ID[it1->first] == it1->first) {
					float x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc_for_canonical_repair(it1->first, terminal, target_v);
					if (x == std::numeric_limits<float>::max()) {
						return x;
					}
					else {
						selected_distance.push_back(x + float(it1->second));
					}
				}
			}
		}
	}
	else
	{
		if (reduction_measures_2019R2_new_ID[terminal] == 2)
		{
			/*"Only terminal reduced"*/
			for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
			{
				if (f_2019R1_new_ID[it2->first] == it2->first) {
					float x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc_for_canonical_repair(source, it2->first, target_v);
					if (x == std::numeric_limits<float>::max()) {
						return x;
					}
					else {
						selected_distance.push_back(x + float(it2->second));
					}
				}
			}
		}
		else
		{
			/*"Nothing happened"*/
			selected_distance.push_back(graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc_for_canonical_repair(source, terminal, target_v));
		}
	}

	float dis = *min_element(selected_distance.begin(), selected_distance.end());

	return dis;
}

float graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_for_canonical_repair(graph_hash_of_mixed_weighted& instance_graph, int source, int terminal, int target_v)
{

	if (source == terminal) {
		return 0;
	}

	if (reduction_measures_2019R1_new_ID[source] == 11)
	{
		/* case 2 */
		if (reduction_measures_2019R1_new_ID[terminal] == 11)
		{
			if (f_2019R1_new_ID[source] == f_2019R1_new_ID[terminal])
			{
				pair<int, double> s_min_adj = min_adjs_new_IDs[source];
				return float(s_min_adj.second * 2);
			}
			else
			{
				return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1_for_canonical_repair
				(f_2019R1_new_ID[source], f_2019R1_new_ID[terminal], target_v);
			}
		}
		/* case 3 */
		else if (reduction_measures_2019R1_new_ID[terminal] == 12)
		{
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1_for_canonical_repair
			(f_2019R1_new_ID[source], f_2019R1_new_ID[terminal], target_v);
		}
		else { /* case 1 */
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1_for_canonical_repair
			(f_2019R1_new_ID[source], terminal, target_v);
		}
	}
	else if (reduction_measures_2019R1_new_ID[source] == 12)
	{
		/* case 5 */
		if (reduction_measures_2019R1_new_ID[terminal] == 11)
		{
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1_for_canonical_repair
			(f_2019R1_new_ID[source], f_2019R1_new_ID[terminal], target_v);
		}
		/* case 6 -- same with case 3 */
		else if (reduction_measures_2019R1_new_ID[terminal] == 12)
		{
			if (f_2019R1_new_ID[source] == f_2019R1_new_ID[terminal])
			{
				pair<int, double> s_min_adj = min_adjs_new_IDs[source];
				float s_t_weight = float(graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, terminal));
				return float(s_min_adj.second * 2) > s_t_weight ? s_t_weight : float(s_min_adj.second * 2);
			}
			else
			{
				return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1_for_canonical_repair
				(f_2019R1_new_ID[source], f_2019R1_new_ID[terminal], target_v);
			}
		}
		/* case 4 */
		else {
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1_for_canonical_repair
			(f_2019R1_new_ID[source], terminal, target_v);
		}
	}
	else {
		/* case 8 -- same with case 1 */
		if (reduction_measures_2019R1_new_ID[terminal] == 11)
		{
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1_for_canonical_repair
			(source, f_2019R1_new_ID[terminal], target_v);
		}
		/* case 9 -- same with case 4 */
		else if (reduction_measures_2019R1_new_ID[terminal] == 12)
		{
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1_for_canonical_repair
			(source, f_2019R1_new_ID[terminal], target_v);
		}
		else {
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1_for_canonical_repair
			(source, terminal, target_v);
		}
	}

}

void sorted_int_vector_binary_insert(std::vector<int>& input_vector, int key) {

	int left = 0, right = input_vector.size() - 1;

	while (left <= right) // it will be skept when input_vector.size() == 0
	{
		int mid = left + ((right - left) / 2); // mid is between left and right (may be equal); 
		if (input_vector[mid] == key) {
			return;
		}
		else if (input_vector[mid] > key) {
			right = mid - 1; // the elements after right are always either empty, or have larger keys than input key
		}
		else {
			left = mid + 1; // the elements before left are always either empty, or have smaller keys than input key
		}
	}
	input_vector.insert(input_vector.begin() + left, key);
}

bool sorted_int_vector_binary_search(std::vector<int>& input_vector, int key) {

	/*return true if key is in vector; time complexity O(log n)*/

	int left = 0, right = input_vector.size() - 1;

	while (left <= right) {
		int mid = left + ((right - left) / 2); // mid is between left and right (may be equal); 
		if (input_vector[mid] == key) {
			return true;
		}
		else if (input_vector[mid] > key) {
			right = mid - 1;
		}
		else {
			left = mid + 1;
		}
	}

	return false;
}

void canonical_repair_element1(graph_hash_of_mixed_weighted* instance_graph_pointer, int target_v) {

	auto& instance_graph = *instance_graph_pointer;

	auto begin = L_temp_595[target_v].begin(), end = L_temp_595[target_v].end();
	for (; begin != end; begin++) {
		int vertex = begin->vertex;
		if (vertex == target_v) {
			continue;
		}
		double distance = begin->distance;
		double query_dis = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_for_canonical_repair
		(instance_graph, target_v, vertex, target_v);

		//cout << "target_v " << target_v << " vertex " << vertex << " query_dis " << query_dis << endl;
		//cout << "reduction_measures_2019R1_new_ID[target_v] " << reduction_measures_2019R1_new_ID[target_v] << endl;
		//cout << "f_2019R1_new_ID[target_v] " << f_2019R1_new_ID[target_v] << endl;
		//cout << "reduction_measures_2019R2_new_ID[target_v] " << reduction_measures_2019R2_new_ID[target_v] << endl;
		//cout << "reduction_measures_2019R1_new_ID[vertex] " << reduction_measures_2019R1_new_ID[vertex] << endl;
		//cout << "reduction_measures_2019R2_new_ID[vertex] " << reduction_measures_2019R2_new_ID[vertex] << endl;

		if (query_dis < distance + 1e-5) { // 1e-5 is precision
			sorted_int_vector_binary_insert(labels_to_be_removed[target_v], vertex); // this is not a canonical label
		}
		else {
			incremental_label_vectors[target_v].push_back(*begin);
		}
	}

}

void canonical_repair_element2(int target_v) {


	auto& it = labels_to_be_removed[target_v];
	vector<two_hop_label_v1> new_label_vector;
	auto begin = L_temp_595[target_v].begin(), end = L_temp_595[target_v].end();
	for (; begin != end; begin++) {
		int vertex = begin->vertex;
		if (sorted_int_vector_binary_search(it, vertex) == false) {
			new_label_vector.push_back(*begin);
		}
	}
	L_temp_595[target_v] = new_label_vector;
	vector<two_hop_label_v1>(L_temp_595[target_v]).swap(L_temp_595[target_v]);

}

void canonical_repair_multi_threads(graph_hash_of_mixed_weighted& instance_graph, long long int& label_size_before_canonical_repair,
	long long int& label_size_after_canonical_repair, double& canonical_repair_remove_label_ratio, int num_of_threads) {

	int max_N_ID = L_temp_595.size();
	incremental_label_vectors.resize(max_N_ID);
	labels_to_be_removed.resize(max_N_ID); // this is like a sorted directed adjacenct list; labels_to_be_removed[i].contain(j) means L[i].(label.vertex == j) should be removed

	ThreadPool pool(num_of_threads);
	std::vector< std::future<int> > results; // return typename: xxx

	auto* instance_graph_pointer = &instance_graph;

	/*find labels_to_be_removed*/
	for (int target_v = 0; target_v < max_N_ID; target_v++) {
		int size = L_temp_595[target_v].size();
		if (size > 0) {		
			results.emplace_back(
				pool.enqueue([instance_graph_pointer, target_v] { // pass const type value j to thread; [] can be empty
					canonical_repair_element1(instance_graph_pointer, target_v);
					return 1; // return to results; the return type must be the same with results
					})
			);
		}
	}
	for (auto&& result : results)
		result.get(); // all threads finish here
	results.clear();


	/*remove labels_to_be_removed*/
	label_size_before_canonical_repair = 0;
	label_size_after_canonical_repair = 0;
	for (int target_v = 0; target_v < max_N_ID; target_v++) {
		label_size_before_canonical_repair = label_size_before_canonical_repair + L_temp_595[target_v].size();
		label_size_after_canonical_repair = label_size_after_canonical_repair + L_temp_595[target_v].size();
		int size = labels_to_be_removed[target_v].size();		
		if (size > 0) {
			label_size_after_canonical_repair = label_size_after_canonical_repair - size;
			results.emplace_back(
				pool.enqueue([target_v] { // pass const type value j to thread; [] can be empty
						canonical_repair_element2(target_v);
						return 1; // return to results; the return type must be the same with results
					})
			);
		}
	}
	for (auto&& result : results)
		result.get(); // all threads finish here

	canonical_repair_remove_label_ratio = (double)(label_size_before_canonical_repair - label_size_after_canonical_repair) / label_size_before_canonical_repair;
}











/*codes for querying distances or paths*/

/*query dis and path after 2019R1*/

float graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(vector<vector<two_hop_label_v1>>& L, int source, int terminal) {

	/*return std::numeric_limits<float>::max() is not connected*/

	if (source == terminal) {
		return 0;
	}

	float distance = std::numeric_limits<float>::max(); // if disconnected, return this large value
	auto vector1_check_pointer = L[source].begin();
	auto vector2_check_pointer = L[terminal].begin();
	auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
	while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end) {
		if (vector1_check_pointer->vertex == vector2_check_pointer->vertex) {
			float dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
			if (distance > dis) {
				distance = dis;
			}
			vector1_check_pointer++;
		}
		else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex) {
			vector2_check_pointer++;
		}
		else {
			vector1_check_pointer++;
		}
	}

	return distance;

}

float graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1 /* we assume that source and terminal are not reduced by 2019R1 here*/
(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& f_2019R1, int source, int terminal)
{

	/* we assume that source and terminal are not reduced by 2019R1*/

	//cout << "source=" << source << endl;
	//cout << "terminal=" << terminal << endl;
	//cout << "reduction_measures_2019R2[source] " << reduction_measures_2019R2[source] << endl;
	//cout << "reduction_measures_2019R2[terminal] " << reduction_measures_2019R2[terminal] << endl;


	if (source == terminal) {
		return 0;
	}

	vector<float> selected_distance = { std::numeric_limits<float>::max() }; // Store all path lengths to be selected; disconnected = std::numeric_limits<float>::max()
	auto s_adj_begin = adjs[source].begin();
	auto s_adj_end = adjs[source].end();
	auto t_adj_begin = adjs[terminal].begin();
	auto t_adj_end = adjs[terminal].end();
	if (reduction_measures_2019R2[source] == 2)
	{
		if (reduction_measures_2019R2[terminal] == 2)
		{
			/*"Both reduced"*/
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
				{
					if (f_2019R1[it1->first] == it1->first && f_2019R1[it2->first] == it2->first) {
						float x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, it2->first);
						if (x == std::numeric_limits<float>::max() && full_two_hop_labels) { // if not full_two_hop_labels, then we may query a max dis while connected
							return x;
						}
						else {
							selected_distance.push_back(x + float(it1->second) + float(it2->second));
						}
					}
				}
			}
		}
		else
		{
			//cout << "here1" << endl;
			/*"Only source reduced"*/
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				//cout << "it1->first: " << it1->first << " it1->second: " << it1->second << endl;
				if (f_2019R1[it1->first] == it1->first) { 
					float x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, terminal);
					if (x == std::numeric_limits<float>::max() && full_two_hop_labels) {
						return x;
					}
					else {
						selected_distance.push_back(x + float(it1->second));
					}
				}
			}
		}
	}
	else
	{
		if (reduction_measures_2019R2[terminal] == 2)
		{
			//cout << "here1" << endl;
			/*"Only terminal reduced"*/
			for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
			{
				//cout << "it2->first: " << it2->first << " it2->second: " << it2->second << endl;
				if (f_2019R1[it2->first] == it2->first) { 
					float x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, source, it2->first);
					//cout << "x " << x << endl;
					if (x == std::numeric_limits<float>::max() && full_two_hop_labels) {
						return x;
					}
					else {
						selected_distance.push_back(x + float(it2->second));
					}
				}
			}
		}
		else
		{
			/*"Nothing happened"*/
			selected_distance.push_back(graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, source, terminal));
		}
	}

	float dis = *min_element(selected_distance.begin(), selected_distance.end());

	return dis;
}

vector<pair<int, int>> graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1 /* we assume that source and terminal are not reduced by 2019R1*/
(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
{

	/* we assume that source and terminal are not reduced by 2019R1*/

	vector<pair<int, int>> paths;
	if (source == terminal)
	{
		return paths;
	}

	float min_dis = std::numeric_limits<float>::max();
	vector<pair<int, int>> partial_edges(2);

	auto s_adj_begin = adjs[source].begin();
	auto s_adj_end = adjs[source].end();
	auto t_adj_begin = adjs[terminal].begin();
	auto t_adj_end = adjs[terminal].end();

	if (reduction_measures_2019R2[source] == 2)
	{
		if (reduction_measures_2019R2[terminal] == 2)
		{
			/*"Both reduced"*/
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
				{
					if (f_2019R1[it1->first] == it1->first && f_2019R1[it2->first] == it2->first) { 
						float x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, it2->first) + float(it1->second) + float(it2->second);
						/*After removing the two edges, it becomes the fourth case: nothing happened*/
						if (min_dis > x)
						{
							min_dis = x;
							partial_edges[0] = { it1->first, source };
							partial_edges[1] = { it2->first, terminal };
						}
					}				
				}
			}

			if (min_dis == std::numeric_limits<float>::max())
			{ // disconnected
				return paths;
			}

			// add partial_edges into paths
			paths.push_back(partial_edges[0]);
			paths.push_back(partial_edges[1]);

			// find new edges
			vector<pair<int, int>> new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, partial_edges[0].first, partial_edges[1].first);
			if (new_edges.size() > 0)
			{
				for (int i = new_edges.size() - 1; i >= 0; i--)
				{
					paths.push_back(new_edges[i]);
				}
			}
		}
		/*"Only source reduced"*/
		else
		{
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				if (f_2019R1[it1->first] == it1->first) { 
					float x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, terminal) + float(it1->second);
					if (min_dis > x)
					{
						min_dis = x;
						partial_edges[0] = { it1->first, source };
					}
				}
			}
			if (min_dis == std::numeric_limits<float>::max())
			{ // disconnected
				return paths;
			}

			// add partial_edges into paths
			paths.push_back(partial_edges[0]);

			// find new edges
			vector<pair<int, int>> new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, partial_edges[0].first, terminal);
			if (new_edges.size() > 0)
			{
				for (int i = new_edges.size() - 1; i >= 0; i--)
				{
					paths.push_back(new_edges[i]);
				}
			}
		}
	}
	else
	{
		/*"Only terminal reduced"*/
		if (reduction_measures_2019R2[terminal] == 2)
		{
			for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
			{
				if (f_2019R1[it2->first] == it2->first) { 
					float x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, source, it2->first) + float(it2->second);
					if (min_dis > x)
					{
						min_dis = x;
						partial_edges[0] = { it2->first, terminal };
					}
				}
			}
			if (min_dis == std::numeric_limits<float>::max())
			{ // disconnected
				return paths;
			}
			// add partial_edges into paths
			paths.push_back(partial_edges[0]);
			// find new edges
			vector<pair<int, int>> new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, partial_edges[0].first);
			if (new_edges.size() > 0)
			{
				for (int i = new_edges.size() - 1; i >= 0; i--)
				{
					paths.push_back(new_edges[i]);
				}
			}
		}
		/*
		"Nothing happened"
		In this case, the problem that the removed vertices appear in the path needs to be solved
		*/
		else
		{
			int vector1_capped_v_parent = 0, vector2_capped_v_parent = 0;
			float distance = std::numeric_limits<float>::max(); // if disconnected, retun this large value
			bool connected = false;
			auto vector1_check_pointer = L[source].begin();
			auto vector2_check_pointer = L[terminal].begin();
			auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
			{
				if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
				{
					connected = true;
					float dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
					if (distance > dis)
					{
						distance = dis;
						vector1_capped_v_parent = vector1_check_pointer->parent_vertex;
						vector2_capped_v_parent = vector2_check_pointer->parent_vertex;
					}
					vector1_check_pointer++;
				}
				else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex)
				{
					vector2_check_pointer++;
				}
				else
				{
					vector1_check_pointer++;
				}
			}

			if (connected)
			{
				/*the following code will not induce redundent edge, since for each vertex v_k, there is a label (v_k,0,v_k) in L[v_k];
				you must ascending from both source and terminal, otherwise you may not be able to extract SP */
				if (source != vector1_capped_v_parent)
				{
					paths.push_back({ source, vector1_capped_v_parent });
					source = vector1_capped_v_parent; // ascending from source
				}
				if (terminal != vector2_capped_v_parent)
				{
					paths.push_back({ terminal, vector2_capped_v_parent });
					terminal = vector2_capped_v_parent; // ascending from terminal
				}
			}
			else
			{
				return paths;
			}

			// find new edges
			vector<pair<int, int>> new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, terminal);

			if (new_edges.size() > 0)
			{
				for (int i = new_edges.size() - 1; i >= 0; i--)
				{
					paths.push_back(new_edges[i]);
				}
			}
		}
	}

	return paths;
}

/*query dis (final)*/
float graph_hash_of_mixed_weighted_two_hop_v1_extract_distance
(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& reduction_measures_2019R1, 
	vector<int>& f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
{

	if (source == terminal) {
		return 0;
	}

	if (reduction_measures_2019R1[source] == 11)
	{
		/* case 2 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			if (f_2019R1[source] == f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = min_adjs[source];
				return float(s_min_adj.second * 2);
			}
			else
			{
				return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, f_2019R1[source], f_2019R1[terminal]);
			}
		}
		/* case 3 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, f_2019R1[source], f_2019R1[terminal]);
		}
		else { /* case 1 */
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, f_2019R1[source], terminal);
		}	
	}
	else if (reduction_measures_2019R1[source] == 12)
	{
		/* case 5 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, f_2019R1[source], f_2019R1[terminal]);
		}
		/* case 6 -- same with case 3 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			if (f_2019R1[source] == f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = min_adjs[source];
				float s_t_weight = float(graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, terminal));
				return float(s_min_adj.second * 2) > s_t_weight ? s_t_weight : float(s_min_adj.second * 2);
			}
			else
			{
				return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, f_2019R1[source], f_2019R1[terminal]);
			}
		}
		/* case 4 */
		else {
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, f_2019R1[source], terminal);
		}	
	}
	else {
		/* case 8 -- same with case 1 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, source, f_2019R1[terminal]);
		}
		/* case 9 -- same with case 4 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			//cout << "here" << endl;
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, source, f_2019R1[terminal]);
		}
		else {
			//cout << "here" << endl;
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, source, terminal);
		}	
	}

}

/*query paths (final)*/
void graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(vector<pair<int, int>>& paths, int old_, int new_)
{
	/*
		this function can change specific vector on paths, used for shortesd paths query
		take care, we assume that the 'old_' vector only appear once
	*/
	for (auto it = paths.begin(); it != paths.end(); it++)
	{
		if (it->first == old_)
		{
			it->first = new_;
			break;
		}
		if (it->second == old_)
		{
			it->second = new_;
			break;
		}
	}
}

vector<pair<int, int>> graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path
(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& reduction_measures_2019R1, 
	vector<int>& f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
{

	vector<pair<int, int>> paths;
	if (source == terminal)
	{
		return paths;
	}
	
	if (reduction_measures_2019R1[source] == 11)
	{
		/* case 2 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			if (f_2019R1[source] == f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = min_adjs[source];
				paths.push_back({ source, s_min_adj.first });
				paths.push_back({ terminal, s_min_adj.first });
				return paths;
			}
			else
			{
				paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
				graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
				graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
				return paths;
			}
		}
		/* case 3 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
			return paths;
		}
		/* case 1 */
		else {
			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], terminal);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
			return paths;
		}	
	}
	else if (reduction_measures_2019R1[source] == 12)
	{
		/* case 5 -- same with case 3 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
			return paths;
		}
		/* case 6 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			if (f_2019R1[source] == f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = min_adjs[source];
				float s_t_weight = float(graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, terminal));
				if (float(s_min_adj.second * 2) < s_t_weight)
				{
					paths.push_back({ source, s_min_adj.first });
					paths.push_back({ terminal, s_min_adj.first });
					return paths;
				}
				else
				{
					paths.push_back({ source, terminal });
					return paths;
				}
			}
			else // f_2019R1[source] != f_2019R1[terminal]
			{
				paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
				graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
				graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
				return paths;
			}
		}
		/* case 4 */
		else {
			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], terminal);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
			return paths;
		}
		
	}
	else {
		/* case 8 -- same with case 1 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, f_2019R1[terminal]);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
			return paths;
		}
		/* case 9 -- same with case 4 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, f_2019R1[terminal]);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
			return paths;
		}
		/* case 7 */
		else {
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, terminal);
		}
		
	}

}







/*query predecessors in instance_graph*/

pair<int, int> graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1 /* we assume that source and terminal are not reduced by 2019R1*/
(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
{

	/* we assume that source and terminal are not reduced by 2019R1*/
	
	if (source == terminal)
	{
		return { source, terminal }; // source_predecessor and terminal_predecessor;
	}

	float min_dis = std::numeric_limits<float>::max();
	vector<pair<int, int>> partial_edges(2);

	auto s_adj_begin = adjs[source].begin();
	auto s_adj_end = adjs[source].end();
	auto t_adj_begin = adjs[terminal].begin();
	auto t_adj_end = adjs[terminal].end();

	if (reduction_measures_2019R2[source] == 2)
	{
		if (reduction_measures_2019R2[terminal] == 2)
		{
			/*"Both reduced"*/
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
				{
					if (f_2019R1[it1->first] == it1->first && f_2019R1[it2->first] == it2->first) {
						float x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, it2->first) + float(it1->second) + float(it2->second);
						/*After removing the two edges, it becomes the fourth case: nothing happened*/
						if (min_dis > x)
						{
							min_dis = x;
							partial_edges[0] = { it1->first, source };
							partial_edges[1] = { it2->first, terminal };
						}
					}
				}
			}

			if (min_dis == std::numeric_limits<float>::max())
			{ // disconnected
				return { source, terminal };
			}

			return { partial_edges[0].first, partial_edges[1].first };
		}
		/*"Only source reduced"*/
		else
		{
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				if (f_2019R1[it1->first] == it1->first) {
					float x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, terminal) + float(it1->second);
					if (min_dis > x)
					{
						min_dis = x;
						partial_edges[0] = { it1->first, source };
					}
				}
			}
			if (min_dis == std::numeric_limits<float>::max())
			{ // disconnected
				return { source, terminal };
			}

			pair<int, int> recursive_two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, partial_edges[0].first, terminal);

			return { partial_edges[0].first , recursive_two_predecessors.second };
		}
	}
	else
	{
		/*"Only terminal reduced"*/
		if (reduction_measures_2019R2[terminal] == 2)
		{
			for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
			{
				if (f_2019R1[it2->first] == it2->first) {
					float x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, source, it2->first) + float(it2->second);
					if (min_dis > x)
					{
						min_dis = x;
						partial_edges[0] = { it2->first, terminal };
					}
				}
			}

			if (min_dis == std::numeric_limits<float>::max())
			{ // disconnected
				return { source, terminal };
			}

			pair<int, int> recursive_two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, partial_edges[0].first);

			return { recursive_two_predecessors.first, partial_edges[0].first };
		}
		/*
		"Nothing happened"
		In this case, the problem that the removed vertices appear in the path needs to be solved
		*/
		else
		{
			int vector1_capped_v_parent = 0, vector2_capped_v_parent = 0;
			float distance = std::numeric_limits<float>::max(); // if disconnected, return this large value
			bool connected = false;
			auto vector1_check_pointer = L[source].begin();
			auto vector2_check_pointer = L[terminal].begin();
			auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
			{
				if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
				{
					connected = true;
					float dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
					if (distance > dis)
					{
						distance = dis;
						vector1_capped_v_parent = vector1_check_pointer->parent_vertex;
						vector2_capped_v_parent = vector2_check_pointer->parent_vertex;
					}
					vector1_check_pointer++;
				}
				else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex)
				{
					vector2_check_pointer++;
				}
				else
				{
					vector1_check_pointer++;
				}
			}

			if (connected)
			{
				/*the following code will not induce redundent edge, since for each vertex v_k, there is a label (v_k,0,v_k) in L[v_k];
				you must ascending from both source and terminal, otherwise you may not be able to extract SP */
				pair<int, int> recursive_two_predecessors;
				if (source != vector1_capped_v_parent) // ascending from source
				{
					recursive_two_predecessors.first = vector1_capped_v_parent;
				}
				else {
					recursive_two_predecessors.first = source;
				}
				if (terminal != vector2_capped_v_parent) // ascending from terminal
				{
					recursive_two_predecessors.second = vector2_capped_v_parent;
				}
				else {
					recursive_two_predecessors.second = terminal;
				}
				return recursive_two_predecessors;
			}
			else
			{
				return { source, terminal };
			}
		}
	}
}

pair<int, int> graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors
(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& reduction_measures_2019R1,
	vector<int>& f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
{

	pair<int, int> two_predecessors; // source_predecessor and terminal_predecessor;  
	/*
	if (source_predecessor == source && terminal_predecessor == terminal){
	no edge is in path
	}

	if (source_predecessor == terminal && terminal_predecessor == source){
	one edge is in path
	}

	*/

	if (source == terminal)
	{
		return { source, terminal };
	}

	if (reduction_measures_2019R1[source] == 11)
	{
		/* case 2 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			if (f_2019R1[source] == f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = min_adjs[source];
				return { s_min_adj.first , s_min_adj.first };
			}
			else
			{
				two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
				if (two_predecessors.first == f_2019R1[source]) {
					two_predecessors.first = source;
				}
				if (two_predecessors.second == f_2019R1[source]) {
					two_predecessors.second = source;
				}
				if (two_predecessors.first == f_2019R1[terminal]) {
					two_predecessors.first = terminal;
				}
				if (two_predecessors.second == f_2019R1[terminal]) {
					two_predecessors.second = terminal;
				}
				return two_predecessors;
			}
		}
		/* case 3 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
			if (two_predecessors.first == f_2019R1[source]) {
				two_predecessors.first = source;
			}
			if (two_predecessors.second == f_2019R1[source]) {
				two_predecessors.second = source;
			}
			if (two_predecessors.first == f_2019R1[terminal]) {
				two_predecessors.first = terminal;
			}
			if (two_predecessors.second == f_2019R1[terminal]) {
				two_predecessors.second = terminal;
			}
			return two_predecessors;
		}
		/* case 1 */
		else {
			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], terminal);
			if (two_predecessors.first == f_2019R1[source]) {
				two_predecessors.first = source;
			}
			if (two_predecessors.second == f_2019R1[source]) {
				two_predecessors.second = source;
			}
			return two_predecessors;
		}
	}
	else if (reduction_measures_2019R1[source] == 12)
	{
		/* case 5 -- same with case 3 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
			if (two_predecessors.first == f_2019R1[source]) {
				two_predecessors.first = source;
			}
			if (two_predecessors.second == f_2019R1[source]) {
				two_predecessors.second = source;
			}
			if (two_predecessors.first == f_2019R1[terminal]) {
				two_predecessors.first = terminal;
			}
			if (two_predecessors.second == f_2019R1[terminal]) {
				two_predecessors.second = terminal;
			}
			return two_predecessors;
		}
		/* case 6 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			if (f_2019R1[source] == f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = min_adjs[source];
				float s_t_weight = float(graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, terminal));
				if (float(s_min_adj.second * 2) < s_t_weight)
				{
					return { s_min_adj.first , s_min_adj.first };
				}
				else
				{
					return { terminal, source };
				}
			}
			else // f_2019R1[source] != f_2019R1[terminal]
			{
				two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
				if (two_predecessors.first == f_2019R1[source]) {
					two_predecessors.first = source;
				}
				if (two_predecessors.second == f_2019R1[source]) {
					two_predecessors.second = source;
				}
				if (two_predecessors.first == f_2019R1[terminal]) {
					two_predecessors.first = terminal;
				}
				if (two_predecessors.second == f_2019R1[terminal]) {
					two_predecessors.second = terminal;
				}
				return two_predecessors;
			}
		}
		/* case 4 */
		else {
			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], terminal);
			if (two_predecessors.first == f_2019R1[source]) {
				two_predecessors.first = source;
			}
			if (two_predecessors.second == f_2019R1[source]) {
				two_predecessors.second = source;
			}
			return two_predecessors;
		}
	}
	else {
		/* case 8 -- same with case 1 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, f_2019R1[terminal]);
			if (two_predecessors.first == f_2019R1[terminal]) {
				two_predecessors.first = terminal;
			}
			if (two_predecessors.second == f_2019R1[terminal]) {
				two_predecessors.second = terminal;
			}
			return two_predecessors;
		}
		/* case 9 -- same with case 4 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, f_2019R1[terminal]);
			if (two_predecessors.first == f_2019R1[terminal]) {
				two_predecessors.first = terminal;
			}
			if (two_predecessors.second == f_2019R1[terminal]) {
				two_predecessors.second = terminal;
			}
			return two_predecessors;
		}
		/* case 7 */
		else {
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, terminal);
		}

	}

}