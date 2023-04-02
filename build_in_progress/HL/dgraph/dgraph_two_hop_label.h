#pragma once
#include <fstream>
#include <string>
#include <dgraph_v_of_v/dgraph_v_of_v.h>

/*
    to include this file, please copy
    #include <dgraph_v_of_v/dgraph_two_hop_label.h>
*/

#define two_hop_weight_type float

class two_hop_label
{
  public:
    int vertex;
    two_hop_weight_type distance;
};

/*common functions shared by PLL and PSL*/
bool compare_two_hop_label_vertex_small_to_large_and_distance_small_to_large(const two_hop_label &i, two_hop_label &j)
{
    return i.vertex == j.vertex? i.distance<j.distance : i.vertex < j.vertex;
}

bool compare_two_hop_label_vertex_small_to_large(const two_hop_label &i, two_hop_label &j)
{
    return i.vertex < j.vertex;
}

//定义全局变量
/* two hop label and its parameters */
vector<vector<two_hop_label>> L_595;
vector<vector<two_hop_label>> L_temp_in;
vector<vector<two_hop_label>> L_temp_out;
long long int max_labal_size_595;
long long int labal_size_595;
/* multi thread parameters */
std::shared_mutex mtx_595_1, mtx_595_2;
int max_N_ID_for_mtx_595 = 1e7;
vector<std::shared_mutex> mtx_595(max_N_ID_for_mtx_595);
/* dijkstra parameters */
dgraph_v_of_v<two_hop_weight_type> ideal_dgraph;
dgraph_v_of_v<two_hop_weight_type> revse_dgraph;
queue<int> Qid_595;
vector<vector<two_hop_weight_type>> dij_dist;
vector<vector<two_hop_weight_type>> L_vk_tmp;
//special for PSL          0 out    1 in
//vector<int> endpos1[2];
vector<int> increment[2];
vector<int> pos_595[2];
vector<int> pos_2_595[2];
queue<int> thread_id;
shared_mutex mtx;
vector<std::vector<int>> dirt;
vector<std::vector<two_hop_weight_type>> dmin;
vector<vector<two_hop_label>> L_PSL[2];
vector<vector<two_hop_label>> L_PSL_temp[2];
//special for canonical repair
vector<vector<two_hop_label>> L_temp_in_for_canonical_repair;
vector<vector<two_hop_label>> L_temp_out_for_canonical_repair;

//清除全局变量
void dgraph_clear_global_values_PLL()
{
    vector<vector<two_hop_label>>().swap(L_595);
    vector<vector<two_hop_label>>().swap(L_temp_in);
    vector<vector<two_hop_label>>().swap(L_temp_out);
    vector<vector<two_hop_weight_type>>().swap(dij_dist);
    vector<vector<two_hop_weight_type>>().swap(L_vk_tmp);
}

void dgraph_clear_global_values_in_canonical_repair()
{
    vector<vector<two_hop_label>>().swap(L_temp_in_for_canonical_repair);
    vector<vector<two_hop_label>>().swap(L_temp_out_for_canonical_repair);
}

//清除全局变量
void dgraph_clear_global_values_PSL()
{   
    for (int k=0;k<2;k++)
    {
        vector<int>().swap(increment[k]);
        vector<int>().swap(pos_595[k]);
        vector<int>().swap(pos_2_595[k]);
        vector<vector<two_hop_label>>().swap(L_PSL[k]);
        vector<vector<two_hop_label>>().swap(L_PSL_temp[k]);
    }
    vector<vector<int>>().swap(dirt);
    vector<std::vector<two_hop_weight_type>>().swap(dmin);
}

void label_output_to_file(std::string instance_name, vector<vector<two_hop_label>> &L_in, vector<vector<two_hop_label>> &L_out)
{
    std::ofstream outputFile;
    outputFile.open(instance_name);

    outputFile << "in label" << std::endl;

    int size1 = L_in.size();
    for (int i = 0; i < size1; i++)
    {
        outputFile << i << ":\t";
        int size2 = L_in[i].size();
        for (int j = 0; j < size2; j++)
        {
            outputFile << "(" << L_in[i][j].vertex << "," << L_in[i][j].distance << ") ";
        }
        outputFile << std::endl;
    }

    outputFile << std::endl << std::endl << std::endl << "out label" << std::endl;;

    size1 = L_out.size();
    for (int i = 0; i < size1; i++)
    {
        outputFile << i << ":\t";
        int size2 = L_out[i].size();
        for (int j = 0; j < size2; j++)
        {
            outputFile << "(" << L_out[i][j].vertex << "," << L_out[i][j].distance << ") ";
        }
        outputFile << std::endl;
    }
}

class dgraph_case_info_v1
{
public:
    //bool use_2019R1 = false;
	//bool use_2019R2 = false;

    //double time_2019R1 = 0;
	//double time_2019R2= 0;

	//double time_initialization = 0;
	//double time_generate_labels = 0;
	
	/*running limits*/
	long long int max_labal_size = 1e12; // 2-hop-label num
	double max_run_time_seconds = 1e12;
    bool use_canonical_repair = false;

    //for canonical repair
    long long int label_size_before_canonical_repair = 0;
    long long int label_size_after_canonical_repair = 0;

    /*labels*/
	//vector<int> reduction_measures_2019R2; // for 2019 R2
	//vector<int> reduction_measures_2019R1; // for 2019 R1;  11 means equivalent_1 relation (no edge between), 12 means equivalent_2 relation (edge between)
	//vector<int> f_2019R1; // for 2019 R1
    vector<vector<two_hop_label>> L_in;
    vector<vector<two_hop_label>> L_out;

    /*compute label size*/
	long long int compute_label_bit_size() {
		long long int size = 0;
		//size = size + reduction_measures_2019R2.size() * 4;
		//size = size + reduction_measures_2019R1.size() * 4;
		//size = size + f_2019R1.size() * 4;
        for (auto it = L_in.begin(); it != L_in.end(); it++)
        {
			size = size + (*it).size() * sizeof(two_hop_label); 
		}

        for (auto it = L_out.begin(); it != L_out.end(); it++)
        {
            size = size + (*it).size() * sizeof(two_hop_label);
        }
		return size;
	}


	/*clear labels*/
	void clear_labels() {
		//vector<int>().swap(reduction_measures_2019R2);
		//vector<int>().swap(reduction_measures_2019R1);
		//vector<int>().swap(f_2019R1);
		vector<vector<two_hop_label>>().swap(L_in);
        vector<vector<two_hop_label>>().swap(L_out);
	}

	/*printing*/
	void print_L() {
		cout << "in label" << std::endl;

        int size1 = L_in.size();
        for (int i = 0; i < size1; i++)
        {
            cout << i << ":\t";
            int size2 = L_in[i].size();
            for (int j = 0; j < size2; j++)
            {
                cout << "(" << L_in[i][j].vertex << "," << L_in[i][j].distance  << ") ";
            }
            cout << std::endl;
        }

        cout << std::endl << std::endl << std::endl << "out label" << std::endl;;

        size1 = L_out.size();
        for (int i = 0; i < size1; i++)
        {
            cout << i << ":\t";
            int size2 = L_out[i].size();
            for (int j = 0; j < size2; j++)
            {
                cout << "(" << L_out[i][j].vertex << "," << L_out[i][j].distance  << ") ";
            }
            cout << std::endl;
        }
	}
};

two_hop_weight_type dgraph_v1_extract_shortest_distance(vector<vector<two_hop_label>> &L_in,
                                                        vector<vector<two_hop_label>> &L_out,
                                           dgraph_v_of_v<two_hop_weight_type> &instance_graph, int source, int terminal)
{
    if (source == terminal)
    {
        return 0;
    }

    two_hop_weight_type distance =
        std::numeric_limits<two_hop_weight_type>::max(); // if disconnected, return this large value
    auto vector1_check_pointer = L_out[source].begin();
    auto vector2_check_pointer = L_in[terminal].begin();
    auto pointer_L_s_end = L_out[source].end(), pointer_L_t_end = L_in[terminal].end();
    while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
    {
        if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
        {
            two_hop_weight_type dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
            if (distance > dis)
            {
                distance = dis;
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

    return distance;
}

/*
vector<pair<int, int>> dgraph_v1_extract_shortest_path(vector<vector<two_hop_label>> &L_in,
                                                       vector<vector<two_hop_label>> &L_out,
                                       dgraph_v_of_v<double> &instance_graph, int source, int terminal)
{   
    vector<pair<int, int>> paths;
    if (source == terminal)
        return paths;

    int vector1_capped_v_parent = 0, vector2_capped_v_parent = 0;
    double distance = std::numeric_limits<double>::max(); // if disconnected, retun this large value
    bool connected = false;
    auto vector1_check_pointer = L_out[source].begin();
    auto vector2_check_pointer = L_in[terminal].begin();
    auto pointer_L_s_end = L_out[source].end(), pointer_L_t_end = L_in[terminal].end();
    while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
    {
        if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
        {
            connected = true;
            double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
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
        if (source != vector1_capped_v_parent)
        {
            paths.push_back({source, vector1_capped_v_parent});
            source = vector1_capped_v_parent; // ascending from source
        }
        if (terminal != vector2_capped_v_parent)
        {
            paths.push_back({vector2_capped_v_parent, terminal});
            terminal = vector2_capped_v_parent; // ascending from terminal
        }
    }
    else
    {
        return paths;
    }

    // find new edges
    vector<pair<int, int>> new_edges =  dgraph_v1_extract_shortest_path(L_in,L_out, instance_graph, source, terminal);

    if (new_edges.size() > 0)
    {
        for (int i = new_edges.size() - 1; i >= 0; i--)
        {
            paths.push_back(new_edges[i]);
        }
    }

    return paths;
}

*/

double dgraph_v1_extract_shortest_distance_for_canonical_repair(int source,int terminal,int in_or_out,int limit)
{
    if (source == terminal)
    {
        return 0;
    }

    two_hop_weight_type distance =
        std::numeric_limits<two_hop_weight_type>::max(); // if disconnected, return this large value
    auto vector1_check_pointer = L_temp_out[source].begin();
    auto vector2_check_pointer = L_temp_in[terminal].begin();
    auto pointer_L_s_end = L_temp_out[source].end(), pointer_L_t_end = L_temp_in[terminal].end();

    while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
    {
        if (in_or_out==1)
        {
            if(vector2_check_pointer->vertex>=limit)
                break;
        }
        else
        {
            if(vector1_check_pointer->vertex>=limit)
                break;
        }


        if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
        {
            two_hop_weight_type dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
            if (distance > dis)
            {
                distance = dis;
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

    return distance;
}

void canonical_repair_in(int target_v)
{
	auto begin = L_temp_in[target_v].begin(), end = L_temp_in[target_v].end();
	for (; begin != end; begin++) {
		int vertex = begin->vertex;
		if (vertex == target_v) {
            L_temp_in_for_canonical_repair[target_v].push_back(*begin);
			continue;
		}
		double distance = begin->distance;
		double query_dis = dgraph_v1_extract_shortest_distance_for_canonical_repair(vertex,target_v,0,vertex);

		if (query_dis >= distance + 1e-5)
        { // 1e-5 is precision
			L_temp_in_for_canonical_repair[target_v].push_back(*begin);
		}
	}
}

void canonical_repair_out(int target_v)
{
    auto begin = L_temp_out[target_v].begin(), end = L_temp_out[target_v].end();
    for (; begin != end; begin++)
    {
        int vertex = begin->vertex;
        if (vertex == target_v)
        {
            L_temp_out_for_canonical_repair[target_v].push_back(*begin);
            continue;
        }
        double distance = begin->distance;
        double query_dis = dgraph_v1_extract_shortest_distance_for_canonical_repair(target_v, vertex, 1, vertex);

        if (query_dis >= distance + 1e-5)
        { // 1e-5 is precision
            L_temp_out_for_canonical_repair[target_v].push_back(*begin);
        }
    }
}


void canonical_repair_multi_threads(int num_of_threads)
{
    int max_N_ID = L_temp_in.size();
    L_temp_in_for_canonical_repair.resize(max_N_ID,vector<two_hop_label>());
    L_temp_out_for_canonical_repair.resize(max_N_ID, vector<two_hop_label>());
    ThreadPool pool(num_of_threads);
    std::vector<std::future<int>> results; // return typename: xxx
    
    /*find labels_to_be_removed_IN_L_temp_in*/
    for (int target_v = 0; target_v < max_N_ID; target_v++)
    {
        int size = L_temp_in[target_v].size();
        if (size > 0)
        {
            results.emplace_back(
                pool.enqueue([target_v] { // pass const type value j to thread; [] can be empty
                    canonical_repair_in(target_v);
                    return 1; // return to results; the return type must be the same with results
                }));
        }
    }

    for (auto&& result : results)
		result.get(); // all threads finish here
	results.clear();
    L_temp_in = L_temp_in_for_canonical_repair;

    for (int target_v = 0; target_v < max_N_ID; target_v++)
    {
        int size = L_temp_out[target_v].size();
        if (size > 0)
        {
            results.emplace_back(pool.enqueue([target_v] { // pass const type value j to thread; [] can be empty
                canonical_repair_out(target_v);
                return 1; // return to results; the return type must be the same with results
            }));
        }
    }

    for (auto &&result : results)
        result.get(); // all threads finish here
    results.clear();
    L_temp_out = L_temp_out_for_canonical_repair;

    return;
}

void canonical_repair_multi_threads_v2(int num_of_threads)
{
    int max_N_ID = L_temp_in.size();
    L_temp_in_for_canonical_repair.resize(max_N_ID, vector<two_hop_label>());
    L_temp_out_for_canonical_repair.resize(max_N_ID, vector<two_hop_label>());
    ThreadPool pool(num_of_threads);
    std::vector<std::future<int>> results; // return typename: xxx

    /*find labels_to_be_removed_IN_L_temp_in*/
    for (int target_v = 0; target_v < max_N_ID; target_v++)
    {
        int size = L_temp_in[target_v].size();
        if (size > 0)
        {
            results.emplace_back(pool.enqueue([target_v] { // pass const type value j to thread; [] can be empty
                canonical_repair_in(target_v);
                canonical_repair_out(target_v);
                return 1; // return to results; the return type must be the same with results
            }));
        }
    }

    for (auto &&result : results)
        result.get(); // all threads finish here
    results.clear();

    L_temp_in = L_temp_in_for_canonical_repair;
    L_temp_out = L_temp_out_for_canonical_repair;

    return;
}