#pragma once

/*the following codes are for testing

---------------------------------------------------
a cpp file (try.cpp) for running the following test code:
----------------------------------------

#include <iostream>
#include <fstream>
using namespace std;
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

#include <build_in_progress/HL/dgraph/exp_dgraph.h>

int main()
{
    cout << "Start running..." << endl;
    auto begin = std::chrono::high_resolution_clock::now();

    main_exp();

    double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
    cout << "END    runningtime: " << runningtime << "s" << endl;
}

------------------------------------------------------------------------------------------
Commends for running the above cpp file on Linux:

g++ -std=c++17 -I/home/boost_1_75_0 -I/home/dengs/dgraph try.cpp -lpthread -Ofast -o A
./A
rm A

(optional to put the above commends in run.sh, and then use the comment: sh run.sh)
*/

#include <unordered_map>
#include <dgraph_v_of_v/dgraph_v_of_v.h>
#include <dgraph_v_of_v/dgraph_generate_random_dgraph.h>
#include <dgraph_v_of_v/dgraph_save_dgraph.h>
#include <dgraph_v_of_v/dgraph_read_dgraph.h>
#include <build_in_progress/HL/dgraph/dgraph_change_IDs.h>
#include <dgraph_v_of_v/dgraph_shortest_paths.h>
#include <build_in_progress/HL/dgraph/dgraph_two_hop_label.h>
#include <build_in_progress/HL/dgraph/dgraph_PLL.h>
#include <build_in_progress/HL/dgraph/dgraph_PSL.h>
#include <build_in_progress/HL/dgraph/dgraph_CT.h>
#include <dgraph_v_of_v/dgraph_binary_save_read_dgraph.h>


double querying_element(dgraph_case_info_v2* case_info, vector<pair<int, int>>* query_list) {

	try {
        auto& list = *query_list;
		auto begin = std::chrono::high_resolution_clock::now();
        for (int i = list.size() - 1; i >= 0; i--) {
            CT_extract_distance(*case_info, list[i].first, list[i].second);
        }
        return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1e9; // s
	}
	catch (...) {
		cout << "querying_element throw error!" << endl;
        exit(1);
	}
}

double querying(dgraph_case_info_v2& case_info, vector<pair<int, int>>& query_list) {

	int query_times = query_list.size();
    auto* case_info_p = &case_info;
	ThreadPool pool(80); // thread_num for querying
	std::vector< std::future<double> > results; // return typename: xxx
    vector<vector<pair<int, int>>> query_lists(80);
    for (int i = 0; i < query_times; i++) {
        int j = i % 80;
        query_lists[j].push_back(query_list[j]);
    }
	for (int i = 0; i < 80; i++) {
        auto* it = &query_lists[i];
		results.emplace_back(
			pool.enqueue([case_info_p, it] { // pass const type value j to thread; [] can be empty
				return querying_element(case_info_p, it);
				})
		);
	}
    double query_dis_avg_time = 0;
	for (auto&& result : results) {
		query_dis_avg_time += result.get();
	}
    return query_dis_avg_time / (double)query_times;
}

void exp_element(string data_name, int ec_type, int thread_num, int d, long long int max_bit_size, double max_run_time_seconds, int query_times) {

	cout << "start indexing " << data_name << " " << ec_type << endl;

    string path = "/home/malu/DHL_exp/" + data_name + "/";   //path = "";
	dgraph_v_of_v<two_hop_weight_type> input_graph;
	if (ec_type == 0) {
		dgraph_binary_read_dgraph(path + data_name + "_Jacard.bin", input_graph);
	}
	else if (ec_type == 1) {
		dgraph_binary_read_dgraph(path + data_name + "_random.bin", input_graph);
	}
    else if (ec_type == 2) {
        dgraph_binary_read_dgraph(path + data_name + "_random.bin", input_graph);
        for (int i = input_graph.INs.size() - 1; i >= 0; i--) {
            for (int j = input_graph.INs[i].size() - 1; j >= 0; j--) {
                input_graph.INs[i][j].second = 1;
            }
        }
        for (int i = input_graph.OUTs.size() - 1; i >= 0; i--) {
            for (int j = input_graph.OUTs[i].size() - 1; j >= 0; j--) {
                input_graph.OUTs[i][j].second = 1;
            }
        }
    }
	vector<pair<int, int>> query_list;
    boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(input_graph.INs.size() - 1) };
    for (int i = 0; i < query_times; i++) {
        query_list.push_back({ dist(boost_random_time_seed) ,dist(boost_random_time_seed) });
    }

	/*output*/
	ofstream outputFile;
	outputFile.precision(8);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);
	string save_name = "exp_" + data_name + "_ectype_" + to_string(ec_type) + "_T_" + to_string(thread_num) + "_d_" + to_string(d);

	outputFile.open(save_name + ".csv");
	outputFile << "data,ec_type,thread_num,d,"
		<< "CT-DPLL_MB,CT-DPLL_time,CT-DPLL_query_dis_time,"
        << "CT-DPSL_MB,CT-DPSL_time,CT-DPSL_query_dis_time" << endl;

	outputFile << data_name << "," << ec_type << "," << thread_num << "," << d << "," << std::flush;

    double DPLL_time = 0;

	/* CT-DPLL */
	if (1) {
        string algo_name = "CT-DPLL";
		cout << "start " + algo_name << endl;
        dgraph_case_info_v2 case_info;
        case_info.use_PLL = 1;
        case_info.two_hop_case_info.use_canonical_repair = 1;
        case_info.two_hop_order_method = 0; // unweighted degree sort
        case_info.d = d;
        case_info.thread_num = thread_num;
        case_info.max_bit_size = max_bit_size;
        case_info.max_run_time_seconds = max_run_time_seconds;

        int catch_error = 0;
        try {
            CT_dgraph(input_graph, case_info);
        }
        catch (string s) {
            cout << "ERROR: should use a dataset that CT-DPLL can conquer!" << endl;
            exit(1);
            if (s == reach_limit_error_string_MB) {
                catch_error = 1;
            }
            else if (s == reach_limit_error_string_time) {
                catch_error = 2;
            }
            clear_gloval_values_CT();
        }
		if (catch_error) {
			if (catch_error == 1) {
				outputFile << "0,0,0," << std::flush;
			}
			else if (catch_error == 2) {
				outputFile << "-1,-1,-1," << std::flush;
			}
		}
		else {
            cout << "start querying" << endl;
            DPLL_time = case_info.time_total;
			double query_avg_time = querying(case_info, query_list);
			outputFile << (double)case_info.total_label_bit_size() / 1024 / 1024 << "," << case_info.time_total - case_info.time5_core_indexs_prepare4 << "," << query_avg_time << "," << std::flush;
		}
		case_info.record_all_details(save_name + "_" + algo_name);
        case_info.clear_labels();
	}

    /* CT-DPSL */
    if (1) {
        string algo_name = "CT-DPSL";
        cout << "start " + algo_name << endl;
        dgraph_case_info_v2 case_info;
        case_info.use_PLL = 0;
        case_info.two_hop_case_info.use_canonical_repair = 1;
        case_info.two_hop_order_method = 0; // unweighted degree sort
        case_info.d = d;
        case_info.thread_num = thread_num;
        case_info.max_bit_size = max_bit_size;
        case_info.max_run_time_seconds = DPLL_time * 2;

        int catch_error = 0;
        try {
            CT_dgraph(input_graph, case_info);
        }
        catch (string s) {
            if (s == reach_limit_error_string_MB) {
                catch_error = 1;
            }
            else if (s == reach_limit_error_string_time) {
                catch_error = 2;
            }
            clear_gloval_values_CT();
        } 
        if (catch_error) {
            if (catch_error == 1) {
                outputFile << "0,0,0," << std::flush;
            }
            else if (catch_error == 2) {
                outputFile << "-1,-1,-1," << std::flush;
            }
        }
        else {
            cout << "start querying" << endl;
            double query_avg_time = querying(case_info, query_list);
            outputFile << (double)case_info.total_label_bit_size() / 1024 / 1024 << "," << case_info.time_total << "," << query_avg_time << "," << std::flush;
            case_info.record_all_details(save_name + "_" + algo_name + "_success");
        }
        case_info.record_all_details(save_name + "_" + algo_name);
        case_info.clear_labels();
    }

    outputFile << endl;
}

void main_exp() {

	vector<string> used_datas = { "soc-Epinions1","citeseer", "web-Stanford","lasagne-yahoo","digg-friends","amazon0601","youtube-links",
            "prosper-loans","higgs-twitter-social","zhishi-baidu-internallink" };

	long long int max_bit_size = pow(1024, 3) * 500;
	double max_run_time_seconds = 3600 * 1;
    int query_times = 1e4;

	/*Jacard & random*/
	if (1) {
		int thread_num = 80;
		int d = 0; // querying is too slow when d is large
		for (int i = 0; i < used_datas.size(); i++) {
			exp_element(used_datas[i], 0, thread_num, d, max_bit_size, max_run_time_seconds, query_times);
			exp_element(used_datas[i], 1, thread_num, d, max_bit_size, max_run_time_seconds, query_times);
            exp_element(used_datas[i], 2, thread_num, d, max_bit_size, max_run_time_seconds, query_times);
		}
	}

	/*Jacard & random*/
	if (0) {
		int thread_num = 50;
		int d = 20; // querying is too slow when d is large
		for (int i = 0; i < used_datas.size(); i++) {
			exp_element(used_datas[i], 0, thread_num, d, max_bit_size, max_run_time_seconds, query_times);
			exp_element(used_datas[i], 1, thread_num, d, max_bit_size, max_run_time_seconds, query_times);
		}
	}

	/*Jacard & random*/
	if (0) {
		int thread_num = 80;
		int d = 50; // querying is too slow when d is large
		for (int i = 0; i < used_datas.size(); i++) {
			exp_element(used_datas[i], 0, thread_num, d, max_bit_size, max_run_time_seconds, query_times);
			exp_element(used_datas[i], 1, thread_num, d, max_bit_size, max_run_time_seconds, query_times);
		}
	}

}

