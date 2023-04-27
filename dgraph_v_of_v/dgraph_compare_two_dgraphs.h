#pragma once
#include <dgraph_v_of_v/dgraph_v_of_v.h>

template <typename weight_type>
bool dgraph_compare_two_dgraphs(dgraph_v_of_v<weight_type>& g1, dgraph_v_of_v<weight_type>& g2) {

    if (g1.INs.size() != g2.INs.size() || g1.OUTs.size() != g2.OUTs.size()) {
        return false;
    }

    for (auto it = g1.INs.begin(), it2 = g2.INs.begin(); it != g1.INs.end(); it++, it2++) {
        if (*it != *it2) {
            return false;
        }
    }
    for (auto it = g1.OUTs.begin(), it2 = g2.OUTs.begin(); it != g1.OUTs.end(); it++, it2++) {
        if (*it != *it2) {
            return false;
        }
    }

    return true;
}


template <typename weight_type>
bool dgraph_compare_two_dgraphs_not_eaxct_same_weight(dgraph_v_of_v<weight_type>& g1, dgraph_v_of_v<weight_type>& g2) {

    if (g1.INs.size() != g2.INs.size() || g1.OUTs.size() != g2.OUTs.size()) {
        cout << "ssss" << endl;
        return false;
    }

    for (auto it = g1.INs.begin(), it2 = g2.INs.begin(); it != g1.INs.end(); it++, it2++) {
        if (it->size() != it2->size()) {
            cout << "1ssss" << endl;
            return false;
        }
        for (auto it3 = it->begin(), it4 = it2->begin(); it3 != it->end(); it3++, it4++) {
            if (it3->first != it4->first || abs(it3->second - it4->second) > 1e-5) {
                cout << "2ssss" << endl;
                return false;
            }
        }
    }
    for (auto it = g1.OUTs.begin(), it2 = g2.OUTs.begin(); it != g1.OUTs.end(); it++, it2++) {
        if (it->size() != it2->size()) {
            cout << "3ssss" << endl;
            return false;
        }
        for (auto it3 = it->begin(), it4 = it2->begin(); it3 != it->end(); it3++, it4++) {
            if (it3->first != it4->first || abs(it3->second - it4->second) > 1e-5) {
                cout << "4ssss" << endl;
                return false;
            }
        }
    }

    return true;
}

