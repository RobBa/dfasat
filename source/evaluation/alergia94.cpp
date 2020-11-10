#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <stdio.h>
//#include <gsl/gsl_cdf.h>

#include "state_merger.h"
#include "evaluate.h"
#include "depth-driven.h"
#include "alergia94.h"

#include "parameters.h"

int alergia94::EVAL_TYPE = 1;


REGISTER_DEF_DATATYPE(alergia94_data);
REGISTER_DEF_TYPE(alergia94);

alergia94_data::alergia94_data(){
};

bool alergia94::alergia_consistency(double right_count, double left_count, double right_total, double left_total){
    double bound = (1.0 / sqrt(left_total) + 1.0 / sqrt(right_total));
    bound = bound * sqrt(0.5 * log(2.0 / CHECK_PARAMETER));
    
    double gamma = (left_count / left_total) - (right_count / right_total);
    
    if(gamma > bound) return false;
    if(-gamma > bound) return false;
    
    return true;
};

bool alergia94::data_consistent(alergia94_data* l, alergia94_data* r){
    /* we ignore low frequency states, decided by input parameter STATE_COUNT */
    if(FINAL_PROBABILITIES) if(r->num_paths() + r->num_final() < STATE_COUNT || l->num_paths() + l->num_final() < STATE_COUNT) return true;
    else if(r->num_paths() < STATE_COUNT || l->num_paths() < STATE_COUNT) return true;

    /* we treat type distributions as independent */
    for(type_num_map::iterator it = l->counts_begin(); it != l->counts_end(); it++) {
        int type = it->first;
        num_map &nm = it->second;

        /* computing the dividers (denominator) */
        double left_divider = (double) l->num_paths(type);
        double right_divider = (double) r->num_paths(type);

        if (FINAL_PROBABILITIES) {
            left_divider += (double) r->num_final(type);
            right_divider += (double) r->num_final(type);
        }

        for (num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2) {
            int symbol = it2->first;
            double left_count = it2->second;
            if (left_count == 0) continue;
            double right_count = r->count(type, symbol);

            if (alergia_consistency(right_count, left_count, right_divider, left_divider) == false) {
                inconsistency_found = true;
                return false;
            }
        }
        /* count the final probabilities */
        if (FINAL_PROBABILITIES) {
            double left_count = l->num_final();
            double right_count = r->num_final();

            if (alergia_consistency(right_count, left_count, right_divider, left_divider) == false) {
                inconsistency_found = true;
                return false;
            }
        }
    }

    /* count the possibly missed right values, with 0 left occurrences */
    for(type_num_map::iterator it = r->counts_begin(); it != r->counts_end(); it++) {
        int type = it->first;
        num_map &nm = it->second;

        /* computing the dividers (denominator) */
        double left_divider = (double) l->num_paths(type);
        double right_divider = (double) r->num_paths(type);

        if (FINAL_PROBABILITIES) {
            left_divider += (double) r->num_final(type);
            right_divider += (double) r->num_final(type);
        }

        for (num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2) {
            int symbol = it2->first;
            double right_count = it2->second;
            if (right_count == 0) continue;
            double left_count = r->count(type, symbol);
            if (left_count != 0) continue;

            if (alergia_consistency(right_count, left_count, right_divider, left_divider) == false) {
                inconsistency_found = true;
                return false;
            }
        }
    }

    return true;
};

/* ALERGIA, consistency based on Hoeffding bound */
bool alergia94::consistent(state_merger *merger, apta_node* left, apta_node* right){
    if(count_driven::consistent(merger, left, right) == false){ inconsistency_found = true; return false; }
    alergia94_data* l = (alergia94_data*) left->data;
    alergia94_data* r = (alergia94_data*) right->data;
    
    return data_consistent(l, r);
};
