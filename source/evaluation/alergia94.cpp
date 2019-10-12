#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <stdio.h>
#include <gsl/gsl_cdf.h>

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
    if(alergia94::EVAL_TYPE == 1){ if(l->pos_paths() + l->pos_final() < STATE_COUNT || r->pos_paths() + l->pos_final() < STATE_COUNT) return true; }
    else if(l->pos_paths() < STATE_COUNT || r->pos_paths() < STATE_COUNT) return true;
    
    double left_count  = 0.0;
    double right_count = 0.0;
    
    double left_total  = (double)l->pos_paths();
    double right_total = (double)r->pos_paths();
    
    if(alergia94::EVAL_TYPE == 1){
        left_total  += (double)l->pos_final();
        right_total += (double)r->pos_final();
    }

    for(type_num_map::iterator it = l->trans_counts.begin(); it != l->trans_counts.end(); it++){
        num_map& l_nm = (*it).second;
        num_map& r_nm = r->trans_counts[(*it).first];
        for(num_map::iterator itnm = l_nm.begin(); itnm != l_nm.end(); ++itnm){
            left_count = (*itnm).second;
            right_count = r_nm[(*itnm).first];
        
            if(left_count >= SYMBOL_COUNT && right_count >= SYMBOL_COUNT){
                if(alergia_consistency(right_count, left_count, right_total, left_total) == false){
                    inconsistency_found = true; return false;
                }
            }
        }
    }
    if(alergia94::EVAL_TYPE == 1){
        left_count = l->pos_final();
        right_count = r->pos_final();
        
        if(left_count >= SYMBOL_COUNT && right_count >= SYMBOL_COUNT){
            if(alergia_consistency(right_count, left_count, right_total, left_total) == false){
                inconsistency_found = true; return false;
            }
        }
    }

    return true;
};

/* ALERGIA, consistency based on Hoeffding bound, only uses positive (type=1) data */
bool alergia94::consistent(state_merger *merger, apta_node* left, apta_node* right){
    if(count_driven::consistent(merger, left, right) == false){ inconsistency_found = true; return false; }
    alergia94_data* l = (alergia94_data*) left->data;
    alergia94_data* r = (alergia94_data*) right->data;
    
    return data_consistent(l, r);
};
