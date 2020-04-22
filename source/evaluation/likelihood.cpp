#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <stdio.h>
#include <gsl/gsl_cdf.h>

#include "state_merger.h"
#include "evaluate.h"
#include "likelihood.h"
#include "parameters.h"

REGISTER_DEF_DATATYPE(likelihood_data);
REGISTER_DEF_TYPE(likelihoodratio);

bool likelihoodratio::consistent(state_merger *merger, apta_node* left, apta_node* right){
    likelihood_data* l = (likelihood_data*) left->data;
    likelihood_data* r = (likelihood_data*) right->data;

    //if(merger->aut->root == left) {inconsistency_found = true; return false;};
    //if(left->depth != right->depth) {inconsistency_found = true; return false;};

    //if(l->num_final() != 0 && r->num_final() == 0) {inconsistency_found = true; return false;};
    //if(l->num_final() == 0 && r->num_final() != 0) {inconsistency_found = true; return false;};

    return count_driven::consistent(merger, left, right);
};

void likelihoodratio::update_likelihood(double left_count, double right_count, double left_divider, double right_divider){    
    if(left_count != 0.0)
        loglikelihood_orig += (left_count + CORRECTION)  * log((left_count + CORRECTION)  / left_divider);
    if(right_count != 0.0)
        loglikelihood_orig += (right_count + CORRECTION) * log((right_count + CORRECTION) / right_divider);
    if(right_count != 0.0 || left_count != 0.0)
        loglikelihood_merged += (left_count + right_count + 2*CORRECTION) * log((left_count + right_count + 2*CORRECTION) / (left_divider + right_divider));
    if(right_count != 0.0 && left_count != 0.0)
        extra_parameters = extra_parameters + 1;

    //cerr << "update: " << left_count << "/" << left_divider << " " << right_count << "/" << right_divider << endl;
    //cerr << loglikelihood_orig << " " << loglikelihood_merged << " " << extra_parameters << endl;
};

/* Likelihood Ratio (LR), computes an LR-test (used in RTI) and uses the p-value as score and consistency */
void likelihoodratio::update_score(state_merger *merger, apta_node* left, apta_node* right){
    likelihood_data* l = (likelihood_data*) left->data;
    likelihood_data* r = (likelihood_data*) right->data;
    
    CORRECTION = 1.0;

    if(r->num_paths() < STATE_COUNT || l->num_paths() < STATE_COUNT) return;
    
    double left_divider = 1.0;
    double right_divider = 1.0;
    double left_count = 0.0;
    double right_count  = 0.0;

    for(type_num_map::iterator it = l->counts_begin(); it != l->counts_end(); it++){
        int type = (*it).first;
        num_map& nm = (*it).second;
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            int symbol = (*it2).first;
            float left_count = (*it2).second;
            float right_count = r->count(type,symbol);
            
            if(left_count >= SYMBOL_COUNT || right_count >= SYMBOL_COUNT){
                left_divider += CORRECTION;
                right_divider += CORRECTION;
            }
        }
    }

    left_divider += (double)l->num_paths();// + (double)l->pos_final();
    right_divider += (double)r->num_paths();// + (double)r->pos_final();
    
    int l1_pool = 0;
    int r1_pool = 0;
    int l2_pool = 0;
    int r2_pool = 0;
    int matching_right = 0;

    for(type_num_map::iterator it = l->counts_begin(); it != l->counts_end(); it++){
        int type = (*it).first;
        num_map& nm = (*it).second;
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            int symbol = (*it2).first;
            float left_count = (*it2).second;
            float right_count = r->count(type,symbol);
    
            matching_right += right_count;
        
            if(left_count >= SYMBOL_COUNT && right_count >= SYMBOL_COUNT)
                update_likelihood(left_count, right_count, left_divider, right_divider);

            if(right_count < SYMBOL_COUNT){
                l1_pool += left_count;
                r1_pool += right_count;
            }
            if(left_count < SYMBOL_COUNT) {
                l2_pool += left_count;
                r2_pool += right_count;
            }
        }
    }
    
    r2_pool += r->num_paths() - matching_right;
    
    left_count = l1_pool;
    right_count = r1_pool;
    
    if(left_count >= SYMBOL_COUNT || right_count >= SYMBOL_COUNT)
        update_likelihood(left_count, right_count, left_divider, right_divider);
    
    left_count = l2_pool;
    right_count = r2_pool;
    
    if(left_count >= SYMBOL_COUNT || right_count >= SYMBOL_COUNT)
        update_likelihood(left_count, right_count, left_divider, right_divider);
    
    left_count = (double)l->pos_final();
    right_count = (double)r->pos_final();
    
    //if(right_count >= SYMBOL_COUNT || right_count >= SYMBOL_COUNT)
    //    update_likelihood(left_count, right_count, left_divider, right_divider);
};

bool likelihoodratio::compute_consistency(state_merger *merger, apta_node* left, apta_node* right){
    double test_statistic = 2.0 * (loglikelihood_orig - loglikelihood_merged);
    double p_value = gsl_cdf_chisq_Q (test_statistic, 1.0 + (double)extra_parameters);
    
    //cerr << loglikelihood_orig << " " << loglikelihood_merged << " " << loglikelihood_orig - loglikelihood_merged << " " << extra_parameters << " " << p_value << endl;
    
    //cerr << "CHECK " << CHECK_PARAMETER << " " << (p_value < CHECK_PARAMETER) << endl;
    
    if (p_value < CHECK_PARAMETER) { return false; }

    if (inconsistency_found) return false;

    return true;
};

double likelihoodratio::compute_score(state_merger *merger, apta_node* left, apta_node* right){
    //if (inconsistency_found) return -1;

    double test_statistic = 2.0 * (loglikelihood_orig - loglikelihood_merged);
    double p_value = gsl_cdf_chisq_Q (test_statistic, (double)extra_parameters);

    //cerr << "merge score: " << p_value << " " << loglikelihood_orig << " " << loglikelihood_merged << " " << extra_parameters << " " << p_value << endl;

    return p_value;
};

void likelihoodratio::reset(state_merger *merger){
    inconsistency_found = false;
    loglikelihood_orig = 0;
    loglikelihood_merged = 0;
    extra_parameters = 0;
};


/*void likelihoodratio::print_dot(iostream& output, state_merger* merger){
    alergia::print_dot(output, merger);
};*/
