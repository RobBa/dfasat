#define STATS_GO_INLINE

#include <math.h>
#include <map>
#include <gsl/gsl_cdf.h>
#include "utility/stats.hpp"

#include "state_merger.h"
#include "evaluate.h"
#include "likelihood.h"
#include "parameters.h"

REGISTER_DEF_DATATYPE(likelihood_data);
REGISTER_DEF_TYPE(likelihoodratio);

void likelihood_data::initialize() {
    alergia_data::initialize();
}

bool likelihoodratio::consistent(state_merger *merger, apta_node* left, apta_node* right){
    likelihood_data* l = (likelihood_data*) left->data;
    likelihood_data* r = (likelihood_data*) right->data;

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
};

/* Likelihood Ratio (LR), computes an LR-test (used in RTI) and uses the p-value as score and consistency */
void likelihoodratio::update_score(state_merger *merger, apta_node* left, apta_node* right){
    likelihood_data* l = (likelihood_data*) left->data;
    likelihood_data* r = (likelihood_data*) right->data;

    /* we ignore low frequency states, decided by input parameter STATE_COUNT */
    if(r->num_paths() < STATE_COUNT || l->num_paths() < STATE_COUNT) return;

    /* computing the dividers (denominator) */
    double left_divider = 0.0;
    double right_divider = 0.0;
    double left_count = 0.0;
    double right_count  = 0.0;

    /* we pool low frequency counts (sum them up in a separate bin), decided by input parameter SYMBOL_COUNT
     * we create pools 1 and 2 separately for left and right low counts
     * in this way, we can detect differences in distributions even if all counts are low (i.e. [0,0,1,1] vs [1,1,0,0])
     * we add correction only to the divider, the correction for counts is in update_likelihood() */
    double l1_pool = 0.0;
    double r1_pool = 0.0;
    double l2_pool = 0.0;
    double r2_pool = 0.0;

    int matching_right = 0;

    for(type_num_map::iterator it = l->counts_begin(); it != l->counts_end(); it++){
        int type = it->first;
        num_map& nm = it->second;
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            int symbol = it2->first;
            left_count = it2->second;
            if(left_count == 0) continue;

            right_count = r->count(type,symbol);
            matching_right += right_count;

            if(left_count >= SYMBOL_COUNT && right_count >= SYMBOL_COUNT){
                left_divider += left_count + CORRECTION;
                right_divider += right_count + CORRECTION;
            } else {
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
    }
    r2_pool += r->num_paths() - matching_right;

    if(l1_pool >= SYMBOL_COUNT || r1_pool >= SYMBOL_COUNT){
        left_divider += l1_pool + CORRECTION;
        right_divider += r1_pool + CORRECTION;
    }
    if(l2_pool >= SYMBOL_COUNT || r2_pool >= SYMBOL_COUNT){
        left_divider += l2_pool + CORRECTION;
        right_divider += r2_pool + CORRECTION;
    }

    /* optionally add final probabilities (input parameter), these are not pooled */
    if(FINAL_PROBABILITIES){
        left_divider += (double)l->num_final() + CORRECTION;
        right_divider += (double)r->num_final() + CORRECTION;
    }

    /* now we have the dividers and pools, we compute the likelihoods */
    for(type_num_map::iterator it = l->counts_begin(); it != l->counts_end(); it++){
        int type = it->first;
        num_map& nm = it->second;
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2) {
            int symbol = it2->first;
            left_count = it2->second;
            right_count = r->count(type, symbol);

            if (left_count >= SYMBOL_COUNT && right_count >= SYMBOL_COUNT)
                update_likelihood(left_count, right_count, left_divider, right_divider);
        }
    }
    
    /* count the pools */
    if(l1_pool >= SYMBOL_COUNT || r1_pool >= SYMBOL_COUNT){
        update_likelihood(l1_pool, r1_pool, left_divider, right_divider);
    }
    if(l2_pool >= SYMBOL_COUNT || r2_pool >= SYMBOL_COUNT){
        update_likelihood(l2_pool, r2_pool, left_divider, right_divider);
    }

    /* and always count final probabilities (are be pooled) */
    if(FINAL_PROBABILITIES){
        update_likelihood(l->num_final(), r->num_final(), left_divider, right_divider);
    }
};

bool likelihoodratio::compute_consistency(state_merger *merger, apta_node* left, apta_node* right){
    double test_statistic = 2.0 * (loglikelihood_orig - loglikelihood_merged);
    //double p_value = gsl_cdf_chisq_Q (test_statistic, 1.0 + (double)extra_parameters);
    double p_value = 1.0 - stats::pchisq(test_statistic, extra_parameters, false);

    //cerr << loglikelihood_orig << " " << loglikelihood_merged << " " << loglikelihood_orig - loglikelihood_merged << " " << extra_parameters << " " << p_value << endl;
    
    //cerr << "CHECK " << CHECK_PARAMETER << " " << (p_value < CHECK_PARAMETER) << endl;
    
    if (p_value < CHECK_PARAMETER) { return false; }

    if (inconsistency_found) return false;

    return true;
};

double likelihoodratio::compute_score(state_merger *merger, apta_node* left, apta_node* right){
    //if (inconsistency_found) return -1;

    double test_statistic = 2.0 * (loglikelihood_orig - loglikelihood_merged);
    double p_value1 = 1.0 - gsl_cdf_chisq_P (test_statistic, (double)extra_parameters);
    double p_value2 = gsl_cdf_chisq_Q (test_statistic, (double)extra_parameters);
    double p_value = 1.0 - stats::pchisq(test_statistic, extra_parameters, false);

    if(p_value2 > 0.05 && (p_value > p_value2 + 0.01 || p_value < p_value2 - 0.01)){
        cerr << endl;
        cerr << " values: " << test_statistic << " " << extra_parameters << " " << p_value << " " << p_value1 << " " << p_value2 << endl;
        cerr << "merge score: " << p_value << " " << loglikelihood_orig << " " << loglikelihood_merged << " " << extra_parameters << " " << p_value << endl;
    }

    return p_value2;
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
