
#define STATS_GO_INLINE

#include <math.h>
#include <vector>
#include <map>

#include "utility/stats.hpp"

#include "state_merger.h"
#include "evaluate.h"
#include "rti+.h"
#include "parameters.h"

REGISTER_DEF_DATATYPE(rtiplus_data);
REGISTER_DEF_TYPE(rtiplus);

vector< vector<double> > rtiplus::attribute_quantiles;

void check_counts(apta_node* node){
    rtiplus_data* dat = (rtiplus_data*)node->data;
    
    int sym_count = 0;
    int fin_count = dat->num_final();
    for(int i = 0; i < inputdata::types.size(); ++i) {
        for(int j = 0; j < inputdata::alphabet.size(); ++j){
            sym_count += dat->count(i,j);
        }
    }
    if(sym_count + fin_count != node->size) cerr << sym_count << " " << fin_count << " " << node->size << endl;
}

rtiplus_data::rtiplus_data() : likelihood_data::likelihood_data() {
    for(int i = 0; i < inputdata::num_attributes; ++i){
        if(inputdata::is_distributionable(i)) quantile_counts.push_back(vector<int>(4,0));
    }
    loglikelihood = 0.0;
};

void rtiplus_data::initialize() {
    likelihood_data::initialize();
    int modifier = 0;
    for(int i = 0; i < inputdata::num_attributes; ++i){
        if(!inputdata::is_distributionable(i)){
            ++modifier;
            continue;
        }
        quantile_counts[i - modifier].assign(4, 0);
    }
    loglikelihood = 0.0;
};

void rtiplus_data::read_from(tail* t){
    alergia_data::read_from(inputdata::get_type(t),
                               inputdata::get_index(t),
                               inputdata::get_length(t),
                               inputdata::get_symbol(t),
                               inputdata::get_data(t));
    int modifier = 0;

    for(int i = 0; i < inputdata::num_attributes; ++i) {
        if(!inputdata::is_distributionable(i)){
            ++modifier;
            continue;
        }
        int attr = i-modifier;
        bool found = false;
        for(int j = 0; j < rtiplus::attribute_quantiles[attr].size(); ++j){
            if(inputdata::get_value(t, i) < rtiplus::attribute_quantiles[attr][j]){
                quantile_counts[attr][j] = quantile_counts[attr][j] + 1;
                found = true;
                break;
            }
        }
        if(!found){
            quantile_counts[attr][rtiplus::attribute_quantiles[attr].size()] = quantile_counts[attr][rtiplus::attribute_quantiles[attr].size()] + 1;
        }
    }
};

void rtiplus_data::print_state_label(iostream& output, apta* aptacontext){
    for(int i = 0; i < inputdata::types.size(); ++i) {
        output << "fin(" << i << "):";
        output << num_final() << endl;
        output << "" << endl;
    }
    /*
    for(int i = 0; i < inputdata::types.size(); ++i) {
        output << "symb(" << i << "):[";
        for(int j = 0; j < inputdata::alphabet.size(); ++j){
            output << count(i,j) << ",";
        }
        output << "]" << endl;
    }
    */
     for(int i = 0; i < quantile_counts.size(); ++i) {
        output << "attr(" << i << "):[";
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size()+1; ++j){
            output << quantile_counts[i][j] << ",";
        }
        output << "]" << endl;
    }
};

void rtiplus_data::update(evaluation_data* right){
    likelihood_data::update(right);
    rtiplus_data* other = (rtiplus_data*)right;
    for(int i = 0; i < quantile_counts.size(); ++i) {
        if(!inputdata::is_distributionable(i)) continue;
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size()+1; ++j){
            quantile_counts[i][j] += other->quantile_counts[i][j];
        }
    }
};

void rtiplus_data::undo(evaluation_data* right){
    likelihood_data::undo(right);
    rtiplus_data* other = (rtiplus_data*)right;
    for(int i = 0; i < quantile_counts.size(); ++i) {
        if(!inputdata::is_distributionable(i)) continue;
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size()+1; ++j){
            quantile_counts[i][j] -= other->quantile_counts[i][j];
        }
    }
};

void rtiplus_data::split_update(evaluation_data* right){
    undo(right);
};

void rtiplus_data::split_undo(evaluation_data* right){
    update(right);
};

void rtiplus_data::split_update_single(evaluation_data* right, tail* t){
    rtiplus_data* other = (rtiplus_data*)right;
    int type  = inputdata::get_type(t);
    int symbol = inputdata::get_symbol(t);
    
    total_paths -= 1;
    other->total_paths += 1;
    trans_counts[type][symbol] = trans_counts[type][symbol] - 1;
    other->trans_counts[type][symbol] = other->trans_counts[type][symbol] + 1;
    for(int i = 0; i < quantile_counts.size(); ++i) {
        bool found = false;
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size(); ++j){
            if(inputdata::get_value(t, i) < rtiplus::attribute_quantiles[i][j]){
                quantile_counts[i][j] = quantile_counts[i][j] - 1;
                other->quantile_counts[i][j] =  other->quantile_counts[i][j] + 1;
                found = true;
                break;
            }
        }
        if(!found){
            quantile_counts[i][rtiplus::attribute_quantiles[i].size()] = quantile_counts[i][rtiplus::attribute_quantiles[i].size()] - 1;
            other->quantile_counts[i][rtiplus::attribute_quantiles[i].size()] = other->quantile_counts[i][rtiplus::attribute_quantiles[i].size()] + 1;
        }
    }
};

void rtiplus_data::del_tail(tail* t){
    int type  = inputdata::get_type(t);
    int symbol = inputdata::get_symbol(t);
    
    if(t->get_index() == -1){
        final_counts[type]--;
        total_final--;
    } else {
        total_paths -= 1;
        trans_counts[type][symbol] = trans_counts[type][symbol] - 1;
        for(int i = 0; i < quantile_counts.size(); ++i) {
            bool found = false;
            for(int j = 0; j < rtiplus::attribute_quantiles[i].size(); ++j){
                if(inputdata::get_value(t, i) < rtiplus::attribute_quantiles[i][j]){
                    quantile_counts[i][j] = quantile_counts[i][j] - 1;
                    found = true;
                    break;
                }
            }
            if(!found){
                quantile_counts[i][rtiplus::attribute_quantiles[i].size()] = quantile_counts[i][rtiplus::attribute_quantiles[i].size()] - 1;
            }
        }
    }
};

void rtiplus_data::set_loglikelihood(){
    loglikelihood = 0.0;

    double divider = 0.0;
    double pool = 0.0;
    for(type_num_map::iterator it = counts_begin(); it != counts_end(); it++){
        num_map& nm = it->second;
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            double count = it2->second;
            if(count >= SYMBOL_COUNT) divider += count + CORRECTION;
            else pool += count;
        }
    }

    if(FINAL_PROBABILITIES){
        double count = num_final();
        if(count >= SYMBOL_COUNT) divider += count + CORRECTION;
        else pool += count;
    }

    if(pool > 0.0) divider += pool + CORRECTION;

    for(type_num_map::iterator it = counts_begin(); it != counts_end(); it++){
        num_map& nm = it->second;
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            double count = it2->second;
            if(count < SYMBOL_COUNT) ;
            else if(count > 0) loglikelihood += (count + CORRECTION) * log((count + CORRECTION) / divider);
        }
    }

    if(FINAL_PROBABILITIES){
        double count = num_final();
        if(count < SYMBOL_COUNT) ;
        else if(count > 0) loglikelihood += (count + CORRECTION) * log((count + CORRECTION) / divider);
    }

    if(pool > 0) loglikelihood += (pool + CORRECTION) * log((pool + CORRECTION) / divider);

    divider = 0.0;
    pool = 0.0;
    for(int i = 0; i < quantile_counts.size(); ++i) {
        for (int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j) {
            double count = quantile_counts[i][j];
            if (count >= SYMBOL_COUNT) divider += count + CORRECTION;
            else pool += count;
        }
    }
    if(pool > 0) divider += pool + CORRECTION;

    for(int i = 0; i < quantile_counts.size(); ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            double count = quantile_counts[i][j];
            if (count >= SYMBOL_COUNT) loglikelihood += (count + CORRECTION) * log((count + CORRECTION) / divider);
        }
    }
    if(pool > 0) loglikelihood += (pool + CORRECTION) * log((pool + CORRECTION) / divider);
};

int rtiplus_data::num_parameters(){
    int result = 1;
    for(type_num_map::iterator it = counts_begin(); it != counts_end(); it++){
        num_map& nm = it->second;
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            double count = it2->second;
            if(count != 0) result++;
        }
    }

    if(FINAL_PROBABILITIES) {
        result += 1;
    }

    for(int i = 0; i < quantile_counts.size(); ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            double count = quantile_counts[i][j];
            if(count != 0) result++;
        }
    }
    return result;
};

/* Likelihood Ratio (LR), computes an LR-test (used in RTI) and uses the p-value as score and consistency */
void rtiplus::update_score(state_merger *merger, apta_node* left, apta_node* right){
    double temp_loglikelihood_orig = loglikelihood_orig;
    double temp_loglikelihood_merged = loglikelihood_merged;
    int temp_extra_parameters = extra_parameters;

    likelihoodratio::update_score(merger, left, right);

    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;


    /* we ignore low frequency states, decided by input parameter STATE_COUNT */
    if(FINAL_PROBABILITIES) if(r->num_paths() + r->num_final() < STATE_COUNT || l->num_paths() + l->num_final() < STATE_COUNT) return;
    else if(r->num_paths() < STATE_COUNT || l->num_paths() < STATE_COUNT) return;

    /* we treat type distributions as independent */
    for(int i = 0; i < l->quantile_counts.size(); ++i) {
        /* computing the dividers (denominator) */
        double left_divider = 0.0;
        double right_divider = 0.0;
        double left_count = 0.0;
        double right_count  = 0.0;

        double l1_pool = 0.0;
        double r1_pool = 0.0;
        double l2_pool = 0.0;
        double r2_pool = 0.0;

        for (int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j) {
            left_count = l->quantile_counts[i][j];
            right_count = r->quantile_counts[i][j];

            update_divider(left_count, right_count, left_divider, right_divider);
            update_left_pool(left_count, right_count, l1_pool, r1_pool);
            update_right_pool(left_count, right_count, l2_pool, r2_pool);
        }

        update_divider_pool(l1_pool, r1_pool, left_divider, right_divider);
        update_divider_pool(l2_pool, r2_pool, left_divider, right_divider);

        if(left_divider < STATE_COUNT || right_divider < STATE_COUNT) continue;

        for (int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j) {
            left_count = l->quantile_counts[i][j];
            right_count = r->quantile_counts[i][j];

            update_likelihood(left_count, right_count, left_divider, right_divider);
        }
        update_likelihood_pool(l1_pool, r1_pool, left_divider, right_divider);
        update_likelihood_pool(l2_pool, r2_pool, left_divider, right_divider);
    }

    r->undo_loglikelihood_orig = loglikelihood_orig - temp_loglikelihood_orig;
    r->undo_loglikelihood_merged = loglikelihood_merged - temp_loglikelihood_merged;
    r->undo_extra_parameters = extra_parameters - temp_extra_parameters;
};

/*
void rtiplus::update_score_after(state_merger *merger, apta_node* left, apta_node* right){
    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;

    l->set_loglikelihood();
    
    if(l->num_paths() < STATE_COUNT) return;

    loglikelihood_merged += l->loglikelihood;
    extra_parameters -= l->num_parameters();
};

void rtiplus::split_update_score_before(state_merger* merger, apta_node* left, apta_node* right, tail* t){
    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;

    l->set_loglikelihood();
    r->set_loglikelihood();

    if(l->num_paths() + r->num_paths() < STATE_COUNT) return;

    if(right->size == 0){
        extra_parameters += l->num_parameters();    
        loglikelihood_merged += l->loglikelihood;
    } else {
        loglikelihood_orig -= l->loglikelihood;
        loglikelihood_orig -= r->loglikelihood;
    }
};
*/

void rtiplus::split_update_score_before(state_merger* merger, apta_node* left, apta_node* right, tail* t) {
    rtiplus_data *l = (rtiplus_data *) left->data;
    rtiplus_data *r = (rtiplus_data *) right->data;

    loglikelihood_orig -= r->undo_loglikelihood_orig;
    loglikelihood_merged -= r->undo_loglikelihood_merged;
    extra_parameters -= r->undo_extra_parameters;

    r->undo_loglikelihood_orig = 0.0;
    r->undo_loglikelihood_merged = 0.0;
    r->undo_extra_parameters = 0;
};

void rtiplus::split_update_score_after(state_merger* merger, apta_node* left, apta_node* right, tail* t) {
    rtiplus_data *l = (rtiplus_data *) left->data;
    rtiplus_data *r = (rtiplus_data *) right->data;

    update_score(merger, left, right);
};
/*

    l->set_loglikelihood();
    r->set_loglikelihood();
    
    if(l->num_paths() + r->num_paths() < STATE_COUNT) return;

    if(left->size == 0){
        extra_parameters -= r->num_parameters();
        loglikelihood_merged -= r->loglikelihood;
    }

    loglikelihood_orig += l->loglikelihood;
    loglikelihood_orig += r->loglikelihood;
};
 */

bool rtiplus::split_compute_consistency(state_merger *, apta_node* left, apta_node* right){
    double test_statistic = 2.0 * (loglikelihood_orig - loglikelihood_merged);
    double p_value = 1.0 - stats::pchisq(test_statistic, extra_parameters, false);

    if(left->size <= STATE_COUNT || right->size <= STATE_COUNT) return false;
    if(USE_SINKS && (left->size <= SINK_COUNT || right->size <= SINK_COUNT)) return false;
    if (p_value > CHECK_PARAMETER) return false;
    
    if (inconsistency_found) return false;
    
    return true;
};

double rtiplus::split_compute_score(state_merger *, apta_node* left, apta_node* right){
    double test_statistic = 2.0 * (loglikelihood_orig - loglikelihood_merged);
    double p_value = 1.0 - stats::pchisq(test_statistic, extra_parameters, false);

    return 1.0 + CHECK_PARAMETER - p_value;
};

void rtiplus::initialize(state_merger* m){
    for(merged_APTA_iterator it = merged_APTA_iterator(m->aut->root); *it != 0; ++it){
        apta_node* node = *it;
        rtiplus_data* n = (rtiplus_data*)node->data;
        n->set_loglikelihood();
    }
};

void rtiplus::initialize_globals(){
    CORRECTION = 0.0;
    for(int a = 0; a < inputdata::num_attributes; ++a){
        if(!inputdata::is_distributionable(a)) continue;
        rtiplus::attribute_quantiles.push_back(vector<double>(3,0.0));
        multiset<double> values;
        for(int i = 0; i < inputdata::get_size(); ++i){
            for(int j = 0; j < inputdata::get_length(i); ++j){
                values.insert(inputdata::get_value(i, j, a));
            }
        }
        
        int Q1 = (int)round((double)values.size() / 4.0);
        int Q2 = Q1 + Q1;
        int Q3 = Q2 + Q1;

        int V1 = 0;
        int V2 = 0;
        int V3 = 0;

        int count = 0;
        for(multiset<double>::iterator it = values.begin(); it != values.end(); ++it){
            if(count == Q1) V1 = *it;
            if(count == Q2) V2 = *it;
            if(count == Q3) V3 = *it;
            count = count + 1;
        }
        
        rtiplus::attribute_quantiles[a][0] = V1;
        rtiplus::attribute_quantiles[a][1] = V2;
        rtiplus::attribute_quantiles[a][2] = V3;
    }
}

void rtiplus::reset_split(state_merger *merger, apta_node* node){
    inconsistency_found = false;
    loglikelihood_orig = 0;
    loglikelihood_merged = 0;
    extra_parameters = 0;
       
    return;
};  
