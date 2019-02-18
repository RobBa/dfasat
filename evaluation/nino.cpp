#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <stdio.h>
#include <cmath>
#include <gsl/gsl_cdf.h>

#include "state_merger.h"
#include "evaluate.h"
#include "rti+.h"
#include "parameters.h"

REGISTER_DEF_DATATYPE(rtiplus_data);
REGISTER_DEF_TYPE(rtiplus);

vector< vector<double> > rtiplus::attribute_quantiles;

rtiplus_data::rtiplus_data() : likelihood_data::likelihood_data() {
    for(int i = 0; i < inputdata::num_attributes; ++i){
        quantile_counts.push_back(vector<int>(4,0));
    }
};

void rtiplus_data::read_from(tail* t){
    likelihood_data::read_from(inputdata::get_type(t),
                               inputdata::get_index(t),
                               inputdata::get_length(t),
                               inputdata::get_symbol(t),
                               inputdata::get_data(t));

    for(int i = 0; i < inputdata::num_attributes; ++i) {
        bool found = false;
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size(); ++j){
            if(inputdata::get_value(t, i) < rtiplus::attribute_quantiles[i][j]){
                quantile_counts[i][j] = quantile_counts[i][j] + 1;
                found = true;
                break;
            }
        }
        if(!found){
            quantile_counts[i][rtiplus::attribute_quantiles[i].size()] = quantile_counts[i][rtiplus::attribute_quantiles[i].size()] + 1;
        }
    }
};

void rtiplus_data::update(evaluation_data* right){
    likelihood_data::update(right);
    rtiplus_data* other = (rtiplus_data*)right;
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size()+1; ++j){
            quantile_counts[i][j] = quantile_counts[i][j] + other->quantile_counts[i][j];
        }
    }
};

void rtiplus_data::undo(evaluation_data* right){
    likelihood_data::undo(right);
    rtiplus_data* other = (rtiplus_data*)right;
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size()+1; ++j){
            quantile_counts[i][j] = quantile_counts[i][j] - other->quantile_counts[i][j];
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
    int symbol = inputdata::get_symbol(t);
    num_pos[symbol] = num_pos[symbol] - 1;
    other->num_pos[symbol] = other->num_pos[symbol] + 1;
    for(int i = 0; i < inputdata::num_attributes; ++i) {
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

bool rtiplus::consistent(state_merger *merger, apta_node* left, apta_node* right){
    return likelihoodratio::consistent(merger, left, right);
};

/* Likelihood Ratio (LR), computes an LR-test (used in RTI) and uses the p-value as score and consistency */

bool rtiplus::split_consistent(state_merger*, apta_node* left, apta_node* right){
    return true;
};

void rtiplus::add_merged_likelihood(double count, double divider){
    if(count != 0.0){
        //cerr << "adding m: " << count << " / " << divider << " " << (count + CORRECTION) * log((count + CORRECTION) / (divider)) << endl;
        //extra_parameters++;
        loglikelihood_merged += (count + CORRECTION) * log((count + CORRECTION) / (divider));
        //cerr << "ll_merged: " << loglikelihood_merged << " " << extra_parameters << endl;
    }
};

void rtiplus::del_merged_likelihood(double count, double divider){
    if(count != 0.0){
        //cerr << "deleting m: " << count << " / " << divider << " " << (count + CORRECTION) * log((count + CORRECTION) / (divider)) << endl;
        //extra_parameters++;
        loglikelihood_merged -= (count + CORRECTION) * log((count + CORRECTION) / (divider));
        //cerr << "ll_merged: " << loglikelihood_merged << " " << extra_parameters << endl;
    }
};

void rtiplus::add_split_likelihood(double count, double divider){
    if(count != 0.0){
        //cerr << "adding: " << count << " / " << divider << " " << (count + CORRECTION) * log((count + CORRECTION) / (divider)) << endl;
        //extra_parameters++;
        loglikelihood_orig += (count + CORRECTION) * log((count + CORRECTION) / (divider));
        //cerr << "ll_orig: " << loglikelihood_orig << " " << extra_parameters << endl;
    }
};

void rtiplus::del_split_likelihood(double count, double divider){
    if(count != 0.0){
        //cerr << "deleting: " << count << " / " << divider << " " << (count + CORRECTION) * log((count + CORRECTION) / (divider)) << endl;
        //extra_parameters--;
        loglikelihood_orig -= (count + CORRECTION) * log((count + CORRECTION) / (divider));
        //cerr << "ll_orig: " << loglikelihood_orig << " " << extra_parameters << endl;
    }
};

void rtiplus::add_parameters(rtiplus_data* l){
    double count = 0.0;
    double pool = 0.0;
    for(num_map::iterator it = l->num_pos.begin(); it != l->num_pos.end(); ++it){
        count = (*it).second;
        if(count >= SYMBOL_COUNT) extra_parameters++;
        else pool += count;
    }
    if(pool != 0) extra_parameters++;

    for(int i = 0; i < inputdata::num_attributes; ++i) {
        pool = 0;
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            count  = l->quantile_counts[i][j];
            if(count >= SYMBOL_COUNT) extra_parameters++;
            else pool += count;
        }
        if(pool != 0) extra_parameters++;
    }
};

void rtiplus::del_parameters(rtiplus_data* l){
    double count = 0.0;
    double pool = 0.0;
    for(num_map::iterator it = l->num_pos.begin(); it != l->num_pos.end(); ++it){
        count = (*it).second;
        if(count >= SYMBOL_COUNT) extra_parameters--;
        else pool += count;
    }
    if(pool != 0) extra_parameters--;
    
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        pool = 0;
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            count  = l->quantile_counts[i][j];
            if(count >= SYMBOL_COUNT) extra_parameters--;
            else pool += count;
        }
        if(pool != 0) extra_parameters--;
    }
};

void rtiplus::add_merged_likelihood(rtiplus_data* l){
    double divider = 0.0;
    double count = 0.0;
    double pool = 0.0;
    for(num_map::iterator it = l->num_pos.begin(); it != l->num_pos.end(); ++it){
        divider += (*it).second + CORRECTION;
    }
    for(num_map::iterator it = l->num_pos.begin(); it != l->num_pos.end(); ++it){
        count = (*it).second;
        if(count >= SYMBOL_COUNT) add_merged_likelihood(count, divider);
        else pool += count;
    }
    //if(pool >= SYMBOL_COUNT)
    add_merged_likelihood(pool, divider);

    for(int i = 0; i < inputdata::num_attributes; ++i) {
        pool = 0;
        divider = 0.0;
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            if(l->quantile_counts[i][j] != 0) divider += l->quantile_counts[i][j] + CORRECTION;
        }
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            count  = l->quantile_counts[i][j];
            if(count >= SYMBOL_COUNT) add_merged_likelihood(count, divider);
            else pool += count;
        }
        //if(pool >= SYMBOL_COUNT)
        add_merged_likelihood(pool, divider);
    }
};

void rtiplus::del_merged_likelihood(rtiplus_data* l){
    double divider = 0.0;
    double count = 0.0;
    double pool = 0.0;
    for(num_map::iterator it = l->num_pos.begin(); it != l->num_pos.end(); ++it){
        divider += (*it).second + CORRECTION;
    }
    for(num_map::iterator it = l->num_pos.begin(); it != l->num_pos.end(); ++it){
        count = (*it).second;
        if(count >= SYMBOL_COUNT) del_merged_likelihood(count, divider);
        else pool += count;
    }
    //if(pool >= SYMBOL_COUNT)
    del_merged_likelihood(pool, divider);

    for(int i = 0; i < inputdata::num_attributes; ++i) {
        pool = 0;
        divider = 0.0;
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            if(l->quantile_counts[i][j] != 0) divider += l->quantile_counts[i][j] + CORRECTION;
        }
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            count  = l->quantile_counts[i][j];
            if(count >= SYMBOL_COUNT) del_merged_likelihood(count, divider);
            else pool += count;
        }
        //if(pool >= SYMBOL_COUNT)
        del_merged_likelihood(pool, divider);
    }
};

void rtiplus::add_split_likelihood(rtiplus_data* l){
    double divider = 0.0;
    double count = 0.0;
    double pool = 0.0;
    for(num_map::iterator it = l->num_pos.begin(); it != l->num_pos.end(); ++it){
        divider += (*it).second + CORRECTION;
    }
    for(num_map::iterator it = l->num_pos.begin(); it != l->num_pos.end(); ++it){
        count = (*it).second;
        if(count >= SYMBOL_COUNT) add_split_likelihood(count, divider);
        else pool += count;
    }
    //if(pool >= SYMBOL_COUNT)
    add_split_likelihood(pool, divider);

    for(int i = 0; i < inputdata::num_attributes; ++i) {
        pool = 0;
        divider = 0.0;
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            if(l->quantile_counts[i][j] != 0) divider += l->quantile_counts[i][j] + CORRECTION;
        }
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            count  = l->quantile_counts[i][j];
            if(count >= SYMBOL_COUNT) add_split_likelihood(count, divider);
            else pool += count;
        }
        //if(pool >= SYMBOL_COUNT)
        add_split_likelihood(pool, divider);
    }
};

void rtiplus::del_split_likelihood(rtiplus_data* l){
    double divider = 0.0;
    double count = 0.0;
    double pool = 0.0;
    for(num_map::iterator it = l->num_pos.begin(); it != l->num_pos.end(); ++it){
        divider += (*it).second + CORRECTION;
    }
    for(num_map::iterator it = l->num_pos.begin(); it != l->num_pos.end(); ++it){
        count = (*it).second;
        if(count >= SYMBOL_COUNT) del_split_likelihood(count, divider);
        else pool += count;
    }
    //if(pool >= SYMBOL_COUNT)
    del_split_likelihood(pool, divider);
    
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        pool = 0;
        divider = 0.0;
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            if(l->quantile_counts[i][j] != 0) divider += l->quantile_counts[i][j] + CORRECTION;
        }
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            count  = l->quantile_counts[i][j];
            if(count >= SYMBOL_COUNT) del_split_likelihood(count, divider);
            else pool += count;
        }
        //if(pool >= SYMBOL_COUNT)
        del_split_likelihood(pool, divider);
    }
};

void rtiplus::update_score(state_merger *merger, apta_node* left, apta_node* right){
    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;
    
    /* cerr << "adding left: " << left->size << endl;
    for(int a = 0; a < inputdata::alphabet.size(); ++a){
        cerr << l->pos(a) << " ";
    }
    cerr << endl;
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            cerr << l->quantile_counts[i][j] << " ";
        }
        cerr << endl;
    }
    cerr << "adding right: " << right->size << endl;
    for(int a = 0; a < inputdata::alphabet.size(); ++a){
        cerr << r->pos(a) << " ";
    }
    cerr << endl;
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            cerr << r->quantile_counts[i][j] << " ";
        }
        cerr << endl;
    }*/
    add_split_likelihood(l);
    add_split_likelihood(r);
    
    return;
};

void rtiplus::update_score_after(state_merger *merger, apta_node* left, apta_node* right){
    rtiplus_data* l = (rtiplus_data*) left->data;

    /*cerr << "adding merged: " << left->size << endl;
    for(int a = 0; a < inputdata::alphabet.size(); ++a){
        cerr << l->pos(a) << " ";
    }
    cerr << endl;
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            cerr << l->quantile_counts[i][j] << " ";
        }
        cerr << endl;
    }*/
    add_merged_likelihood(l);
    add_parameters(l);

    return;
};

void rtiplus::split_update_score_before(state_merger*, apta_node* left, apta_node* right, tail* t){
    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;
    
    if(right->size == 0){
        add_parameters(l);
    }

    /*cerr << "deleting left:" << endl;
    for(int a = 0; a < inputdata::alphabet.size(); ++a){
        cerr << l->pos(a) << " ";
    }
    cerr << endl;
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            cerr << l->quantile_counts[i][j] << " ";
        }
        cerr << endl;
    }
    cerr << "deleting right:" << endl;
    for(int a = 0; a < inputdata::alphabet.size(); ++a){
        cerr << r->pos(a) << " ";
    }
    cerr << endl;
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            cerr << r->quantile_counts[i][j] << " ";
        }
        cerr << endl;
    }*/

    //cerr << "loglikelihood " << loglikelihood_orig << endl;
    del_split_likelihood(l);
    //cerr << "loglikelihood " << loglikelihood_orig << endl;
    del_split_likelihood(r);
    //cerr << "loglikelihood " << loglikelihood_orig << endl;
};

void rtiplus::split_update_score_after(state_merger*, apta_node* left, apta_node* right, tail* t){
    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;
    
    if(left->size == 0){
        del_parameters(r);
    }
    
    /*cerr << "adding left:" << endl;
    for(int a = 0; a < inputdata::alphabet.size(); ++a){
        cerr << l->pos(a) << " ";
    }
    cerr << endl;
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            cerr << l->quantile_counts[i][j] << " ";
        }
        cerr << endl;
    }
    cerr << "adding right:" << endl;
    for(int a = 0; a < inputdata::alphabet.size(); ++a){
        cerr << r->pos(a) << " ";
    }
    cerr << endl;
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            cerr << r->quantile_counts[i][j] << " ";
        }
        cerr << endl;
    }*/
    
    //cerr << "loglikelihood " << loglikelihood_orig << endl;
    add_split_likelihood(l);
    //cerr << "loglikelihood " << loglikelihood_orig << endl;
    add_split_likelihood(r);
    //cerr << "loglikelihood " << loglikelihood_orig << endl;
};

bool rtiplus::split_compute_consistency(state_merger *, apta_node* left, apta_node* right){
    if (inconsistency_found) return false;
    
    double test_statistic = 2.0 * (loglikelihood_orig - loglikelihood_merged);
    double p_value = gsl_cdf_chisq_Q (test_statistic, (double)extra_parameters);
    
    if (p_value > CHECK_PARAMETER) { inconsistency_found = true; return false; }
    
    return true;
};

double  rtiplus::split_compute_score(state_merger *, apta_node* left, apta_node* right){
    if (inconsistency_found) return -1;
    
    double test_statistic = 2.0 * (loglikelihood_orig - loglikelihood_merged);
    double p_value = gsl_cdf_chisq_Q (test_statistic, (double)extra_parameters);
    
    cerr << "split score: " << loglikelihood_orig << " " << loglikelihood_merged << " " << extra_parameters << " " << p_value << endl;

    return 1.0 - p_value;
};


void rtiplus::initialize_globals(){
    CORRECTION = 1.0;
    //rtiplus::attribute_quantiles = vector< vector<double> >(inputdata::num_attributes, vector<double>(4,0.0));
    for(int a = 0; a < inputdata::num_attributes; ++a){
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
        int Q4 = values.size();
        
        int V1 = 0;
        int V2 = 0;
        int V3 = 0;
        int V4 = 0;
        
        int count = 0;
        for(multiset<double>::iterator it = values.begin(); it != values.end(); ++it){
            if(count == Q1) V1 = *it;
            if(count == Q2) V2 = *it;
            if(count == Q3) V3 = *it;
            if(count == Q4) V4 = *it;
            count = count + 1;
        }
        
        rtiplus::attribute_quantiles[a][0] = Q1;
        rtiplus::attribute_quantiles[a][1] = Q2;
        rtiplus::attribute_quantiles[a][2] = Q3;
    }
    /*cerr << "quantiles:" << endl;
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size(); ++j){
            cerr << rtiplus::attribute_quantiles[i][j] << " ";
        }
        cerr << endl;
    }*/
};

void rtiplus::reset_split(state_merger *merger, apta_node* node){
    inconsistency_found = false;
    loglikelihood_orig = 0;
    loglikelihood_merged = 0;
    extra_parameters = 0;
    
    return;
    
    rtiplus_data* l = (rtiplus_data*) node->data;

    double left_divider = 1.0;
    double left_count = 0.0;
    double pool = 0;
    
    for(int a = 0; a < inputdata::alphabet.size(); ++a) left_divider += CORRECTION;
    left_divider += (double)l->accepting_paths;
    
    for(num_map::iterator it = l->num_pos.begin(); it != l->num_pos.end(); ++it){
        left_count = (*it).second;
        
        if(left_count >= SYMBOL_COUNT) add_merged_likelihood(left_count, left_divider);
        else pool += left_count;
    }
    if(pool >= SYMBOL_COUNT) add_merged_likelihood(pool, left_divider);
    
    for(int i = 0; i < inputdata::num_attributes; ++i){
        left_divider = 0;
        left_count = 0;
        pool = 0;
        
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size(); ++j) left_divider += CORRECTION;
        left_divider += (double)l->accepting_paths;// + (double)l->num_accepting;
        
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size(); ++j){
            left_count = l->quantile_counts[i][j];
            
            if(left_count >= SYMBOL_COUNT) add_merged_likelihood(left_count, left_divider);
            else pool += left_count;
        }
        if(pool >= SYMBOL_COUNT) add_merged_likelihood(pool, left_divider);
    }
    loglikelihood_orig = loglikelihood_merged;
};


/*void likelihoodratio::print_dot(iostream& output, state_merger* merger){
    alergia::print_dot(output, merger);
};*/
