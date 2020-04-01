#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <stdio.h>
#include <cmath>
#include <gsl/gsl_cdf.h>
#include <cassert>

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
    assert(sym_count + fin_count == node->size);
}

rtiplus_data::rtiplus_data() : likelihood_data::likelihood_data() {
    for(int i = 0; i < inputdata::num_attributes; ++i){
        quantile_counts.push_back(vector<int>(4,0));
    }
    loglikelihood = 0.0;
};

void rtiplus_data::read_from(tail* t){
    alergia_data::read_from(inputdata::get_type(t),
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

void rtiplus_data::print_state_label(iostream& output, apta* aptacontext){
    for(int i = 0; i < inputdata::types.size(); ++i) {
        output << "fin(" << i << "):";
        output << num_final() << endl;
        output << "" << endl;
    }    for(int i = 0; i < inputdata::types.size(); ++i) {
        output << "symb(" << i << "):[";
        for(int j = 0; j < inputdata::alphabet.size(); ++j){
            output << count(i,j) << ",";
        }
        output << "]" << endl;
    }
    for(int i = 0; i < inputdata::num_attributes; ++i) {
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
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size()+1; ++j){
            quantile_counts[i][j] += other->quantile_counts[i][j];
        }
    }
};

void rtiplus_data::undo(evaluation_data* right){
    likelihood_data::undo(right);
    rtiplus_data* other = (rtiplus_data*)right;
    for(int i = 0; i < inputdata::num_attributes; ++i) {
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

void rtiplus_data::del_tail(tail* t){
    int type  = inputdata::get_type(t);
    int symbol = inputdata::get_symbol(t);
    
    if(t->index == -1){
        final_counts[type]--;
        total_final--;
    } else {
        total_paths -= 1;
        trans_counts[type][symbol] = trans_counts[type][symbol] - 1;
        for(int i = 0; i < inputdata::num_attributes; ++i) {
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
    double divider = (double)num_paths();
    
    if(divider == 0) return;
    
    for(type_num_map::iterator it = counts_begin(); it != counts_end(); it++){
        num_map& nm = (*it).second;
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            double count = (*it2).second;
            if(count != 0) loglikelihood += count * log(count / divider);
        }
    }

    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            double count = quantile_counts[i][j];
            if(count != 0) loglikelihood += count * log(count / divider);
        }
    }
};

int rtiplus_data::num_parameters(){
    int result = 0;
    for(type_num_map::iterator it = counts_begin(); it != counts_end(); it++){
        num_map& nm = (*it).second;
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            double count = (*it2).second;
            if(count != 0) result++;
        }
    }

    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            double count = quantile_counts[i][j];
            if(count != 0) result++;
        }
    }
    return result;
};

/* Likelihood Ratio (LR), computes an LR-test (used in RTI) and uses the p-value as score and consistency */
void rtiplus::update_score(state_merger *merger, apta_node* left, apta_node* right){
    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;

    l->set_loglikelihood();
    r->set_loglikelihood();

    check_counts(left);
    check_counts(right);
    
    if(l->num_paths() + r->num_paths() < STATE_COUNT) return;

    /*cerr << "adding left: " << left->size << endl;
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

    loglikelihood_orig += l->loglikelihood;
    loglikelihood_orig += r->loglikelihood;
    extra_parameters += l->num_parameters();
    extra_parameters += r->num_parameters();    
};

void rtiplus::update_score_after(state_merger *merger, apta_node* left, apta_node* right){
    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;

    check_counts(left);
    check_counts(right);

    l->set_loglikelihood();
    
    if(l->num_paths() < STATE_COUNT) return;

    loglikelihood_merged += l->loglikelihood;
    extra_parameters -= l->num_parameters();
    //cerr << "merge score " << compute_score(merger, left, right) << " " << loglikelihood_orig << " " << loglikelihood_merged << " " << extra_parameters << endl;
};

void rtiplus::split_update_score_before(state_merger* merger, apta_node* left, apta_node* right, tail* t){
    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;

    check_counts(left);
    check_counts(right);

    l->set_loglikelihood();
    r->set_loglikelihood();

    if(l->num_paths() + r->num_paths() < STATE_COUNT) return;

    if(right->size == 0){
        extra_parameters += l->num_parameters();    
        //extra_parameters += 5;    
        loglikelihood_merged += l->loglikelihood;
    } else {
        //cerr << "before1 " << loglikelihood_orig << endl;
        loglikelihood_orig -= l->loglikelihood;
        loglikelihood_orig -= r->loglikelihood;
        //cerr << "before2 " << loglikelihood_orig << endl;
        //extra_parameters -= l->num_parameters();
        //extra_parameters -= r->num_parameters();
    }
};

void rtiplus::split_update_score_after(state_merger* merger, apta_node* left, apta_node* right, tail* t){
    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;

    check_counts(left);
    check_counts(right);

    l->set_loglikelihood();
    r->set_loglikelihood();
    
    if(l->num_paths() + r->num_paths() < STATE_COUNT) return;

    if(left->size == 0){
        //extra_parameters -= 5;    
        extra_parameters -= r->num_parameters();    
        loglikelihood_merged -= r->loglikelihood;
    }

    /*
     * cerr << "splitting left: " << left->size << endl;
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
    cerr << "splitting right: " << right->size << endl;
    for(int a = 0; a < inputdata::alphabet.size(); ++a){
        cerr << r->pos(a) << " ";
    }
    cerr << endl;
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            cerr << r->quantile_counts[i][j] << " ";
        }
        cerr << endl;
    }
     * */
    
    //cerr << "after1 " << loglikelihood_orig << endl;
    loglikelihood_orig += l->loglikelihood;
    loglikelihood_orig += r->loglikelihood;
    //cerr << "after2 " << loglikelihood_orig << endl;
    //extra_parameters += l->num_parameters();
    //extra_parameters += r->num_parameters();
    //cerr << "split score " << split_compute_score(merger, left, right) << " " << loglikelihood_orig << " " << loglikelihood_merged << " " << extra_parameters << endl;
};

/*    likelihoodratio::split_update_score_before(merger, left, right, t);
    
    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;

    if(r->num_paths() < STATE_COUNT || l->num_paths() < STATE_COUNT || r->num_paths() < 1 || l->num_paths() < 1) return;
    
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        double left_divider = (double)l->num_paths();
        double right_divider = (double)r->num_paths();
        double left_count = 0.0;
        double right_count  = 0.0;
        
        int l1_pool = 0;
        int r1_pool = 0;
        int l2_pool = 0;
        int r2_pool = 0;
        int matching_right = 0;

        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            left_count = l->quantile_counts[i][j];
            right_count = r->quantile_counts[i][j];
            
            if(left_count >= SYMBOL_COUNT || right_count >= SYMBOL_COUNT){
                left_divider += CORRECTION;
                right_divider += CORRECTION;
            }
        }

        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            left_count = l->quantile_counts[i][j];
            right_count = r->quantile_counts[i][j];
    
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
        r2_pool += r->num_paths() - matching_right;
    
        left_count = l1_pool;
        right_count = r1_pool;
    
        if(left_count >= SYMBOL_COUNT || right_count >= SYMBOL_COUNT)
            update_likelihood(left_count, right_count, left_divider, right_divider);
    
        left_count = l2_pool;
        right_count = r2_pool;
    
        if(left_count >= SYMBOL_COUNT || right_count >= SYMBOL_COUNT)
            update_likelihood(left_count, right_count, left_divider, right_divider);   
    }
}*/

/*    //undo_update_score(merger, left, right);
    
    //return;
    
    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;
    
    inconsistency_found = false;
    
    if(right->size == 0){
        add_parameters(l);
    }

    cerr << "deleting left:" << endl;
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
    }

    //cerr << "loglikelihood " << loglikelihood_orig << endl;
    del_split_likelihood(l);
    //cerr << "loglikelihood " << loglikelihood_orig << endl;
    del_split_likelihood(r);
    //cerr << "loglikelihood " << loglikelihood_orig << endl;
};*/

/*void rtiplus::split_update_score_after(state_merger* merger, apta_node* left, apta_node* right, tail* t){
    loglikelihood_orig += left->loglikelihood;
    loglikelihood_orig += right->loglikelihood;
    return;
    //update_score(merger, left, right);
    
    //return;
    
    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;
    
    if(left->size == 0){
        del_parameters(r);
    }
    
    cerr << "adding left:" << endl;
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
    }
    
    //cerr << "loglikelihood " << loglikelihood_orig << endl;
    add_split_likelihood(l);
    //cerr << "loglikelihood " << loglikelihood_orig << endl;
    add_split_likelihood(r);
    //cerr << "loglikelihood " << loglikelihood_orig << endl;
    //cerr << "intermediate ll: " << loglikelihood_orig << " " << loglikelihood_merged << " " << extra_parameters << endl;
};*/

bool rtiplus::split_compute_consistency(state_merger *, apta_node* left, apta_node* right){
    double test_statistic = 2.0 * (loglikelihood_orig - loglikelihood_merged);
    double p_value = gsl_cdf_chisq_Q (test_statistic, (double)extra_parameters);
    
    //cerr << "split score: " << p_value << " " << loglikelihood_orig << " " << loglikelihood_merged << " " << extra_parameters << " " << p_value << endl;

    if(left->size < STATE_COUNT || right->size < STATE_COUNT) return false;

    if (p_value > CHECK_PARAMETER) { return false; }
    
    //cerr << "true " << inconsistency_found << endl;
    if (inconsistency_found) return false;
    
    return true;
};

double  rtiplus::split_compute_score(state_merger *, apta_node* left, apta_node* right){
    double test_statistic = 2.0 * (loglikelihood_orig - loglikelihood_merged);
    double p_value = gsl_cdf_chisq_Q (test_statistic, (double)extra_parameters);
    
    //cerr << "split score: " << p_value << " " << loglikelihood_orig << " " << loglikelihood_merged << " " << extra_parameters << " " << p_value << endl;

    //if (inconsistency_found) return -1;
    
    if(left->size < STATE_COUNT || right->size < STATE_COUNT) return 0.0;
    
    //cerr << "score " << 1.0 + CHECK_PARAMETER - p_value << endl;

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
        
        rtiplus::attribute_quantiles[a][0] = V1;
        rtiplus::attribute_quantiles[a][1] = V2;
        rtiplus::attribute_quantiles[a][2] = V3;
    }
    cerr << "quantiles:" << endl;
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size(); ++j){
            cerr << rtiplus::attribute_quantiles[i][j] << " ";
        }
        cerr << endl;
    }
}

void rtiplus::reset_split(state_merger *merger, apta_node* node){
    inconsistency_found = false;
    loglikelihood_orig = 0;
    loglikelihood_merged = 0;
    extra_parameters = 0;
       
    return;
};  
/*    rtiplus_data* l = (rtiplus_data*) node->data;

    double left_divider = 1.0;
    double left_count = 0.0;
    double pool = 0;
    
    for(int a = 0; a < inputdata::alphabet.size(); ++a) left_divider += CORRECTION;
    left_divider += (double)l->pos_paths();
    
    for(num_map::iterator it = l->pos_begin(); it != l->pos_end(); ++it){
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
        left_divider += (double)l->pos_paths();// + (double)l->pos_final();
        
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size(); ++j){
            left_count = l->quantile_counts[i][j];
            
            if(left_count >= SYMBOL_COUNT) add_merged_likelihood(left_count, left_divider);
            else pool += left_count;
        }
        if(pool >= SYMBOL_COUNT) add_merged_likelihood(pool, left_divider);
    }
    loglikelihood_orig = loglikelihood_merged;
};*/


/*    likelihoodratio::update_score(merger, left, right);
    
    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;

    if(r->num_paths() < STATE_COUNT || l->num_paths() < STATE_COUNT || r->num_paths() < 1 || l->num_paths() < 1) return;
    
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        double left_divider = (double)l->num_paths();
        double right_divider = (double)r->num_paths();
        double left_count = 0.0;
        double right_count  = 0.0;
        
        int l1_pool = 0;
        int r1_pool = 0;
        int l2_pool = 0;
        int r2_pool = 0;
        int matching_right = 0;

        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            left_count = l->quantile_counts[i][j];
            right_count = r->quantile_counts[i][j];
            
            if(left_count >= SYMBOL_COUNT || right_count >= SYMBOL_COUNT){
                left_divider += CORRECTION;
                right_divider += CORRECTION;
            }
        }

        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            left_count = l->quantile_counts[i][j];
            right_count = r->quantile_counts[i][j];
    
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
        r2_pool += r->num_paths() - matching_right;
    
        left_count = l1_pool;
        right_count = r1_pool;
    
        if(left_count >= SYMBOL_COUNT || right_count >= SYMBOL_COUNT)
            update_likelihood(left_count, right_count, left_divider, right_divider);
    
        left_count = l2_pool;
        right_count = r2_pool;
    
        if(left_count >= SYMBOL_COUNT || right_count >= SYMBOL_COUNT)
            update_likelihood(left_count, right_count, left_divider, right_divider);   
    }
};*/

/*void rtiplus::update_score(state_merger *merger, apta_node* left, apta_node* right){
    rtiplus_data* l = (rtiplus_data*) left->data;
    rtiplus_data* r = (rtiplus_data*) right->data;
    
    //loglikelihood_merged = 0;
    //loglikelihood_orig = 0;
    //extra_parameters = 0;

    cerr << "adding left: " << left->size << endl;
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
    }
    add_split_likelihood(l);
    add_split_likelihood(r);
    
    add_parameters(l);
    add_parameters(r);

    //cerr << loglikelihood_orig << " " << loglikelihood_merged << " " << extra_parameters << endl;
    
    return;
};*/

/*void rtiplus::update_score_after(state_merger *merger, apta_node* left, apta_node* right){
    rtiplus_data* l = (rtiplus_data*) left->data;

    cerr << "adding merged: " << left->size << endl;
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
    add_merged_likelihood(l);
    del_parameters(l);

    //cerr << loglikelihood_orig << " " << loglikelihood_merged << " " << extra_parameters << endl;    
    //compute_score(merger, left, right);
   
    return;
};*/

/*void likelihoodratio::print_dot(iostream& output, state_merger* merger){
    alergia::print_dot(output, merger);
};*/

/*void rtiplus::add_split_likelihood(double count, double divider){
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
    for(type_num_map::iterator it = l->trans_counts.begin(); it != l->trans_counts.end(); ++it){
        num_map& nm = (*it).second;
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            count = (*it2).second;
            if(count >= SYMBOL_COUNT && count > 0) extra_parameters++;
            else pool += count;
        }
    }
    if(pool != 0) extra_parameters++;

    for(int i = 0; i < inputdata::num_attributes; ++i) {
        pool = 0;
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            count  = l->quantile_counts[i][j];
            if(count >= SYMBOL_COUNT && count > 0) extra_parameters++;
            else pool += count;
        }
        if(pool != 0) extra_parameters++;
    }
};

void rtiplus::del_parameters(rtiplus_data* l){
    double count = 0.0;
    double pool = 0.0;
    for(type_num_map::iterator it = l->trans_counts.begin(); it != l->trans_counts.end(); ++it){
        num_map& nm = (*it).second;
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            count = (*it2).second;
            if(count >= SYMBOL_COUNT && count > 0) extra_parameters--;
            else pool += count;
        }
    }
    if(pool != 0) extra_parameters--;
    
    for(int i = 0; i < inputdata::num_attributes; ++i) {
        pool = 0;
        for(int j = 0; j < rtiplus::attribute_quantiles[i].size() + 1; ++j){
            count  = l->quantile_counts[i][j];
            if(count >= SYMBOL_COUNT && count > 0) extra_parameters--;
            else pool += count;
        }
        if(pool != 0) extra_parameters--;
    }
};

void rtiplus::add_merged_likelihood(rtiplus_data* l){
    for(type_num_map::iterator it = l->trans_counts.begin(); it != l->trans_counts.end(); ++it){
        num_map& nm = (*it).second;
        
        
        //cerr << (*it).first << endl;

        double divider = 0.0;
        double count = 0.0;
        double pool = 0.0;

        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            count = (*it2).second;
            divider += count + CORRECTION;
        }
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            count = (*it2).second;
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
        return;
    }
};

void rtiplus::del_merged_likelihood(rtiplus_data* l){
    for(type_num_map::iterator it = l->trans_counts.begin(); it != l->trans_counts.end(); ++it){
        num_map& nm = (*it).second;

        double divider = 0.0;
        double count = 0.0;
        double pool = 0.0;
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            count = (*it2).second;
            divider += count + CORRECTION;
        }
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            count = (*it2).second;
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
        return;
    }
};

void rtiplus::add_split_likelihood(rtiplus_data* l){
    double divider = 0.0;
    double count = 0.0;
    double pool = 0.0;
    for(num_map::iterator it = l->pos_begin(); it != l->pos_end(); ++it){
        divider += (*it).second + CORRECTION;
    }
    for(num_map::iterator it = l->pos_begin(); it != l->pos_end(); ++it){
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
    for(num_map::iterator it = l->pos_begin(); it != l->pos_end(); ++it){
        divider += (*it).second + CORRECTION;
    }
    for(num_map::iterator it = l->pos_begin(); it != l->pos_end(); ++it){
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
};*/
