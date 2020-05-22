#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <stdio.h>
#include <gsl/gsl_cdf.h>

#include "state_merger.h"
#include "evaluate.h"
#include "kldistance.h"

#include "parameters.h"

REGISTER_DEF_DATATYPE(kl_data);
REGISTER_DEF_TYPE(kldistance);

void kl_data::update(evaluation_data* right){
    alergia_data::update(right);
    kl_data* other = (kl_data*)right;
    for(type_prob_map::iterator it = other->original_probability_count.begin(); it != other->original_probability_count.end(); ++it){
        for(prob_map::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2){
            original_probability_count[it->first][it2->first] = opc(it->first, it2->first) + it2->second;
        }
    }
};

void kl_data::undo(evaluation_data* right){
    alergia_data::undo(right);
    kl_data* other = (kl_data*)right;
    for(type_prob_map::iterator it = other->original_probability_count.begin(); it != other->original_probability_count.end(); ++it){
        for(prob_map::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2){
            original_probability_count[it->first][it2->first] = opc(it->first, it2->first) - it2->second;
        }
    }
};

bool kldistance::consistent(state_merger *merger, apta_node* left, apta_node* right){
    return count_driven::consistent(merger, left, right);
};

void kldistance::update_perplexity(apta_node* left, float op_count_left, float op_count_right, float count_left, float count_right, float left_divider, float right_divider){
    //cerr << "\t" << op_count_left << " " << op_count_right << " " << count_left << " " << count_right << endl;
    //cerr << "\t" << count_left / left_divider << " " << count_right / right_divider << endl;
    if(already_merged(left) == false){
        if(count_left != 0){
            perplexity += op_count_left * log(count_left / left_divider);
            perplexity -= op_count_left * log((count_left + count_right) / (left_divider + right_divider));
        }
        if(count_right != 0){
            perplexity += op_count_right * log(count_right / right_divider);
            perplexity -= op_count_right * log((count_left + count_right) / (left_divider + right_divider));
        }
    } else {
        if(count_left != 0){
            perplexity += op_count_left * log(count_left / left_divider);
            perplexity -= op_count_left * log((count_left + count_right) / (left_divider + right_divider));
        }
        if(count_right != 0){
            perplexity += op_count_right * log(count_right / right_divider);
            perplexity -= op_count_right * log((count_left + count_right) / (left_divider + right_divider));
        }
    }
    if(count_left > 0.0 && count_right > 0.0) extra_parameters = extra_parameters + 1;
    //cerr << perplexity << endl;
};

/* Kullback-Leibler divergence (KL), MDI-like, computes the KL value/extra parameters and uses it as score and consistency */
void kldistance::update_score(state_merger *merger, apta_node* left, apta_node* right){
    evaluation_function::update_score(merger, left, right);
    kl_data* l = (kl_data*) left->data;
    kl_data* r = (kl_data*) right->data;

    if(r->pos_paths() < 0 || l->pos_paths() < 0) return;

    float left_divider = (float)l->pos_paths();
    float right_divider = (float)r->pos_paths();
    
    for(type_num_map::iterator it = l->counts_begin(); it != l->counts_end(); it++){
        int type = it->first;
        num_map& nm = it->second;
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            int symbol = it2->first;
            int left_count = it2->second;
            int right_count = r->count(type, symbol);
            update_perplexity(left, l->opc(type, symbol), r->opc(type, symbol), left_count, right_count, left_divider, right_divider);
        }        
    }    
};

bool kldistance::compute_consistency(state_merger *merger, apta_node* left, apta_node* right){
  if (inconsistency_found) return false;
  if (extra_parameters == 0) return false;

  if ((perplexity / (float)extra_parameters) > CHECK_PARAMETER) return false;

  return true;
};

double kldistance::compute_score(state_merger *merger, apta_node* left, apta_node* right){
  if (inconsistency_found == true) return false;
  if (extra_parameters == 0) return -1;

  double val = (perplexity / (double)extra_parameters);
  
  return 100000 - (int)(val * 100.0);
};

void kldistance::reset(state_merger *merger){
  alergia::reset(merger);
  inconsistency_found = false;
  perplexity = 0;
  extra_parameters = 0;
};

/*void kldistance::print_dot(iostream& output, state_merger* merger){
    count_driven::print_dot(output, merger);
};*/

void kldistance::initialize(state_merger* merger){
    for(merged_APTA_iterator Ait = merged_APTA_iterator(merger->aut->root); *Ait != 0; ++Ait){
        apta_node* node = *Ait;
        kl_data* l = (kl_data*) node->data;
        
        if(l->pos_paths() < STATE_COUNT) continue;

        for(type_num_map::iterator it = l->counts_begin(); it != l->counts_end(); it++){
            int type = it->first;
            num_map& nm = it->second;
            for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
                int symbol = it2->first;
                float count = it2->second;
                l->original_probability_count[type][symbol] = count * (count / (float)l->pos_paths());
            }
        }
    }
};

