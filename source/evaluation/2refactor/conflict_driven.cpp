#include "state_merger.h"
#include "evaluate.h"
#include "evaluation_factory.h"
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <stdio.h>
#include <gsl/gsl_cdf.h>


count_data::count_data(){
    pos_final() = 0;
    neg_final() = 0;
    pos_paths() = 0;
    neg_paths() = 0;
};

void count_data::update(evaluation_data* right){
    pos_final() += right->pos_final();
    neg_final() += right->neg_final();
    pos_paths() += right->pos_paths();
    neg_paths() += right->neg_paths();
};

void count_data::undo(evaluation_data* right){
    pos_final() -= right->pos_final();
    neg_final() -= right->neg_final();
    pos_paths() -= right->pos_paths();
    neg_paths() -= right->neg_paths();
};

void count_data::from_string(int type, int index, int length, int symbol, string data){
    if(type == 1){
        pos_paths()++;
        if(length == index+1){
            pos_final()++;
        }
    } else {
        neg_paths()++;
        if(length == index+1){
            neg_final()++;
        }
    }
};

string count_data::to_string(){
        int node_size = (float)(pos_paths() + neg_paths())/2.0;

        if(pos_final() != 0) fprintf(output,"[shape=ellipse style=\"filled\" color=\"green\" label=\"%i\"]", pos_final());
        else if(neg_final() != 0) fprintf(output,"[shape=ellipse style=\"filled\" color=\"red\" label=\"%i\"]", neg_final());
        else fprintf(output,"[shape=ellipse\"]");
};

            else {
                if(positive){
                    node->num_pos[c] = node->pos(c) + 1;
                    node->pos_paths()++;
                } else {
                    node->num_neg[c] = node->neg(c) + 1;
                    node->neg_paths()++;
                }
            }
            //node->input_output
            if(occ >= 0)
                node->occs.push_front(occ);
                //node->child(c)->occs.push_front(occ);
        if(positive) node->pos_final()++;
        else node->neg_final()++;

    num_pos = num_map();
    num_neg = num_map();

    //left->old_depth = left->depth;
    //left->depth = min(left->depth, right->depth);

    right->merge_point = left->conflicts.end();
    --(right->merge_point);
    left->conflicts.splice(left->conflicts.end(), right->conflicts);
    ++(right->merge_point);

    right->occ_merge_point = left->occs.end();
    --(right->occ_merge_point);
    left->occs.splice(left->occs.end(), right->occs);
    ++(right->occ_merge_point);

    for(num_map::iterator it = right->num_pos.begin();it != right->num_pos.end(); ++it){
        left->num_pos[(*it).first] = left->pos((*it).first) + (*it).second;
    }
    for(num_map::iterator it = right->num_neg.begin();it != right->num_neg.end(); ++it){
        left->num_neg[(*it).first] = left->neg((*it).first) + (*it).second;
    }

    //left->depth = left->old_depth;

    for(num_map::iterator it = right->num_pos.begin();it != right->num_pos.end(); ++it){
        left->num_pos[(*it).first] = left->pos((*it).first) - (*it).second;
    }
    for(num_map::iterator it = right->num_neg.begin();it != right->num_neg.end(); ++it){
        left->num_neg[(*it).first] = left->neg((*it).first) - (*it).second;
    }

    right->conflicts.splice(right->conflicts.begin(), left->conflicts, right->merge_point, left->conflicts.end());
    right->occs.splice(right->occs.begin(), left->occs, right->occ_merge_point, left->occs.end());


    inline int pos(int i){
        num_map::iterator it = num_pos.find(i);
        if(it == num_pos.end()) return 0;
        return (*it).second;
    }

    inline int neg(int i){
        num_map::iterator it = num_neg.find(i);
        if(it == num_neg.end()) return 0;
        return (*it).second;
    }

/* default evaluation, count number of performed merges */
bool evaluation_function::consistent(state_merger *merger, apta_node* left, apta_node* right){
  if(inconsistency_found) return false;
  
  if(left->pos_final() != 0 && right->neg_final() != 0){ inconsistency_found = true; return false; }
  if(left->neg_final() != 0 && right->pos_final() != 0){ inconsistency_found = true; return false; }
    
  return true;
};

void evaluation_function::update_score(state_merger *merger, apta_node* left, apta_node* right){
  num_merges += 1;
};

bool evaluation_function::compute_consistency(state_merger *merger, apta_node* left, apta_node* right){
  return inconsistency_found == false && compute_score(merger, left, right) > LOWER_BOUND;
};

int evaluation_function::compute_score(state_merger *merger, apta_node* left, apta_node* right){
  return num_merges;
};

void evaluation_function::reset(state_merger *merger){
  inconsistency_found = false;
  num_merges = 0;
};

void evaluation_function::update(state_merger *merger){
};

void evaluation_function::initialize(state_merger *merger){
    compute_before_merge = false;
};

