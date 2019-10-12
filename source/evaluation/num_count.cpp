#include "state_merger.h"
#include "evaluate.h"
#include "evaluation_factory.h"
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <set>
#include <stdio.h>
#include <gsl/gsl_cdf.h>

#include "parameters.h"
#include "num_count.h"

REGISTER_DEF_TYPE(count_driven);
REGISTER_DEF_DATATYPE(count_data);

set<int> types;

count_data::count_data() : evaluation_data() {
    total_paths = 0;
    total_final = 0;
};

void count_data::print_transition_label(iostream& output, int symbol, apta* aptacontext){
    for(num_map::iterator it = path_counts.begin(); it != path_counts.end(); ++it){
        output << (*it).first << ":" << (*it).second << " - ";
    }
};

void count_data::print_state_label(iostream& output, apta* aptacontext){
    for(num_map::iterator it = final_counts.begin(); it != final_counts.end(); ++it){
        output << (*it).first << ":" << (*it).second << " - ";
    }
};

void count_data::read_from(int type, int index, int length, int symbol, string data){
    if(path_counts.find(type) == path_counts.end()){
        path_counts[type] = 1;
    } else {
        path_counts[type]++;
    }
    total_paths++;
    types.insert(type);
};

void count_data::read_to(int type, int index, int length, int symbol, string data){
    if(index < length - 1) return;
    
    if(final_counts.find(type) == final_counts.end()){
        final_counts[type] = 1;
    } else {
        final_counts[type]++;
    }
    total_final++;
};

void count_data::update(evaluation_data* right){
    count_data* other = reinterpret_cast<count_data*>(right);
    for(num_map::iterator it = other->final_counts.begin(); it != other->final_counts.end(); ++it){
        int type = (*it).first;
        int count = (*it).second;
        if(final_counts.find(type) != final_counts.end()){
            final_counts[type] += count;
        } else {
            final_counts[type] = count;
        }
    }
    for(num_map::iterator it = other->path_counts.begin(); it != other->path_counts.end(); ++it){
        int type = (*it).first;
        int count = (*it).second;
        if(path_counts.find(type) != path_counts.end()){
            path_counts[type] += count;
        } else {
            path_counts[type] = count;
        }
    }
    total_paths += other->total_paths;
    total_final += other->total_final;
};

void count_data::undo(evaluation_data* right){
    count_data* other = reinterpret_cast<count_data*>(right);

    for(num_map::iterator it = other->final_counts.begin(); it != other->final_counts.end(); ++it){
        int type = (*it).first;
        int count = (*it).second;
        final_counts[type] -= count;
    }
    for(num_map::iterator it = other->path_counts.begin(); it != other->path_counts.end(); ++it){
        int type = (*it).first;
        int count = (*it).second;
        path_counts[type] -= count;
    }
    total_paths -= other->total_paths;
    total_final -= other->total_final;
};

/*bool count_driven::consistency_check(evaluation_data* left, evaluation_data* right){
    count_data* l = reinterpret_cast<count_data*>(left);
    count_data* r = reinterpret_cast<count_data*>(right);
    
    if(l->pos_final() != 0 && r->neg_final() != 0){ return false; }
    if(l->neg_final() != 0 && r->pos_final() != 0){ return false; }
    
    return true;
};*/

/* default evaluation, count number of performed merges */
bool count_driven::consistent(state_merger *merger, apta_node* left, apta_node* right){
    if(inconsistency_found) return false;
  
    count_data* l = (count_data*)left->data;
    count_data* r = (count_data*)right->data;

    //if(l->pos_final() != 0 && r->neg_final() != 0){ inconsistency_found = true; return false; }
    //if(l->neg_final() != 0 && r->pos_final() != 0){ inconsistency_found = true; return false; }
    
    for(num_map::iterator it = l->final_counts.begin(); it != l->final_counts.end(); ++it){
        int type = (*it).first;
        int count = (*it).second;
        if(count != 0){
            for(num_map::iterator it2 = r->final_counts.begin(); it2 != r->final_counts.end(); ++it2){
                int type2 = (*it2).first;
                int count2 = (*it2).second;
                if(count2 != 0 && type2 != type){
                    inconsistency_found = true;
                    return false;
                }
            }
        }
    }
    
    return true;
};

void count_driven::update_score(state_merger *merger, apta_node* left, apta_node* right){
  num_merges += 1;
};

double count_driven::compute_score(state_merger *merger, apta_node* left, apta_node* right){
  return num_merges;
};

void count_driven::reset(state_merger *merger){
  num_merges = 0;
  evaluation_function::reset(merger);
  compute_before_merge=false;
};


// sinks for evaluation data type

/*bool is_rejecting_sink(apta_node* node){
    count_data* d = reinterpret_cast<count_data*>(node->data);

    node = node->find();
    return d->pos_paths() == 0 && d->pos_final() == 0;
};*/

int count_data::sink_type(apta_node* node){
    if(!USE_SINKS) return -1;

    count_data* d = reinterpret_cast<count_data*>(node->data);
    
    int type = -1;
    for(num_map::iterator it = d->final_counts.begin(); it != d->final_counts.end(); ++it){
        int count = (*it).second;
        if(count != 0){
            if(type == -1) type = (*it).second;
            else return -1;
        }
    }
    //return d->neg_paths() == 0 && d->neg_final() == 0;
    return type;
};

bool count_data::sink_consistent(apta_node* node, int type){
    if(!USE_SINKS) return true;
    
    //if(type == 0) return is_rejecting_sink(node);
    //if(type == 1) return is_accepting_sink(node);
    //return true;
    
    return sink_type(node) == type;
};

int count_data::num_sink_types(){
    if(!USE_SINKS) return 0;
    
    // accepting or rejecting
    return types.size();
};

