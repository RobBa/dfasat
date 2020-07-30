#include "state_merger.h"
#include "evaluate.h"
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <stdio.h>
//#include <gsl/gsl_cdf.h>

#include "overlap-driven.h"
#include "parameters.h"

REGISTER_DEF_DATATYPE(overlap_data);
REGISTER_DEF_TYPE(overlap_driven);
//DerivedRegister<overlap_driven> overlap_driven::reg("overlap_driven");

void overlap_data::print_transition_label(iostream& output, int symbol, apta* apta_context){
    output << pos(symbol) << " ";
};

void overlap_data::print_state_label(iostream& output, apta* aptacontext){
    count_data::print_state_label(output, aptacontext);
    output << "\n" << num_paths() << " " << num_final();
};

/* Overlap driven, count overlap in positive transitions, used in Stamina winner */
bool overlap_driven::consistent(state_merger *merger, apta_node* left, apta_node* right){
    if(count_driven::consistent(merger, left, right) == false){
        inconsistency_found = true;
        return false;
    }
    
    //return true;
    
    
    /*if(merger->aut->root == left){// || left->label != right->label){
        inconsistency_found = true;
        return false;
    }*/
    //return true;

    overlap_data* l = (overlap_data*) left->data;
    overlap_data* r = (overlap_data*) right->data;

    /*if(l->pos_final() != 0 && r->pos_final() == 0){
        inconsistency_found = true; return false;
    }
    if(r->pos_final() != 0 && l->pos_final() == 0){
        inconsistency_found = true; return false;
    }*/

    //cerr << l->pos_paths() << " " << r->pos_paths() << endl;

    if(l->pos_paths() >= STATE_COUNT){
        for(num_map::iterator it = r->pos_begin(); it != r->pos_end(); ++it){
            if(it->second >= SYMBOL_COUNT && l->pos(it->first) == 0){
                inconsistency_found = true;
                return false;
            }
        }
        /*if(l->pos_final() >= SYMBOL_COUNT && r->pos_final() == 0){
            inconsistency_found = true;
            return false;
        }*/
        /*for(num_map::iterator it = r->num_neg.begin(); it != r->num_neg.end(); ++it){
            if(it->second >= SYMBOL_COUNT & l->neg(it->first) == 0){
                inconsistency_found = true;
                return false;
            }
        }*/
    }

    if(r->pos_paths() >= STATE_COUNT){
        for(num_map::iterator it = l->pos_begin(); it != l->pos_end(); ++it){
            if(it->second >= SYMBOL_COUNT && r->pos(it->first) == 0){
                inconsistency_found = true;
                return false;
            }
        }
        /*if(l->pos_final() == 0 && r->pos_final() >= SYMBOL_COUNT){
            inconsistency_found = true;
            return false;
        }*/
        /*for(num_map::iterator it = l->num_neg.begin();it != l->num_neg.end(); ++it){
            if(it->second >= SYMBOL_COUNT & r->neg(it->first) == 0){
                inconsistency_found = true;
                return false;
            }
        }*/
    }
    return true;
};

void overlap_driven::update_score(state_merger *merger, apta_node* left, apta_node* right){
    overlap_data* l = (overlap_data*) left->data;
    overlap_data* r = (overlap_data*) right->data;

    //cerr << "update score" << endl;

    if (inconsistency_found) return;
    if (consistent(merger, left, right) == false) return;
    
    double num_matched = 0;
    double num_unmatched = 0;
    
    for(int i = 0; i < inputdata::alphabet.size(); ++i){
        //cerr << l->pos(i) << " " << r->pos(i) << endl;
        if(l->pos(i) != 0 && r->pos(i) != 0){
            num_matched++;
            //if(l->pos(i) > r->pos(i)) overlap += r->pos(i);
            //if(r->pos(i) > l->pos(i)) overlap += l->pos(i);
        } else {
            if(l->pos(i) != 0){
                num_unmatched++;//overlap -= l->pos(i);
            }
            if(r->pos(i) != 0){
                num_unmatched++;//overlap -= r->pos(i);
            }
        }
        /*if(l->pos(i) != 0 && r->pos(i) != 0){
            overlap += 1;
        }*/
        /*if(l->neg(i) != 0 && r->neg(i) != 0){
            overlap += 1;
        }*/
    }
    overlap += 1;//num_matched;// - num_unmatched;
    return;
    
    if(l->pos_final() != 0 && r->pos_final() != 0){
        if(l->pos_final() > r->pos_final()) overlap += r->pos_final();
        if(r->pos_final() > l->pos_final()) overlap += l->pos_final();
    } else {
        if(l->pos_final() != 0){
            overlap -= l->pos_final();
        }
        if(r->pos_final() != 0){
            overlap -= r->pos_final();
        }
    }
    /*if(l->pos_final() != 0 && r->pos_final() != 0){
        overlap += 1;
    }*/
};


double overlap_driven::compute_score(state_merger *merger, apta_node* left, apta_node* right){
    //cerr << overlap << endl;
  if(overlap > 0) return (int) (overlap);
  return 0;
};

void overlap_driven::reset(state_merger *merger){
  inconsistency_found = false;
  overlap = 0;
};

void overlap_data::print_transition_label(iostream& output, int symbol){
    output << (pos(symbol));
};



