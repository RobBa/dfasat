#include "state_merger.h"
#include "tsp.h"
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <stdio.h>
#include <gsl/gsl_cdf.h>
#include "parameters.h"

REGISTER_DEF_DATATYPE(tsp_data);
REGISTER_DEF_TYPE(tsp);

void tsp_data::print_state_label(iostream& output, apta* apta_context){
    /*
     for(set<int>::iterator it = done_tasks.begin(); it != done_tasks.end(); ++it){
        output << apta_context->alph_str(*it) << " ";
    }
    output << "\n";
    for(set<int>::iterator it = future_tasks.begin(); it != future_tasks.end(); ++it){
        output << apta_context->alph_str(*it) << " ";
    }
    */
    //output << "\n";
    output << "past: " << exp_pastlength << "\n" << "exp length: " << exp_futurelength << "\n" << "min length: " << min_futurelength;
};

tsp_data::tsp_data(){
    min_futurelength = -1.0;
    undo_futurelength = -1.0;
    min_pastlength = -1.0;
    undo_pastlength = -1.0;
    exp_pastlength = 0.0;
    exp_futurelength = 0.0;
    num_past = 0;
    num_future = 0;
};

void tsp_data::read_from(tail* t){
    while(t->past() != 0){
        t = t->past();
    }
    double length = std::stod(inputdata::get_data(t));
    while(t->future() != 0){
        t = t->future();
        length = length + std::stod(inputdata::get_data(t));
    }
    if(min_futurelength == -1.0 || min_futurelength > length){
        min_futurelength = length;
    }
    exp_futurelength = (double)((exp_futurelength * ((double)num_future)) + length) / ((double)(num_future + 1));
    num_future = num_future + 1;
};

void tsp_data::read_to(tail* t){
    double length = std::stod(inputdata::get_data(t));
    while(t->past() != 0){
        t = t->past();
        length = length + std::stod(inputdata::get_data(t));
    }
    if(min_pastlength == -1.0 || min_pastlength > length){
        min_pastlength = length;
    }
    exp_pastlength = (double)((exp_pastlength * ((double)num_past)) + length) / ((double)(num_past + 1));
    num_past = num_past + 1;
};

void tsp_data::update(evaluation_data* right){
    tsp_data* r = (tsp_data*) right;
    
    if(min_futurelength > r->min_futurelength){
        r->undo_futurelength = min_futurelength;
        min_futurelength = r->min_futurelength;
    }

    if(min_pastlength > r->min_pastlength){
        r->undo_pastlength = min_pastlength;
        min_pastlength = r->min_pastlength;
    }
    
    if(num_past + r->num_past != 0){
        exp_pastlength = ((exp_pastlength * ((double)num_past) + (r->exp_pastlength * ((double)r->num_past)))) / ((double)num_past + (double)r->num_past);
    }
    if(num_future + r->num_future != 0){
        exp_futurelength = ((exp_futurelength * ((double)num_future) + (r->exp_futurelength * ((double)r->num_future)))) / ((double)num_future + (double)r->num_future);
    }
    
    num_past = num_past + r->num_past;
    num_future = num_future + r->num_future;

    set<int>::iterator it  = future_tasks.begin();
    set<int>::iterator it2 = r->future_tasks.begin();
    r->undo_info.clear();
    
    while(it != future_tasks.end() && it2 != r->future_tasks.end()){
        if(*it == *it2){
            it++;
            it2++;
        } else if(*it < *it2){
            it++;
        } else {
            r->undo_info.insert(*it2);
            it2++;
        }
    }
    while(it2 != r->future_tasks.end()){
        r->undo_info.insert(*it2);
        it2++;
    }
    
    for(it = r->undo_info.begin(); it != r->undo_info.end(); ++it){
        future_tasks.insert(*it);
    }

    it  = done_tasks.begin();
    it2 = r->done_tasks.begin();
    r->undo_dinfo.clear();
    
    while(it != done_tasks.end() && it2 != r->done_tasks.end()){
        if(*it == *it2){
            it++;
            it2++;
        } else if(*it < *it2){
            it++;
        } else {
            r->undo_dinfo.insert(*it2);
            it2++;
        }
    }
    while(it2 != r->done_tasks.end()){
        r->undo_dinfo.insert(*it2);
        it2++;
    }
    
    for(it = r->undo_dinfo.begin(); it != r->undo_dinfo.end(); ++it){
        done_tasks.insert(*it);
    }
};

void tsp_data::undo(evaluation_data* right){
    tsp_data* r = (tsp_data*) right;

    if(r->undo_futurelength != -1.0){
        min_futurelength = r->undo_futurelength;
        r->undo_futurelength = -1.0;
    }
    if(r->undo_pastlength != -1.0){
        min_pastlength = r->undo_pastlength;
        r->undo_pastlength = -1.0;
    }
    
    num_past = num_past - r->num_past;
    num_future = num_future - r->num_future;

    if(num_past != 0)
        exp_pastlength = ((exp_pastlength * ((double)num_past + (double)r->num_past)) - (r->exp_pastlength * ((double)r->num_past))) / ((double)num_past);
    else
        exp_pastlength = 0.0;

    if(num_future != 0)
        exp_futurelength = ((exp_futurelength * ((double)num_future + (double)r->num_future)) - (r->exp_futurelength * ((double)r->num_future))) / ((double)num_future);
    else
        exp_futurelength = 0.0;
    
    for(set<int>::iterator it = r->undo_info.begin(); it != r->undo_info.end(); ++it){
        future_tasks.erase(*it);
    }
    for(set<int>::iterator it = r->undo_dinfo.begin(); it != r->undo_dinfo.end(); ++it){
        done_tasks.erase(*it);
    }
};

bool tsp::consistent(state_merger *merger, apta_node* left, apta_node* right){
    if(left->depth != right->depth){ inconsistency_found = true; return false; }

    if(merger->aut->root == left){// || left->label != right->label){
        inconsistency_found = true;
        return false;
    }

    tsp_data* l = (tsp_data*) left->data;
    tsp_data* r = (tsp_data*) right->data;
    
    //if((l->min_futurelength < 80 || r->min_futurelength < 80) && l->min_futurelength - r->min_futurelength > 10){ inconsistency_found = true; return false; }
    //if((l->min_futurelength < 80 || r->min_futurelength < 80) && l->min_futurelength - r->min_futurelength <-10){ inconsistency_found = true; return false; }
    
    //return true;

    if(left->size > 5 && right->size > 5 && l->exp_futurelength - r->exp_futurelength > 20){ inconsistency_found = true; return false; }
    if(left->size > 5 && right->size > 5 && l->exp_futurelength - r->exp_futurelength <-20){ inconsistency_found = true; return false; }
//    if(l->num_past > 10 && r->num_past > 10 && l->exp_pastlength - r->exp_pastlength > 10){ inconsistency_found = true; return false; }
//    if(l->num_past > 10 && r->num_past > 10 && l->exp_pastlength - r->exp_pastlength <-10){ inconsistency_found = true; return false; }
    //if(l->min_pastlength - l->min_pastlength > 20){ inconsistency_found = true; return false; }
    //if(r->min_pastlength - l->min_pastlength > 20){ inconsistency_found = true; return false; }

    return true;

    if(l->min_pastlength - r->min_pastlength >= 0){
        for(set<int>::iterator it = l->done_tasks.begin(); it != l->done_tasks.end(); ++it){
            if(r->done_tasks.find(*it) == r->done_tasks.end()){
                inconsistency_found = true;
                return false;
            }
        }
    }
    if(r->min_pastlength - l->min_pastlength >= 0){
        for(set<int>::iterator it = r->done_tasks.begin(); it != r->done_tasks.end(); ++it){
            if(l->done_tasks.find(*it) == l->done_tasks.end()){
                inconsistency_found = true;
                return false;
            }
        }
    }
    
    return true;
    
    if(l->min_futurelength - r->min_futurelength > 20){ inconsistency_found = true; return false; }
    if(r->min_futurelength - l->min_futurelength > 20){ inconsistency_found = true; return false; }
    //if(l->min_pastlength - l->min_pastlength > 20){ inconsistency_found = true; return false; }
    //if(r->min_pastlength - l->min_pastlength > 20){ inconsistency_found = true; return false; }

    return true;
};

void tsp::update_score(state_merger *merger, apta_node* left, apta_node* right){
    tsp_data* l = (tsp_data*) left->data;
    tsp_data* r = (tsp_data*) right->data;

    if(l->exp_pastlength - r->exp_pastlength < 20 && r->exp_pastlength - l->exp_pastlength >-20){ total_merges++; }
    return;
    
    double diff = l->exp_futurelength - r->exp_futurelength;
    if(diff < 0) diff = r->exp_futurelength - l->exp_futurelength;

    total_merges += 50.0 - diff;
};

double tsp::compute_score(state_merger *merger, apta_node* left, apta_node* right){
    return total_merges;
};

void tsp::reset(state_merger *merger ){
    inconsistency_found = false;
    total_merges = 0;
};

void tsp::initialize(state_merger* m){
    for(merged_APTA_iterator it = merged_APTA_iterator(m->aut->root); *it != 0; ++it){
        apta_node* n = *it;
        tsp_data* d = (tsp_data*)n->data;
        d->min_futurelength = -1.0;
        d->exp_futurelength = 0.0;
        d->num_future = 0;
        for(tail_iterator it = tail_iterator(n); *it != 0; ++it){
            tail* t = *it;
            d->read_from(t);
        }
    }

    for(merged_APTA_iterator it = merged_APTA_iterator(m->aut->root); *it != 0; ++it){
        apta_node* n = *it;
        apta_node* p = n;
        tsp_data* data = (tsp_data*)n->data;
        while(p->source != 0){
            data->done_tasks.insert(p->label);
            p = p->source;
        }
    }

    for(merged_APTA_iterator it = merged_APTA_iterator(m->aut->root); *it != 0; ++it){
        apta_node* n = *it;
        apta_node* p = n;
        while(p->source != 0){
            tsp_data* data = (tsp_data*)p->source->data;
            data->future_tasks.insert(n->label);
            p = p->source;
        }
    }
};


