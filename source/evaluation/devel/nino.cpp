#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <stdio.h>
#include <cmath>
#include <gsl/gsl_cdf.h>

#include "state_merger.h"
#include "evaluate.h"
#include "nino.h"
#include "parameters.h"

REGISTER_DEF_DATATYPE(nino_data);
REGISTER_DEF_TYPE(nino);

nino_data::nino_data() {
};

void nino_data::update(evaluation_data* right){
};

void nino_data::undo(evaluation_data* right){
};

void nino_data::split_update(evaluation_data* right){
    undo(right);
};

void nino_data::split_undo(evaluation_data* right){
    update(right);
};

void nino_data::split_update_single(evaluation_data* right, tail* t){
};

bool nino::consistent(state_merger *merger, apta_node* left, apta_node* right){
    return likelihoodratio::consistent(merger, left, right);
};

/* Likelihood Ratio (LR), computes an LR-test (used in RTI) and uses the p-value as score and consistency */

bool nino::split_consistent(state_merger*, apta_node* left, apta_node* right){
    return true;
};

void nino::update_score(state_merger *merger, apta_node* left, apta_node* right){
    nino_data* l = (nino_data*) left->data;
    nino_data* r = (nino_data*) right->data;
    
    return;
};

void nino::update_score_after(state_merger *merger, apta_node* left, apta_node* right){
    nino_data* l = (nino_data*) left->data;
    nino_data* r = (nino_data*) right->data;

    return;
};

void nino::split_update_score_before(state_merger*, apta_node* left, apta_node* right, tail* t){
    nino_data* l = (nino_data*) left->data;
    nino_data* r = (nino_data*) right->data;
    
};

void nino::split_update_score_after(state_merger*, apta_node* left, apta_node* right, tail* t){
    nino_data* l = (nino_data*) left->data;
    nino_data* r = (nino_data*) right->data;
    
};

bool nino::compute_consistency(state_merger *, apta_node* left, apta_node* right){
    nino_data* l = (nino_data*) left->data;
    nino_data* r = (nino_data*) right->data;
    
    return true;
};

bool nino::compute_score(state_merger *, apta_node* left, apta_node* right){
    nino_data* l = (nino_data*) left->data;
    nino_data* r = (nino_data*) right->data;
    
    return true;
};

double nino::split_compute_consistency(state_merger *, apta_node* left, apta_node* right){
    nino_data* l = (nino_data*) left->data;
    nino_data* r = (nino_data*) right->data;
    
    return true;
};

double nino::split_compute_score(state_merger *, apta_node* left, apta_node* right){
    nino_data* l = (nino_data*) left->data;
    nino_data* r = (nino_data*) right->data;
    
    return true;
};
