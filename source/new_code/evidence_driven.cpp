#include "evidence_driven.h"

/* Evidence driven state merging, count number of pos-pos and neg-neg merges */
void evidence_driven::update_score(state_merger *merger, apta_node* left, apta_node* right){
  if(left->pos_final() > 0 && right->pos_final() > 0) num_pos += 1;
  if(left->neg_final() > 0 && right->neg_final() > 0) num_neg += 1;
};

double evidence_driven::compute_score(state_merger *merger, apta_node* left, apta_node* right){
  return num_pos + num_neg;
};

void evidence_driven::reset(state_merger *merger){
  inconsistency_found = false;
  num_pos = 0;
  num_neg = 0;
};