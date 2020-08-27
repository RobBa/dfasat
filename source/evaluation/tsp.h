#ifndef __MSEERROR__
#define __MSEERROR__

#include "evaluate.h"

typedef list<double> double_list;

/* The data contained in every node of the prefix tree or DFA */
class tsp_data: public evaluation_data {

protected:
  REGISTER_DEC_DATATYPE(tsp_data);
  
public:
    set<int> done_tasks;
    set<int> future_tasks;
    set<int> undo_info;
    set<int> undo_dinfo;
    
    double exp_pastlength;
    double exp_futurelength;
    int num_future;
    int num_past;

    double min_pastlength;
    double min_futurelength;
    double undo_pastlength;
    double undo_futurelength;

    tsp_data();

    virtual void print_state_label(iostream& output, apta* apta_context);

    virtual void read_from(tail* t);
    virtual void read_to(tail* t);
    virtual void update(evaluation_data* right);
    virtual void undo(evaluation_data* right);
};

class tsp: public evaluation_function{

protected:
  REGISTER_DEC_TYPE(tsp);

public:
  double total_merges = 0;
  
  virtual bool consistent(state_merger *merger, apta_node* left, apta_node* right);
  virtual void update_score(state_merger *merger, apta_node* left, apta_node* right);
  virtual double compute_score(state_merger*, apta_node* left, apta_node* right);
  virtual void reset(state_merger *merger);
  virtual void initialize(state_merger *merger);
};

#endif
