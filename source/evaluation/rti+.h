#ifndef __RTIPLUS__
#define __RTIPLUS__

#include "likelihood.h"

typedef vector< vector<int> > quantile_map;

/* The data contained in every node of the prefix tree or DFA */
class rtiplus_data: public likelihood_data {
protected:
  REGISTER_DEC_DATATYPE(rtiplus_data);

public:

    quantile_map quantile_counts;
    double loglikelihood;
    
    rtiplus_data();
    
    virtual void read_from(tail* t);
    virtual void update(evaluation_data* right);
    virtual void undo(evaluation_data* right);
    
    virtual void split_update_single(evaluation_data* other, tail* t);
    virtual void split_update(evaluation_data* other);
    virtual void split_undo(evaluation_data* other);
    
    virtual void del_tail(tail* t);

    virtual void set_loglikelihood();
    
    int num_parameters();

    virtual void print_state_label(iostream& output, apta* aptacontext);
};

class rtiplus: public likelihoodratio {

protected:
  REGISTER_DEC_TYPE(rtiplus);
  
  rtiplus_data temp_data;
  
public:
    
  static vector< vector<double> > attribute_quantiles;
    
    //void add_merged_likelihood(double count, double divider);
    //void del_merged_likelihood(double count, double divider);
    //void add_split_likelihood(double count, double divider);
    //void del_split_likelihood(double count, double divider);
    
  //void del_merged_likelihood(rtiplus_data* l);
  //void add_merged_likelihood(rtiplus_data* l);
  //void del_split_likelihood(rtiplus_data* l);
  //void add_split_likelihood(rtiplus_data* l);

    //void del_parameters(rtiplus_data* l);
    //void add_parameters(rtiplus_data* l);

  //virtual bool consistent(state_merger *merger, apta_node* left, apta_node* right);
  virtual void update_score(state_merger *merger, apta_node* left, apta_node* right);
  virtual void update_score_after(state_merger *merger, apta_node* left, apta_node* right);

  //virtual bool split_consistent(state_merger*, apta_node* left, apta_node* right);
  virtual void split_update_score_before(state_merger*, apta_node* left, apta_node* right, tail* t);
  virtual void split_update_score_after(state_merger*, apta_node* left, apta_node* right, tail* t);
    
  virtual bool split_compute_consistency(state_merger *, apta_node* left, apta_node* right);
  virtual double split_compute_score(state_merger *, apta_node* left, apta_node* right);
    
  virtual void initialize(state_merger* m);
  virtual void initialize_globals();
  virtual void reset_split(state_merger *, apta_node *);
};

#endif
