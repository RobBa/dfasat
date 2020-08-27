#ifndef __LIKELIHOOD__
#define __LIKELIHOOD__

#include "alergia.h"

/* The data contained in every node of the prefix tree or DFA */
class likelihood_data: public alergia_data {
protected:
  REGISTER_DEC_DATATYPE(likelihood_data);
    virtual void initialize();
};

class likelihoodratio: public alergia {

protected:
  REGISTER_DEC_TYPE(likelihoodratio);
  

public:
  void update_likelihood(double,double,double,double);

  double loglikelihood_orig;
  double loglikelihood_merged;
  double extra_parameters;

  virtual bool consistent(state_merger *merger, apta_node* left, apta_node* right);
  virtual void update_score(state_merger *merger, apta_node* left, apta_node* right);
  virtual bool compute_consistency(state_merger *merger, apta_node* left, apta_node* right);
  virtual double  compute_score(state_merger*, apta_node* left, apta_node* right);
  virtual void reset(state_merger *merger);
  //virtual void print_dot(iostream&, state_merger *);
};

#endif
