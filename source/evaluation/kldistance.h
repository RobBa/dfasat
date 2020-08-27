#ifndef __KLDISTANCE__
#define __KLDISTANCE__

#include "alergia.h"

typedef map<int, float> prob_map;
typedef map<int, prob_map> type_prob_map;

/* The data contained in every node of the prefix tree or DFA */
class kl_data: public alergia_data {
protected:
  REGISTER_DEC_DATATYPE(kl_data);
public:

    type_prob_map original_probability_count;
    
    inline double opc(int type, int symbol){
        type_prob_map::iterator it = original_probability_count.find(type);
        if(it == original_probability_count.end()) return 0;
        prob_map::iterator it2 = it->second.find(symbol);
        if(it2 == it->second.end()) return 0;
        return it2->second;
    };
    
    virtual void update(evaluation_data* right);
    virtual void undo(evaluation_data* right);
};

class kldistance: public alergia {

protected:
  REGISTER_DEC_TYPE(kldistance);

  void update_perplexity(apta_node*,float,float,float,float,float,float);

public:
  double perplexity;
  int extra_parameters;

  virtual bool consistent(state_merger *merger, apta_node* left, apta_node* right);
  virtual void update_score(state_merger *merger, apta_node* left, apta_node* right);
  virtual bool compute_consistency(state_merger *merger, apta_node* left, apta_node* right);
  virtual double compute_score(state_merger*, apta_node* left, apta_node* right);
  virtual void reset(state_merger *merger);
  //virtual void print_dot(iostream&, state_merger *);
  virtual void initialize(state_merger *);
};

#endif
