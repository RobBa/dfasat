#ifndef __COUNT__
#define __COUNT__

#include "evaluate.h"

typedef map<int, int> num_map;

/* The data contained in every node of the prefix tree or DFA */
class count_data: public evaluation_data {

protected:
  REGISTER_DEC_DATATYPE(count_data);

public:
    //int pos_final();
    //int neg_final();
    //int pos_paths();
    //int neg_paths();
    
    num_map final_counts;
    num_map path_counts;
    
    int total_final;
    int total_paths;

    count_data();
    
    virtual void read_from(int type, int index, int length, int symbol, string data);
    virtual void read_to(int type, int index, int length, int symbol, string data);

    virtual void print_transition_label(iostream& output, int symbol, apta* aptacontext);
    virtual void print_state_label(iostream& output, apta* aptacontext);
 
    virtual void update(evaluation_data* right);
    virtual void undo(evaluation_data* right);

    virtual int sink_type(apta_node* node);
    virtual bool sink_consistent(apta_node* node, int type);
    virtual int num_sink_types();
    
    inline int num_paths(){
        return total_paths;
    }
    
    inline int num_final(){
        return total_final;
    }
    
    inline int pos_paths(){
        return path_counts[1];
    }

    inline int neg_paths(){
        return path_counts[0];
    }

    inline int pos_final(){
        return final_counts[1];
    }

    inline int neg_final(){
        return final_counts[0];
    }
};

class count_driven: public evaluation_function {

protected:
  REGISTER_DEC_TYPE(count_driven);

public:
  int num_merges;

  virtual void update_score(state_merger *merger, apta_node* left, apta_node* right);
  virtual double  compute_score(state_merger*, apta_node* left, apta_node* right);
  virtual void reset(state_merger *merger);
  //virtual bool consistency_check(evaluation_data* l, evaluation_data* r);
  virtual bool consistent(state_merger *merger, apta_node* left, apta_node* right);

  //virtual void print_dot(iostream&, state_merger *);
};

#endif
