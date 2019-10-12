#ifndef __NINO__
#define __NINO__

typedef neighbor_map map< tail*, pair< tail*, double >  >;

/* The data contained in every node of the prefix tree or DFA */
class nino_data: public evaluation_data {
protected:
  REGISTER_DEC_DATATYPE(nino_data);

public:

    neighbor_map neighbors;
    
    nino_data();
    
    virtual void read_from(tail* t);
    virtual void update(evaluation_data* right);
    virtual void undo(evaluation_data* right);
    
    virtual void split_update_single(evaluation_data* other, tail* t);
    virtual void split_update(evaluation_data* other);
    virtual void split_undo(evaluation_data* other);
};

class nino: public evaluation_function {

protected:
  REGISTER_DEC_TYPE(nino);
  
public:
    
  virtual bool consistent(state_merger *merger, apta_node* left, apta_node* right);
  virtual void update_score(state_merger *merger, apta_node* left, apta_node* right);
  virtual void update_score_after(state_merger *merger, apta_node* left, apta_node* right);

  virtual bool compute_consistency(state_merger *, apta_node* left, apta_node* right);
  virtual double compute_score(state_merger *, apta_node* left, apta_node* right);

  virtual bool split_consistent(state_merger*, apta_node* left, apta_node* right);
  virtual void split_update_score_before(state_merger*, apta_node* left, apta_node* right, tail* t);
  virtual void split_update_score_after(state_merger*, apta_node* left, apta_node* right, tail* t);
    
  virtual bool split_compute_consistency(state_merger *, apta_node* left, apta_node* right);
  virtual double split_compute_score(state_merger *, apta_node* left, apta_node* right);
    
  virtual void reset(state_merger *);
  virtual void reset_split(state_merger *, apta_node *);
};

#endif
