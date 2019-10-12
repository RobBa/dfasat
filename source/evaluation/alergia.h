#ifndef __ALERGIA__
#define __ALERGIA__

#include "evaluate.h"
#include "num_count.h"

typedef map<int, num_map> type_num_map;

/* The data contained in every node of the prefix tree or DFA */
class alergia_data: public count_data {

protected:
  REGISTER_DEC_DATATYPE(alergia_data);

public:
    /* counts of positive and negative transition uses */
    type_num_map trans_counts;
    //num_map num_pos;
    //num_map num_neg;
    
    inline int count(int type, int symbol){
        type_num_map::iterator it = trans_counts.find(type);
        if(it != trans_counts.end()){
            num_map::iterator it2 = (*it).second.find(symbol);
            if(it2 == (*it).second.end()) return 0;
            return (*it2).second;
        }
        return 0;
    }

    inline int count_type(int type){
        type_num_map::iterator it = trans_counts.find(type);
        if(it != trans_counts.end()){
            int sum = 0;
            for(num_map::iterator it2 = (*it).second.begin(); it2 != (*it).second.end(); it2++){
                sum += (*it2).second;
            }
            return sum;
        }
        return 0;
    }

    inline int count_symbol(int symbol){
        int sum = 0;
        for(type_num_map::iterator it = trans_counts.begin(); it != trans_counts.end(); it++){
            num_map::iterator it2 = (*it).second.find(symbol);
            if(it2 != (*it).second.end()){
                sum += (*it2).second;
            }
        }
        return sum;
    }
    
    inline type_num_map::iterator counts_begin(){
        return trans_counts.begin();
    }

    inline type_num_map::iterator counts_end(){
        return trans_counts.end();
    }

    inline num_map::iterator type_counts_begin(int type){
        num_map& nm = trans_counts[type];
        return nm.begin();
    }

    inline num_map::iterator type_counts_end(int type){
        num_map& nm = trans_counts[type];
        return nm.end();
    }
    
    inline int pos(int symbol){
        return count(1,symbol);
    }

    inline int neg(int symbol){
        return count(0,symbol);
    }
    
    inline num_map::iterator pos_begin(){
        return type_counts_begin(1);
    }

    inline num_map::iterator pos_end(){
        return type_counts_begin(1);
    }
    
    inline num_map::iterator neg_begin(){
        return type_counts_begin(0);
    }

    inline num_map::iterator neg_end(){
        return type_counts_begin(0);
    }

    /*inline int neg(int i){
        num_map::iterator it = trans_counts.find(0);
        if(it != trans_counts.end()){
            num_map::iterator it2 = (*it).find(i);
            if(it2 == (*it).end()) return 0;
            return (*it2).second;
        }
        return 0;
    }*/
    
    alergia_data();

    virtual void read_from(int type, int index, int length, int symbol, string data);
    virtual void print_transition_label(iostream& output, int symbol, apta* apta_context);
    virtual void print_state_label(iostream& output);
    virtual void update(evaluation_data* right);
    virtual void undo(evaluation_data* right);
    
    virtual bool is_low_count_sink();
    virtual bool is_stream_sink(apta_node* node);
    virtual int sink_type(apta_node* node);
    virtual bool sink_consistent(int type);
};


class alergia: public count_driven {

protected:
  REGISTER_DEC_TYPE(alergia);

public:
  static bool alergia_consistency(double right_count, double left_count, double right_total, double left_total);

  virtual bool data_consistent(alergia_data* l, alergia_data* r);
  virtual bool consistent(state_merger *merger, apta_node* left, apta_node* right);
  //virtual void print_dot(iostream&, state_merger *);

  //virtual int sink_type(apta_node* node);
  //virtual bool sink_consistent(apta_node* node, int type);
  virtual int num_sink_types();
};

#endif
