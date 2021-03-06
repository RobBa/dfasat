
#ifndef _INPUTDATA_H_
#define _INPUTDATA_H_

class inputdata;
class tail;
class tail_wrapper;

#include <istream>
#include <sstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>

#include "json.hpp"

// for convenience
using json = nlohmann::json;
using namespace std;

// sequence-index pairs
//typedef pair<int,int> tail;
//typedef vector<tail> tail_list;
//typedef set<tail> tail_set;

#include "apta.h"

class tail_data{
public:
    int sequence;
    int index;

    int type;
    int length;
    int symbol;
    int* attr;
    int* trace_attr;
    string data;
};

class tail{
public:
    tail(tail *ot);
    void initialize(tail* ot);

    tail_data* td;

    tail* future_tail;
    tail* past_tail;
    tail* next_in_list;
    tail* split_from;
    tail* split_to;
    
    tail(int seq, int i, tail* past_tail);
    ~tail();

    void split(tail* t);
    void undo_split();
    tail* next();
    tail* future();
    tail* past();
    tail* splitted();

    inline int get_index(){
        return td->index;
    };
    inline int get_type(){
        return td->type;
    };
    inline int get_length(){
        return td->length;
    };
    inline int get_sequence(){
        return td->sequence;
    };
    inline int get_symbol(){
        return td->symbol;
    }

    const string to_string();
};

/**
 * @brief Wrapper class for the input data. Supports functionalities 
 * such as alphabet functions, file transformations and data added to the APTA.
 * 
 */
class attribute{
public:
    bool discrete;
    bool splittable;
    bool distributionable;
    bool target;

    vector<string> values;
    map<string, int> r_values;

    attribute(string);

    inline float get_value(string val){
        if(discrete){
            if(r_values.find(val) == r_values.end()) {
                r_values[val] = values.size();
                values.push_back(val);
            }
            return (float) r_values[val];
        } else {
            return stof(val);
        }
    };
};

class inputdata{
public:
    // TODO: not all public. Can we circumvent this class to have cleaner code?
    static json all_data;

    static vector<string> alphabet;
    static map<string, int> r_alphabet;

    static vector<string> types;
    static map<string, int> r_types;

    static vector<attribute> trace_attributes;
    static vector<attribute> symbol_attributes;

    static int num_sequences;
    static int alphabet_size;

    static int num_trace_attributes;
    static int num_symbol_attributes;
    static int num_attributes;

    static int node_number;

    void read_json_file(istream &input_stream);
    void read_abbadingo_file(istream &input_stream);
    void read_abbadingo_sequence(istream &input_stream);

    static inline bool is_splittable(int attr){
        return get_attribute(attr)->splittable;
    };
    static inline bool is_distributionable(int attr){
        return get_attribute(attr)->distributionable;
    };
    static inline bool is_discrete(int attr){
        return get_attribute(attr)->discrete;
    };
    static inline bool is_target(int attr){
        return get_attribute(attr)->target;
    };
    static inline attribute* get_attribute(int attr){
        if(attr < inputdata::num_symbol_attributes){
            return &inputdata::symbol_attributes[attr];
        }
        return &inputdata::trace_attributes[attr - inputdata::num_symbol_attributes];
    };
    static inline int get_size(){
        return inputdata::all_data.size();
    };
    static inline int get_type(int seq_nr){
        return inputdata::all_data[seq_nr]["T"];
    };
    static inline int get_length(int seq_nr){
        return inputdata::all_data[seq_nr]["L"];
    };
    static inline int get_symbol(int seq_nr, int index){
        if(index > -1)
            return inputdata::all_data[seq_nr]["S"][index];
        return -1;
    };
    static inline float get_value(int seq_nr, int index, int attr){
        if(index > -1){
            if(attr < inputdata::num_symbol_attributes){
                return inputdata::all_data[seq_nr]["V" + to_string(attr)][index];
            }
            return inputdata::all_data[seq_nr]["VT"][attr - inputdata::num_symbol_attributes];
        }
        return -1;
    };
    static inline float get_trace_value(int seq_nr, int index, int attr){
        if(index > -1)
            return inputdata::all_data[seq_nr]["VT"][attr];
        return -1;
    };
    static inline string get_data(int seq_nr, int index){
        if(index > 0 && inputdata::all_data[seq_nr].find("D") != inputdata::all_data[seq_nr].end())
            return inputdata::all_data[seq_nr]["D"][index];
        return "";
    };
    
    static inline int get_type(tail* t){
        return t->td->type;
        //return inputdata::all_data[t->sequence]["T"];
    };
    static inline int get_length(tail* t){
        return t->td->length;
        //return inputdata::all_data[t->sequence]["L"];
    };
    static inline int get_symbol(tail* t){
        return t->td->symbol;
        //if(t->index > -1)
        //    return inputdata::all_data[t->sequence]["S"][t->index];
        //return -1;
    };
    static inline int get_index(tail* t){
        return t->td->index;
    };
    static inline int get_value(tail* t, int a){
        if(a < inputdata::num_symbol_attributes)
            return t->td->attr[a];
        return t->td->trace_attr[a - inputdata::num_symbol_attributes];
        //if(t->index > -1)
        //    return inputdata::all_data[t->sequence]["V" + to_string(a)][t->index];
        //return -1;
    };
    static inline string get_data(tail* t){
        return t->td->data;
        //if(t->index > -1 && inputdata::all_data[t->sequence].find("D") != inputdata::all_data[t->sequence].end())
        //    return inputdata::all_data[t->sequence]["D"][t->index];
        //return "";
    };
	
    void add_data_to_apta(apta* the_apta);
    void add_sequence_to_apta(apta* the_apta, int seq_nr);

    const string to_json_str() const;
	//const string to_abbadingo_str() const;

    // to init counters etc
    inputdata();

    void abbadingo_init(istream &input_stream);
};

#endif /* _INPUTDATA_H_*/
