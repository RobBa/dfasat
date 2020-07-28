#include <math.h>
#include <queue>
#include "refinement.h"
#include "parameters.h"

using namespace std;

merge_refinement::merge_refinement(state_merger* m, double s, apta_node* l, apta_node* r){
    red = m->get_tail_from_state(l);
    blue = m->get_tail_from_state(r);
    tempnode = l;
    tempblue = r;
    score = s;
}

split_refinement::split_refinement(state_merger* m, double s, apta_node* r, tail* t, int a){
    split_point = t;
    red = m->get_tail_from_state(r);
    score = s;
    attribute = a;
    tempnode = r;
}

extend_refinement::extend_refinement(state_merger* m, apta_node* r){
    //cerr << r << endl;
    red = m->get_tail_from_state(r);
    score = LOWER_BOUND;
    tempnode = r;
}

inline void refinement::print() const{
    cerr << "score " << score << endl;
};
	
inline void refinement::print_short() const{
    cerr << score;
};

inline void refinement::doref(state_merger* m){
};
	
inline void refinement::undo(state_merger* m){
};

inline bool refinement::testref(state_merger* m){
    return true;
};

inline void merge_refinement::print() const{
    //cerr << "merge( " << score << " " << left->number << " " << right->number << " )" << endl;
};
	
inline void merge_refinement::print_short() const{
    cerr << "m" << score;
};

inline void merge_refinement::doref(state_merger* m){
    apta_node* left = m->get_state_from_tail(red);
    apta_node* right = m->get_state_from_tail(blue);
    //left = tempnode;
    //right = tempblue;
    m->perform_merge(left, right);
};
	
inline void merge_refinement::undo(state_merger* m){
    //cerr << "undo merge" << endl;
    apta_node* left = m->get_state_from_tail(red);
    apta_node* right = m->get_state_from_tail(blue);
    //left = tempnode;
    //right = tempblue;
    m->undo_perform_merge(left, right);
};

inline bool merge_refinement::testref(state_merger* m){
    apta_node* left = m->get_state_from_tail(red);
    apta_node* right = m->get_state_from_tail(blue);
    //left = tempnode;
    //right = tempblue;
    if(left->red == false || right->red == true || right->source->find()->red == false) return false;
    if(left->representative != 0 || right->representative != 0) return false;
    refinement* ref = m->test_merge(left, right);
    if(ref != 0){
        score = ref->score;
        return true;
    }
    return false;
};

inline void split_refinement::print() const{
    //cerr << "split( " << score << " q:" << right->number << " s:" << inputdata::get_symbol(split_point) << " a:" <<attribute << " v:" << inputdata::get_value(split_point, attribute) << " )";
};
	
inline void split_refinement::print_short() const{
    cerr << "s" << score;
};

inline void split_refinement::doref(state_merger* m){
    apta_node* right = m->get_state_from_tail(red);
    //right = tempnode;
    m->perform_split(right, split_point, attribute);
};
	
inline void split_refinement::undo(state_merger* m){
    //cerr << "undo split" << endl;
    apta_node* right = m->get_state_from_tail(red);
    //right = tempnode;
    m->undo_perform_split(right, split_point, attribute);
};

inline bool split_refinement::testref(state_merger* m){
    //return false;
    apta_node* right = m->get_state_from_tail(red);
    //right = tempnode;
    if(right->red == true || right->source->find()->red == false) return false;
    if(right->representative != 0) return false;
    refinement* ref = m->test_split(right, split_point, attribute);
    if(ref != 0){
        score = ref->score;
        return true;
    }
    return false;
};

inline void extend_refinement::print() const{
    //cerr << "extend( " << right->number << " )" << endl;
};
	
inline void extend_refinement::print_short() const{
    cerr << "x";
};

inline void extend_refinement::doref(state_merger* m){
    apta_node* right = m->get_state_from_tail(red);
    //right = tempnode;
    if(right->red == true){ cerr << "error3" << endl; assert(false); };
        if(right->source->find()->red == false){ cerr << "error" << endl; };
    if(right->representative != 0){ cerr << "error2" << endl; };
    m->extend(right);
};
	
inline void extend_refinement::undo(state_merger* m){
    //cerr << "undo extend" << endl;
    apta_node* right = m->get_state_from_tail(red);
    //right = tempnode;
    m->undo_extend(right);
};

inline bool extend_refinement::testref(state_merger* m){
    apta_node* right = m->get_state_from_tail(red);
    //right = tempnode;
    if(right->red == true || right->source->find()->red == false) return false;
    if(right->representative != 0) return false;
    return true;
};
