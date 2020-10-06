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

void merge_refinement::initialize(state_merger* m, double s, apta_node* l, apta_node* r){
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

void split_refinement::initialize(state_merger* m, double s, apta_node* r, tail* t, int a){
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

void extend_refinement::initialize(state_merger* m, apta_node* r){
    //cerr << r << endl;
    red = m->get_tail_from_state(r);
    score = LOWER_BOUND;
    tempnode = r;
}

inline void refinement::print() const{
    cerr << "score " << score << endl;
};

inline void refinement::print_json(iostream& output) const{
    output << "\t\t[\n";
    output << "score " << score << endl;
    output << "\t\t]\n";
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

inline void refinement::erase(){
};

inline void merge_refinement::print() const{
    //cerr << "merge( " << score << " " << left->number << " " << right->number << " )" << endl;
};
	
inline void merge_refinement::print_short() const{
    cerr << "m" << score;
};

inline void merge_refinement::print_json(iostream& output) const{
    output << "\t\t{\n";
    output << "\t\t\t\"type\" : \"merge\", " << endl;
    if(red->past() != 0){
        output << "\t\t\t\"red\" : " << red->past()->to_string() << "," << endl;
    } else {
        output << "\t\t\t\"red\" : " << 0 << "," << endl;
    }
    output << "\t\t\t\"blue\" : " << blue->past()->to_string() << "," << endl;
    output << "\t\t\t\"score\" : " << score << endl;
    output << "\t\t}\n";
};

inline void merge_refinement::doref(state_merger* m){
    apta_node* left = m->get_state_from_tail(red);
    apta_node* right = m->get_state_from_tail(blue);
    //left = tempnode;
    //right = tempblue;
    //cerr << "merge " << left << " " << right << endl;
    if(left->red == false){
        m->extend(left);
        right->red = true;
    }
    m->perform_merge(left, right);
};
	
inline void merge_refinement::undo(state_merger* m){
    //cerr << "undo merge" << endl;
    apta_node* left = m->get_state_from_tail(red);
    apta_node* temp = m->get_state_from_tail(blue->past_tail);
    apta_node* right =  temp->guard(blue->past_tail)->target;
    right = right->find_until(left);
    //left = tempnode;
    //right = tempblue;
    //cerr << "undo merge " << left << " " << right << endl;
    m->undo_perform_merge(left, right);
    if(right->red == true){
        right->red = false;
        m->undo_extend(left);
    }
};

inline bool merge_refinement::testref(state_merger* m){
    apta_node* left = m->get_state_from_tail(red);
    apta_node* right = m->get_state_from_tail(blue);
    //left = tempnode;
    //right = tempblue;
    if((left->red == false && left->source->find()->red == false) || right->red == true || right->source->find()->red == false) return false;
    if(left->representative != 0 || right->representative != 0) return false;
    refinement* ref = m->test_merge(left, right);
    if(ref != 0){
        score = ref->score;
        return true;
    }
    return false;
};

inline void merge_refinement::erase(){
    mem_store::delete_merge_refinement(this);
};

inline void split_refinement::print() const{
    //cerr << "split( " << score << " q:" << right->number << " s:" << inputdata::get_symbol(split_point) << " a:" <<attribute << " v:" << inputdata::get_value(split_point, attribute) << " )";
};
	
inline void split_refinement::print_short() const{
    cerr << "s" << score;
};

inline void split_refinement::print_json(iostream& output) const{
    output << "\t\t{\n";
    output << "\t\t\t\"type\" : \"split\", " << endl;
    output << "\t\t\t\"red\" : " << red->past()->to_string() << "," << endl;
    output << "\t\t\t\"point\" : " << split_point->past()->to_string() << "," << endl;
    output << "\t\t\t\"attribute\" : " << attribute << "," << endl;
    output << "\t\t\t\"score\" : " << score << endl;
    output << "\t\t}\n";
};

inline void split_refinement::doref(state_merger* m){
    apta_node* right = m->get_state_from_tail(red);
    //right = tempnode;
    //cerr << "split " << right << endl;
    m->perform_split(right, split_point, attribute);
};
	
inline void split_refinement::undo(state_merger* m){
    //cerr << "undo split" << endl;
    apta_node* right = m->get_state_from_tail(red);
    //right = tempnode;
    //cerr << "undo split " << right << endl;
    m->undo_perform_split(right, split_point, attribute);
};

inline bool split_refinement::testref(state_merger* m){
    //return false;
    apta_node* right = m->get_state_from_tail(red);
    //right = tempnode;
    if(right->red == false) return false;
    if(right->guard(split_point) == 0) return false;
    if(right->guard(split_point)->target == 0) return false;
    if(right->guard(split_point)->target->representative != 0) return false;
    if(right->guard(split_point)->target->red == true) return false;
    refinement* ref = m->test_split(right, split_point, attribute);
    if(ref != 0){
        score = ref->score;
        return true;
    }
    return false;
};

inline void split_refinement::erase(){
    mem_store::delete_split_refinement(this);
};

inline void extend_refinement::print() const{
    //cerr << "extend( " << right->number << " )" << endl;
};
	
inline void extend_refinement::print_short() const{
    cerr << "x";
};

inline void extend_refinement::print_json(iostream& output) const{
    output << "\t\t{\n";
    output << "\t\t\t\"type\" : \"extend\", " << endl;
    output << "\t\t\t\"red\" : " << red->past()->to_string() << "," << endl;
    output << "\t\t\t\"score\" : " << score << endl;
    output << "\t\t}\n";
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

inline void extend_refinement::erase(){
    mem_store::delete_extend_refinement(this);
};

void refinement::print_refinement_list_json(iostream& output, refinement_list* list){
    output << "[\n";

    refinement_list::iterator it = list->begin();
    while(it != list->end()){
        (*it)->print_json(output);
        it++;
        if(it != list->end())
            output << ",\n";
    }

    output << "]\n";
};


