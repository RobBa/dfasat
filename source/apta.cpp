#include "apta.h"
#include "state_merger.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <unordered_map>
#include <string>
#include <sstream>
#include <iterator>
#include <cassert>

#include "evaluators.h"
#include "parameters.h"
#include "evaluation_factory.h"

using namespace std;

bool apta_guard::bounds_satisfy(tail* t){
    for(bound_map::iterator it = min_attribute_values.begin(); it != min_attribute_values.end(); ++it){
        if(inputdata::get_value(t, (*it).first) < (*it).second) return false;
    }
    for(bound_map::iterator it = max_attribute_values.begin(); it != max_attribute_values.end(); ++it){
        if(inputdata::get_value(t, (*it).first) >= (*it).second) return false;
    }
    return true;
};


/* constructors and destructors */
apta::apta(){
    root = new apta_node(this);
    root->red = true;
    max_depth = 0;
    merge_count = 0;
}

apta::~apta(){
    delete root;
}

/* for batch mode */
// i want tis to map back to the sample by sample read functions from streaming
void apta::read_file(istream &input_stream){
    int num_words;
    int num_alph = 0;
    map<string, int> seen;
    //apta::node_number = 1;
    input_stream >> num_words >> alphabet_size;

    for(int line = 0; line < num_words; line++){
        int type;
        int length;
        string id;
        
        // I need to refactor this...
        if(EXCEPTION4OVERLAP) {
           input_stream >> id >> type >> length;
        }
        else {
           input_stream >> type >> length;
        }

        int depth = 0;
        apta_node* node = root;
        for(int index = 0; index < length; index++){
            depth++;
            string tuple;
            input_stream >> tuple;

            std::stringstream lineStream;
            lineStream.str(tuple);

            string symbol;
            std::getline(lineStream,symbol,'/');
            string data;
            std::getline(lineStream,data);

            if(seen.find(symbol) == seen.end()){
                alphabet[num_alph] = symbol;
                seen[symbol] = num_alph;
                num_alph++;
            }
            int c = seen[symbol];
            if(node->child(c) == 0){
                apta_node* next_node = new apta_node(this);
                node->set_child(c,next_node);
                next_node->source = node;
                next_node->label  = c;
                //next_node->number = apta::node_number++;
                next_node->depth = depth;
            }
            node->size = node->size + 1;
            node->data->read_from(type, index, length, c, data);
            if(EXCEPTION4OVERLAP)
                reinterpret_cast<overlap4logs_data*>(node->data)->store_id(id);
            node = node->child(c);
            node->data->read_to(type, index, length, c, data);
        }
        if(depth > max_depth) max_depth = depth;
        node->type = type;
        node->size = node->size + 1;
    }
};

bool is_sink(apta_node* node){
    //cerr << node->size << endl;
    return node->data->sink_type(node) != -1;
}

int apta::sink_type(apta_node* node) {
    return node->data->sink_type(node);
}

void apta::print_dot(iostream& output){
    output << "digraph DFA {\n";
    output << "\t" << root->find()->number << " [label=\"root\" shape=box];\n";
    output << "\t\tI -> " << root->find()->number << ";\n";
    int ncounter = 0;
    //for(merged_APTA_iterator_func Ait = merged_APTA_iterator_func(root, is_sink); *Ait != 0; ++Ait){
    //for(APTA_iterator Ait = APTA_iterator(root); *Ait != 0; ++Ait){
    for(merged_APTA_iterator Ait = merged_APTA_iterator(root); *Ait != 0; ++Ait){
        apta_node* n = *Ait;
        n->number = ncounter++;
    }
    for(merged_APTA_iterator_func Ait = merged_APTA_iterator_func(root, is_sink); *Ait != 0; ++Ait){
    //for(APTA_iterator Ait = APTA_iterator(root); *Ait != 0; ++Ait){
    //for(merged_APTA_iterator Ait = merged_APTA_iterator(root); *Ait != 0; ++Ait){
        apta_node* n = *Ait;
        output << "\t" << n->number << " [ label=\"";
        output << n->number << ":#" << n->size << "\n";
        //output << n << "\n";
        //output << n->representative << "\n";
        /*if(inputdata::num_attributes > 0){
            for(tail_iterator it = tail_iterator(n); *it != 0; ++it){
                tail* t = *it;
                if(t->past_tail != 0) output << " " << inputdata::get_value(t->past_tail, 0);
            }
        }*/
        output << "\n";
        n->data->print_state_label(output, this);
        output << "\" ";
        n->data->print_state_style(output, this);
        if(n->red == false) output << " style=dotted";
        output << ", penwidth=" << log(n->size);
        output << "];\n";

        // transition labels for item
        map<apta_node*, set<int>> childlabels;
        map<int, set<int>> sinklabels;

        for(guard_map::iterator it = n->guards.begin(); it != n->guards.end(); ++it){
            /*if((*it).second->target == 0) continue;
            apta_node* child = (*it).second->target->find();
            if(sink_type(child) != -1){
                if(sinklabels.find(sink_type(child)) == sinklabels.end())
                    sinklabels[sink_type(child)] = set<int>();
                sinklabels[sink_type(child)].insert( it->first );
            } else {
                if(childlabels.find(child) == childlabels.end())
                    childlabels[child] = set<int>();
                childlabels[child].insert( it->first );
            }
        }
        for(map<apta_node*, set<int>>::iterator it2 = childlabels.begin(); it2 != childlabels.end(); ++it2){
            apta_node* child = (*it2).first;
            set<int> labels  = (*it2).second;
            */
            if((*it).second->target == 0) continue;
            int symbol = (*it).first;
            apta_guard* g = (*it).second;
            apta_node* child = (*it).second->target->find();
            if(sink_type(child) != -1) continue;
            output << "\t\t" << n->number << " -> " << child->number << " [label=\"";
            /*for(set<int>::iterator its = labels.begin(); its != labels.end(); its++){
                output << " " << inputdata::alphabet[*its];
            }*/
            
            output << alph_str((*it).first) << endl;
            //output << (*it).first << endl;

            n->data->print_transition_label(output, (*it).first, this);
            
            
            for(bound_map::iterator it2 = g->min_attribute_values.begin(); it2 != g->min_attribute_values.end(); ++it2){
                output << " " << (*it2).first << " >= " << (*it2).second;
            }
            for(bound_map::iterator it2 = g->max_attribute_values.begin(); it2 != g->max_attribute_values.end(); ++it2){
                output << " " << (*it2).first << " < " << (*it2).second;
            }
            
             
            /*for(set<int>::iterator it3 = labels.begin(); it3 != labels.end(); ++it3){
                output << alph_str(*it3) << ":";
                n->data->print_transition_label(output, *it3, this);
                if(std::next(it3) != labels.end()) output << ",";
            }*/

            output << "\" ";
            if(alph_str((*it).first) == "13" || alph_str((*it).first) == "31" || alph_str((*it).first) == "9" || alph_str((*it).first) == "10" || alph_str((*it).first) == "32") output << ", color=red";
            //n->data->print_transition_style(output, labels, this);
            output << " ];\n";
        }
        /*
        for(map<int, set<int>>::iterator it2 = sinklabels.begin(); it2 != sinklabels.end(); ++it2){
            int stype = (*it2).first;
            set<int> labels  = (*it2).second;
            
            output << "\tS" << n->number << "t" << stype << " [ label=\"";
            for(set<int>::iterator it3 = labels.begin(); it3 != labels.end(); ++it3){
                output << n->get_child(*it3)->size << " ";
                n->get_child(*it3)->data->print_state_label(output, this);
            }
            output << "\"\n ";
            n->get_child(*(labels.begin()))->data->print_state_style(output, this);
            if(n->red == false) output << " shape=box";
            output << " ];\n";

            output << "\t\t" << n->number << " -> S" << n->number << "t" << stype << " [ label=\"";
            
            for(set<int>::iterator it3 = labels.begin(); it3 != labels.end(); ++it3){
                output << alph_str(*it3) << ":";
                n->data->print_transition_label(output, *it3, this);
                output << "\n";
            }

            output << "\" ";
            n->data->print_transition_style(output, labels, this);
            output << " ];\n";
        }
        */
    }
    output << "}\n";
};


void apta::print_json(iostream& output){
    output << "{\n";
    output << "\t\"nodes\" : [\n";

    int count = 0;

    for(merged_APTA_iterator Ait = merged_APTA_iterator(root); *Ait != 0; ++Ait){
        apta_node* n = *Ait;
        n->number = count++;
    }
    // output states
    for(merged_APTA_iterator_func Ait = merged_APTA_iterator_func(root, is_sink); *Ait != 0; ++Ait){

        apta_node* n = *Ait;
        if(n->number > 0)
            output << ",\n";

        output << "\t\t{\n";
        output << "\t\t\t\"id\" : " << n->number << ",\n";
        output << "\t\t\t\"label\" : \"";
        n->data->print_state_label(output, this);
        output  << "\",\n"; 
        output << "\t\t\t\"size\" : " << n->size << ",\n";
        output<< "\t\t\t\"level\" : " << n->depth << ",\n";
        output << "\t\t\t\"style\" : \"";
        n->data->print_state_style(output, this);
        output  << "\",\n";
        output << "\t\t\t\"isred\" :  ";
        
        if(n->red == false) 
            output << "\"false\"";
        else
            output << "\"true\"";

        output << "\t\t\t\"traces\" : { ";
        for(tail_iterator it = tail_iterator(n); *it != 0; ++it) {
            tail* t = *it;
            output << t->sequence << ",";
        }
        output << "}\" :  ";

        output << "\n\t\t}";
 
       count++;
    }

    output << "\n\t],\n";
    output << "\t\"edges\" : [\n";

    // output transitions
    count = 0;

    for(merged_APTA_iterator_func Ait = merged_APTA_iterator_func(root, is_sink); *Ait != 0; ++Ait){

        apta_node* n = *Ait;

        // transition labels for item
        map<apta_node*, set<int>> childlabels;
        map<int, set<int>> sinklabels;

        for(guard_map::iterator it = n->guards.begin(); it != n->guards.end(); ++it){
            if((*it).second->target == 0) continue;
            apta_node* child = (*it).second->target->find();
            if(sink_type(child) != -1){
                continue;
                if(sinklabels.find(sink_type(child)) == sinklabels.end())
                    sinklabels[sink_type(child)] = set<int>();
                sinklabels[sink_type(child)].insert( it->first );
            } else {
                if(childlabels.find(child) == childlabels.end())
                    childlabels[child] = set<int>();
                childlabels[child].insert( it->first );
            }

            if(count > 0)
                output << ",\n";
            else
                output << "\n";

            output << "\t\t{\n";
            output << "\t\t\t\"id\" : \"" << n->number << "_" << child->number << "\",\n";
            output << "\t\t\t\"source\" : \"" << n->number << "\",\n";
            output << "\t\t\t\"target\" : \"" << child->number << "\",\n";

            output << "\t\t\t\"name\": \"" << (*it).first << "\",\n";
            output << "\t\t\t\"appearances\": \"";
            n->data->print_transition_label(output, (*it).first, this);
            output << "\"}\n";
            count++;
        }

//        for(map<apta_node*, set<int>>::iterator it2 = childlabels.begin(); it2 != childlabels.end(); ++it2){
//            apta_node* child = (*it2).first;
//            set<int> labels  = (*it2).second;
//
//            if(count > 0)
//                output << ",\n";
//            else
//                output << "\n";
//
//            output << "\t\t{\n";
//            output << "\t\t\t\"id\" : \"" << n->number << "_" << child->number << "\",\n";
//            output << "\t\t\t\"source\" : \"" << n->number << "\",\n";
//            output << "\t\t\t\"target\" : \"" << child->number << "\",\n";
//
//
//            output << "\t\t\t\"label\" : \"";
//	    // transition label information
//            for(set<int>::iterator it3 = labels.begin(); it3 != labels.end(); ++it3){
//                output << alph_str(*it3);
//                //n->data->print_transition_label(output, *it3, this);
//                if(std::next(it3) != labels.end()) output << ",";
//            }
//            output << "\",\n";
//
//            output << "\t\t\t\"labelinfo\" : [\n";
//	    // transition label information
//            for(set<int>::iterator it3 = labels.begin(); it3 != labels.end(); ++it3){
//                if (it3!=labels.begin()){
//                    output << ",";
//                    output << "\n";
//                }
//                output << "\t\t\t\t { " << "\"symbol\" : \"" << alph_str(*it3) << "\"}";
//                n->data->print_transition_properties(output, *it3, this);
//            }
//	        output << "\t\t\t]\n";
//            output << "\n";
//
//            // n->data->print_transition_style(output, labels, this);
//            output << "\t\t}";
//
//           count++;
//        }

    }
    // n->data->print_transition_style(output, labels, this);
    output << "\n\t]\n}\n";
};


string apta::alph_str(int i){
    return alp[i];
}

apta_guard::apta_guard(){
    target = 0;
    undo = 0;
}

apta_guard::apta_guard(apta_guard* g){
    target = 0;
    undo = 0;
    
    min_attribute_values = bound_map(g->min_attribute_values);
    max_attribute_values = bound_map(g->max_attribute_values);
}

void apta_node::add_tail(tail* t){
    t->next_in_list = tails_head;
    tails_head = t;
    //data->add_tail(t);
};

apta_node::apta_node(apta *context) {
    this->context = context;
    source = 0;
    representative = 0;
    next_merged_node = 0;
    representative_of = 0;
    
    //children = child_map();
    //det_undo = child_map();
    tails_head = 0;

    label = 0;
    number = -1;
    satnumber = 0;
    colour = 0;
    size = 0;
    final = 0;
    depth = 0;
    type = -1;

    age = 0;
    
    red = false;

    try {
       data = (DerivedDataRegister<evaluation_data>::getMap())->at(eval_string)();
       data->node = this;
    } catch(const std::out_of_range& oor ) {
       std::cerr << "No data type found..." << std::endl;
    }

}

apta_node::apta_node(){
    source = 0;
    representative = 0;
    next_merged_node = 0;
    representative_of = 0;
    
    //children = child_map();
    //det_undo = child_map();
    tails_head = 0;
    
    label = 0;
    number = -1;
    satnumber = 0;
    colour = 0;
    size = 1;
    final = 0;
    depth = 0;
    type = -1;

    red = false;

    try {
       data = (DerivedDataRegister<evaluation_data>::getMap())->at(eval_string)();
       data->node = this;
    } catch(const std::out_of_range& oor ) {
       std::cerr << "No data type found..." << std::endl;
    }
}

apta_node* apta_node::child(tail* t){
        int symbol = inputdata::get_symbol(t);
        for(guard_map::iterator it = guards.lower_bound(symbol); it != guards.upper_bound(symbol); ++it){
            if((*it).first != symbol) break;
            apta_guard* g = (*it).second;
            bool outside_range = false;
            for(bound_map::iterator it2 = g->min_attribute_values.begin(); it2 != g->min_attribute_values.end(); ++it2){
                if(inputdata::get_value(t,(*it2).first) < (*it2).second){
                    outside_range = true;
                    break;
                }
            }
            if(outside_range) continue;
            for(bound_map::iterator it2 = g->max_attribute_values.begin(); it2 != g->max_attribute_values.end(); ++it2){
                if(inputdata::get_value(t,(*it2).first) > (*it2).second){
                    outside_range = true;
                    break;
                }
            }
            if(outside_range) continue;
            return g->target;
        }
        return 0;
};

apta_guard* apta_node::guard(int symbol, apta_guard* g){
        for(guard_map::iterator it = guards.lower_bound(symbol); it != guards.upper_bound(symbol); ++it){
            if((*it).first != symbol) break;
            apta_guard* g2 = (*it).second;
            bool outside_range = false;
            for(bound_map::iterator it2 = g->min_attribute_values.begin(); it2 != g->min_attribute_values.end(); ++it2){
                bound_map::iterator it3 = g2->min_attribute_values.find((*it2).first);
                if(it3 == g2->min_attribute_values.end() || (*it3).second != (*it2).second){
                    outside_range = true;
                    break;
                }
            }
            if(outside_range) continue;
            for(bound_map::iterator it2 = g->max_attribute_values.begin(); it2 != g->max_attribute_values.end(); ++it2){
                bound_map::iterator it3 = g2->max_attribute_values.find((*it2).first);
                if(it3 == g2->max_attribute_values.end() || (*it3).second != (*it2).second){
                    outside_range = true;
                    break;
                }
            }
            if(outside_range) continue;
            return g2;
        }
        return 0;
};

apta_guard* apta_node::guard(tail* t){
    int symbol = inputdata::get_symbol(t);
    guard_map::iterator it = guards.lower_bound(symbol);
    guard_map::iterator it_end = guards.upper_bound(symbol);
    for(;it != it_end; ++it){
        if((*it).second->bounds_satisfy(t)){
            return (*it).second;
        }
    }
    return 0;
};

/* iterators for the APTA and merged APTA */
APTA_iterator::APTA_iterator(apta_node* start){
    base = start;
    current = start;
}

apta_node* APTA_iterator::next_forward() {
    guard_map::iterator it;
    for(it = current->guards.begin();it != current->guards.end(); ++it){
        apta_node* target = (*it).second->target;
        if(target != 0 && target->source == current){
            return target;
        }
    }
    return 0;
}

apta_node* APTA_iterator::next_backward() {
    guard_map::iterator it;
    apta_node* source = current;
    while(source != base){
        current = source;
        source = source->source;
        it = source->guards.begin();//source->guards.find(current->label);
        while((*it).second->target != current) ++it;
        ++it;
        for(; it != source->guards.end(); ++it){
            apta_node* target = (*it).second->target;
            if(target != 0 && target->source == source){
                return target;
            }
        }
    }
    return 0;
}

void APTA_iterator::increment() {
    apta_node* next = next_forward();
    if(next != 0){ current = next; return; }
    next = next_backward();
    if(next != 0){ current = next; return; }
    current = 0;
}

merged_APTA_iterator::merged_APTA_iterator(apta_node* start){
    base = start;
    current = start;
}

apta_node* merged_APTA_iterator::next_forward() {
    guard_map::iterator it;
    for(it = current->guards.begin();it != current->guards.end(); ++it){
        apta_node* target = (*it).second->target;
        if(target != 0 && target->representative == 0){
            return target;
        }
    }
    return 0;
}

apta_node* merged_APTA_iterator::next_backward() {
    guard_map::iterator it;
    apta_node* source = current;
    while(source != base){
        current = source->find();
        source = source->source->find();
        it = source->guards.find(current->label);
        it = source->guards.begin();
        while((*it).second->target != current){
            ++it;
            if(it == source->guards.end()) return 0;
        }
        ++it;
        for(; it != source->guards.end(); ++it){
            apta_node* target = (*it).second->target;
            if(target != 0 && target->representative == 0){
                return target;
            }
        }
    }
    return 0;
}

void merged_APTA_iterator::increment() {
    apta_node* next = next_forward();
    if(next != 0){ current = next; return; }
    next = next_backward();
    if(next != 0){ current = next; return; }
    current = 0;
}

void blue_state_iterator::increment() {
    apta_node* next = next_backward();
    while(next != 0){
        while(next != 0){
            current = next;
            if(!current->red) return;
            next = next_forward();
        }
        next = next_backward();
    }
    current = 0;
}

blue_state_iterator::blue_state_iterator(apta_node* start) :
    merged_APTA_iterator(start) {
    apta_node* next = start;
    while(next != 0){
        while(next != 0){
            current = next;
            if(!current->red) return;
            next = next_forward();
        }
        next = next_backward();
    }
    current = 0;
}

void red_state_iterator::increment() {
    apta_node* next = next_forward();
    if(next != 0){
        current = next;
        if(current->red) return;
    }
    next = next_backward();
    while(next != 0){
        current = next;
        if(current->red) return;
        next = next_backward();
    }
    current = 0;
}

red_state_iterator::red_state_iterator(apta_node* start) :
    merged_APTA_iterator(start) {
}

void merged_APTA_iterator_func::increment() {
    apta_node* next = next_forward();
    if(next != 0){
        current = next;
        if(!check_function(current)) return;
    }
    next = next_backward();
    while(next != 0){
        current = next;
        if(!check_function(current)) return;
        next = next_backward();
    }
    current = 0;
}

merged_APTA_iterator_func::merged_APTA_iterator_func(apta_node* start, bool(*node_check)(apta_node*)) : merged_APTA_iterator(start){
    check_function = node_check;
}

tail_iterator::tail_iterator(apta_node* start){
    base = start;
    current = start;
    current_tail = current->tails_head;
    while(current_tail == 0){
        if(current == 0) return;
        next_node();
        if(current == 0) return;
        current_tail = current->tails_head;
    }
    if(current_tail != 0 && current_tail->split_to != 0) increment();
}

void tail_iterator::next_node(){
    if(current->representative_of != 0) current = current->representative_of;
    else if (current->next_merged_node != 0) current = current->next_merged_node;
    else {
        while(current != 0 && current->next_merged_node == 0) current = current->representative;
        if(current != 0) current = current->next_merged_node;
    }
}

void tail_iterator::increment() {
    current_tail = current_tail->next_in_list;
    if(current_tail != 0 && current_tail->split_to != 0) increment();
    
    while(current_tail == 0){
        if(current == 0) return;
        next_node();
        if(current == 0) return;
        current_tail = current->tails_head;
    }
    if(current_tail != 0 && current_tail->split_to != 0) increment();
}

/*apta_node::add_target(int symbol){
    if(node->child(symbol) == 0){
        apta_node* next_node = new apta_node();
        node->children[symbol] = next_node;
        next_node->source = this;
        next_node->label  = symbol;
        next_node->number = apta::node_number++;
        next_node->depth = depth+1;
    }
    size = size + 1;
}*/

std::set<void*> freed;

apta_node::~apta_node(){
    for(guard_map::iterator it = guards.begin();it != guards.end(); ++it){
        if (freed.find((*it).second->target) != freed.end()) {
            freed.insert((*it).second->target);
            //delete (*it).second->target;
            delete (*it).second;
        }
    }
    tail* t = tails_head;
    tail* n = 0;
    while(t != 0){
        n = t->next_in_list;
        delete t;
        t = n;
    }
    delete data;
}
