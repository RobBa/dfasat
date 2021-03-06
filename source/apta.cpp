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

#include "utility/loguru.hpp"

using namespace std;

int apta_node::count_tails(){
    int count = 0;
    int final_count = 0;
    for(tail_iterator it = tail_iterator(this); *it != 0; ++it) {
        count++;
        //cerr << *it << endl;
        if((*it)->future() == 0) final_count++;
    }
    //cerr << this << " " << count << " " << final_count << " " << size << endl;
    return count;
};

tail* apta_node::get_tail_from_state(){
    tail_iterator it = tail_iterator(this);
    tail* t = *it;
    while(t->split_from != 0){
        t = t->split_from;
    }
    return t;
};


bool apta_guard::bounds_satisfy(tail* t){
    for(bound_map::iterator it = min_attribute_values.begin(); it != min_attribute_values.end(); ++it){
        if(inputdata::get_value(t, it->first) < it->second) return false;
    }
    for(bound_map::iterator it = max_attribute_values.begin(); it != max_attribute_values.end(); ++it){
        if(inputdata::get_value(t, it->first) >= it->second) return false;
    }
    return true;
};


/* constructors and destructors */
apta::apta(){

    LOG_S(INFO) << "Creating APTA data structure";
    root = new apta_node(this);
    root->red = true;
    max_depth = 0;
    merge_count = 0;
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
    //for(merged_APTA_iterator_func Ait = merged_APTA_iterator_func(root, is_sink); *Ait != 0; ++Ait){
    //for(APTA_iterator Ait = APTA_iterator(root); *Ait != 0; ++Ait){
    for(merged_APTA_iterator Ait = merged_APTA_iterator(root); *Ait != 0; ++Ait) {
        apta_node *n = *Ait;

        if (n->data->print_state_true(this) == false) {
            continue;
        }

        if (!PRINT_WHITE && n->red == false) {
            if (n->source != 0) {
                if (n->source->find()->red == false)
                    continue;
                if (!PRINT_BLUE)
                    continue;
            }
        }

        output << "\t" << n->number << " [ label=\"";
        output << n->number << ":#" << n->size << "\n";
        n->data->print_state_label(output, this);
        output << "\" ";
        n->data->print_state_style(output, this);
        if (n->red == false) output << " style=dotted";
        output << ", penwidth=" << log(1 + n->size);
        output << "];\n";

        for(guard_map::iterator it = n->guards.begin(); it != n->guards.end(); ++it){
            if(it->second->target == 0) continue;
            int symbol = it->first;
            apta_guard* g = it->second;
            apta_node* child = it->second->target->find();

            if (!PRINT_WHITE && child->red == false) {
                if (n->red == false)
                    continue;
                if (!PRINT_BLUE)
                    continue;
            }

            output << "\t\t" << n->number << " -> " << child->number << " [label=\"";
            output << inputdata::alphabet[it->first] << endl;

            n->data->print_transition_label(output, it->first, this);
            
            
            for(bound_map::iterator it2 = g->min_attribute_values.begin(); it2 != g->min_attribute_values.end(); ++it2){
                output << " " << it2->first << " >= " << it2->second;
            }
            for(bound_map::iterator it2 = g->max_attribute_values.begin(); it2 != g->max_attribute_values.end(); ++it2){
                output << " " << it2->first << " < " << it2->second;
            }

            output << "\" ";
            output << " ];\n";
        }
    }
    output << "}\n";
};

void apta_node::print_json(iostream& output){
    output << "\t\t{\n";
    output << "\t\t\t\"id\" : " << number << ",\n";
    //output << "\t\t\t\"id\" : " << get_tail_from_state()->to_string() << ",\n";
    output << "\t\t\t\"label\" : \"";
    data->print_state_label(output, context);
    output  << "\",\n";
    output << "\t\t\t\"size\" : " << size << ",\n";
    output<< "\t\t\t\"level\" : " << compute_depth() << ",\n";
    output << "\t\t\t\"style\" : \"";
    data->print_state_style(output, context);
    output  << "\",\n";
    output << "\t\t\t\"isred\" :  " << red << ",\n";
    if(source != 0 && red == false) output << "\t\t\t\"isblue\" :  " << source->find()->red << ",\n";
    else output << "\t\t\t\"isblue\" :  " << false << ",\n";

    output << "\t\t\t\"traces\" : [ ";
    for(tail_iterator it = tail_iterator(this); *it != 0; ++it) {
        tail* t = *it;

        tail_iterator secIter = it;

        apta_node* next = secIter.next_forward();
        if (next != 0) {
            output << t->get_sequence() << ",";
        }else{
            output << t->get_sequence();
        }

    }
    output << "]";
    output << "\n\t\t}";
};

void apta_node::print_json_transitions(iostream& output){
    bool first = true;
    for(guard_map::iterator it = guards.begin(); it != guards.end(); ++it){
        if(it->second->target == 0) continue;

        if(!first) output << ",\n";
        else first = false;

        int symbol = it->first;
        apta_guard* g = it->second;
        apta_node* child = it->second->target->find();

        output << "\t\t{\n";
        //output << "\t\t\t\"id\" : " << get_tail_from_state()->to_string() << "_" << child->get_tail_from_state()->to_string() << ",\n";
        //output << "\t\t\t\"source\" : " << get_tail_from_state()->to_string() << ",\n";
        //output << "\t\t\t\"target\" : " << child->get_tail_from_state()->to_string() << ",\n";
        output << "\t\t\t\"id\" : \"" << number << "_" << child->number << "\",\n";
        output << "\t\t\t\"source\" : \"" << number << "\",\n";
        output << "\t\t\t\"target\" : \"" << child->number << "\",\n";

        output << "\t\t\t\"name\": \"" << context->alph_str(it->first) << "\",\n";
        output << "\t\t\t\"appearances\": \"";
        data->print_transition_label(output, it->first, context);
        output << "\"}\n";
    }
};

void apta::print_json(iostream& output){
    int count = 0;
    root->depth = 0;
    for(merged_APTA_iterator Ait = merged_APTA_iterator(root); *Ait != 0; ++Ait){
        apta_node* n = *Ait;
        n->number = count++;
        if(n != root) n->depth = n->source->find()->depth + 1;
    }

    output << "{\n";
    output << "\t\"nodes\" : [\n";
    bool first = true;
    for(red_state_iterator Ait = red_state_iterator(root); *Ait != 0; ++Ait){
        apta_node* n = *Ait;
        if(!first)
            output << ",\n";
        else
            first = false;
        n->print_json(output);
    }
    for(blue_state_iterator Ait = blue_state_iterator(root); *Ait != 0; ++Ait){
        apta_node* n = *Ait;
        if(n->number > 0)
            output << ",\n";
        n->print_json(output);
    }
    output << "\n\t],\n";

    output << "\t\"edges\" : [\n";
    first = true;
    for(red_state_iterator Ait = red_state_iterator(root); *Ait != 0; ++Ait){
        apta_node* n = *Ait;
        if(!first)
            output << ",\n";
        else
            first = false;
        n->print_json_transitions(output);
    }
    output << "\n\t]\n}\n";
};

void apta::print_sinks_json(iostream& output){
    /*int count = 0;
    for(merged_APTA_iterator Ait = merged_APTA_iterator(root); *Ait != 0; ++Ait){
        apta_node* n = *Ait;
        n->number = count++;
    }*/

    output << "{\n";
    output << "\t\"nodes\" : [\n";
    bool first = true;
    for(merged_APTA_iterator Ait = merged_APTA_iterator(root); *Ait != 0; ++Ait) {
        apta_node *n = *Ait;
        if(n->red == true) continue;

        if (!first) output << ",\n";
        else first = false;

        n->print_json(output);
    }
    output << "\n\t],\n";

    output << "\t\"edges\" : [\n";
    first = true;
    for(merged_APTA_iterator Ait = merged_APTA_iterator(root); *Ait != 0; ++Ait) {
        apta_node *n = *Ait;
        if(n->red == true) continue;

        if (!first) output << ",\n";
        else first = false;

        n->print_json_transitions(output);
    }
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

void apta_guard::initialize(apta_guard* g){
    target = 0;
    undo = 0;

    min_attribute_values.clear();
    max_attribute_values.clear();

    min_attribute_values.insert(g->min_attribute_values.begin(), g->min_attribute_values.end());
    max_attribute_values.insert(g->max_attribute_values.begin(), g->max_attribute_values.end());
}

void apta_node::add_tail(tail* t){
    t->next_in_list = tails_head;
    tails_head = t;
    //data->add_tail(t);
};

apta_node::apta_node(apta *context) {
    this->context = context;
    source = 0;
    original_source = 0;
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
    original_source = 0;
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

void apta_node::initialize(apta_node* n){
    source = 0;
    original_source = 0;
    representative = 0;
    next_merged_node = 0;
    representative_of = 0;
    tails_head = 0;
    label = 0;
    number = -1;
    satnumber = 0;
    colour = 0;
    size = 0;
    final = 0;
    depth = 0;
    type = -1;
    red = false;
    data->initialize();
    for(guard_map::iterator it = guards.begin(); it != guards.end(); ++it){
        mem_store::delete_guard(it->second);
    }
    guards.clear();
}

apta_node* apta_node::child(tail* t){
        int symbol = inputdata::get_symbol(t);
        for(guard_map::iterator it = guards.lower_bound(symbol); it != guards.upper_bound(symbol); ++it){
            if(it->first != symbol) break;
            apta_guard* g = it->second;
            bool outside_range = false;
            for(bound_map::iterator it2 = g->min_attribute_values.begin(); it2 != g->min_attribute_values.end(); ++it2){
                if(inputdata::get_value(t,it2->first) < it2->second){
                    outside_range = true;
                    break;
                }
            }
            if(outside_range) continue;
            for(bound_map::iterator it2 = g->max_attribute_values.begin(); it2 != g->max_attribute_values.end(); ++it2){
                if(inputdata::get_value(t,it2->first) >= it2->second){
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
            if(it->first != symbol) break;
            apta_guard* g2 = it->second;
            bool outside_range = false;
            for(bound_map::iterator it2 = g->min_attribute_values.begin(); it2 != g->min_attribute_values.end(); ++it2){
                bound_map::iterator it3 = g2->min_attribute_values.find(it2->first);
                if(it3 == g2->min_attribute_values.end() || (*it3).second != it2->second){
                    outside_range = true;
                    break;
                }
            }
            if(outside_range) continue;
            for(bound_map::iterator it2 = g->max_attribute_values.begin(); it2 != g->max_attribute_values.end(); ++it2){
                bound_map::iterator it3 = g2->max_attribute_values.find(it2->first);
                if(it3 == g2->max_attribute_values.end() || (*it3).second != it2->second){
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
        if(it->second->bounds_satisfy(t)){
            return it->second;
        }
    }
    return 0;
};

int apta_node::compute_depth(){
    int max_depth = 1;
    apta_node* n = source;
    while(n != 0){
        n = n->find();
        max_depth += 1;
        n = n->source;
    }
    if(representative_of != 0 && source != 0){
        for(n = representative_of; n!= 0; n = n->next_merged_node){
            if(n->source != 0 && n->source->find() != source->find()){
                int dist = n->compute_depth();
                if(dist > max_depth) max_depth = dist;
            }
        }
    }
    return max_depth;
};

int apta_node::depth_distance(apta_node* right){
    int min_depth = 0;
    apta_node* l = this;
    apta_node* r = right;
    while(l != 0 && r != 0){
        l = l->find();
        r = r->find();

        if(l == r) break;

        if(r->depth > l->depth) r = r->source;
        else l = l->source;
        min_depth++;
    }
    if(representative_of != 0 && source != 0){
        for(apta_node* n = representative_of; n!= 0; n = n->next_merged_node){
            if(n->source != 0 && n->source->find() != source->find()){
                int dist = n->depth_distance(right);
                if(dist < min_depth) min_depth = dist;
            }
        }
    }
    return min_depth;
};

int apta_node::num_distinct_sources(){
    int num = 0;
    if(representative_of != 0 && source != 0){
        for(apta_node* n = representative_of; n!= 0; n = n->next_merged_node){
            if(n->source != 0 && n->source->find() != source->find()){
                num++;
            }
        }
    }
    return num;
};

/* iterators for the APTA and merged APTA */

APTA_iterator::APTA_iterator(apta_node* start){
    base = start;
    current = start;
}

void APTA_iterator::increment() {
    guard_map::iterator it;
    for(it = current->guards.begin();it != current->guards.end(); ++it) {
        apta_node *target = it->second->target;
        if (target != 0 && target->source == current) {
            q.push(target);
        }
    }
    if(!q.empty()) {
        current = q.front();
        q.pop();
    } else {
        current = 0;
    }
}

merged_APTA_iterator::merged_APTA_iterator(apta_node* start){
    base = start;
    current = start;
}

void merged_APTA_iterator::increment() {
    guard_map::iterator it;
    for(it = current->guards.begin();it != current->guards.end(); ++it) {
        apta_node *target = it->second->target;
        if (target != 0 && target->source->find() == current && target->representative == 0) {
            q.push(target);
        }
    }
    if(!q.empty()) {
        current = q.front();
        q.pop();
    } else {
        current = 0;
    }
}

merged_APTA_iterator_func::merged_APTA_iterator_func(apta_node* start, bool(*node_check)(apta_node*)) : merged_APTA_iterator(start){
    check_function = node_check;
}

void merged_APTA_iterator_func::increment() {
    guard_map::iterator it;
    for(it = current->guards.begin();it != current->guards.end(); ++it) {
        apta_node *target = it->second->target;
        if (target != 0 && target->source->find() == current && target->representative == 0) {
            if(!check_function(current)) q.push(target);
        }
    }
    if(!q.empty()) {
        current = q.front();
        q.pop();
    } else {
        current = 0;
    }
}

blue_state_iterator::blue_state_iterator(apta_node* start) : merged_APTA_iterator(start) {
    if(current->red == true) increment();
}

void blue_state_iterator::increment() {
    if(current->red == true) {
        guard_map::iterator it;
        for (it = current->guards.begin(); it != current->guards.end(); ++it) {
            apta_node *target = it->second->target;
            if (target != 0 && target->source->find() == current && target->representative == 0) {
                q.push(target);
            }
        }
    }
    if(!q.empty()) {
        current = q.front();
        q.pop();
        if(current->red == true) increment();
    } else {
        current = 0;
    }
}

red_state_iterator::red_state_iterator(apta_node* start) : merged_APTA_iterator(start) { }

void red_state_iterator::increment() {
    if(current->red == true) {
        guard_map::iterator it;
        for (it = current->guards.begin(); it != current->guards.end(); ++it) {
            apta_node *target = it->second->target;
            if (target != 0 && target->source->find() == current && target->representative == 0) {
                if(target->red == true) q.push(target);
            }
        }
    }
    if(!q.empty()) {
        current = q.front();
        q.pop();
    } else {
        current = 0;
    }
}

/* OLD
APTA_iterator::APTA_iterator(apta_node* start){
    base = start;
    current = start;
}

apta_node* APTA_iterator::next_forward() {
    guard_map::iterator it;
    for(it = current->guards.begin();it != current->guards.end(); ++it){
        apta_node* target = it->second->target;
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
        while(it->second->target != current) ++it;
        ++it;
        for(; it != source->guards.end(); ++it){
            apta_node* target = it->second->target;
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
        int d = current->depth;
        apta_node* target = it->second->target;
        if(target != 0 && target->representative == 0){
            target->depth = d + 1;
            return target;
        }
    }
    return 0;
}

apta_node* merged_APTA_iterator::next_backward() {
    guard_map::iterator it;
    apta_node* source = current;
    while(source != base){
        int d = current->depth;
        current = source->find();
        if(current->depth < d - 1) current->depth = d - 1;
        source = source->source->find();
        it = source->guards.find(current->label);
        it = source->guards.begin();
        while(it->second->target != current){
            ++it;
            if(it == source->guards.end()) return 0;
        }
        ++it;
        for(; it != source->guards.end(); ++it){
            apta_node* target = it->second->target;
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
*/

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
    else if (current->next_merged_node != 0 && base != current) current = current->next_merged_node;
    else {
        while(current != 0 && current->next_merged_node == 0) current = current->representative;
        if(current != 0 && base != current) current = current->next_merged_node;
        else current = 0;
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


apta_node* tail_iterator::next_forward(){
    current_tail = current_tail->next_in_list;
    if(current_tail != 0 && current_tail->split_to != 0) increment();

    while(current_tail == 0){
        if(current == 0) return 0;
        next_node();
        if(current == 0) return 0;
        current_tail = current->tails_head;
    }
    if(current_tail != 0 && current_tail->split_to != 0) increment();

    return current;
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

//std::set<void*> freed;

apta::~apta(){
    delete root;
}

apta_node::~apta_node(){
    for(guard_map::iterator it = guards.begin();it != guards.end(); ++it){
        //if (freed.find(it->second->target) != freed.end()) {
        //    freed.insert(it->second->target);
        if(it->second->target != 0 && it->second->target->source == this)
            delete it->second->target;
        delete it->second;
        //}
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