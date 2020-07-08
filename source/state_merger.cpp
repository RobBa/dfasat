#include "state_merger.h"
#include "evaluate.h"
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <stdio.h>
#include <gsl/gsl_cdf.h>

#include "parameters.h"
#include "common.h"

state_merger::state_merger()
: context(){
    context.merger = this;
}

state_merger::state_merger(evaluation_function* e, apta* a)
: context(){
    context.merger = this;
    aut = a;
    eval = e;
    reset();
}

void state_merger::reset(){
    red_states.clear();
    blue_states.clear();
    red_states.insert(aut->root);
    update_red_blue();
}

/* BEGIN get special state sets, these are used by the SAT encoding
 * red and blue sets can be accessed directly
 */
state_set& state_merger::get_candidate_states(){
    state_set states = blue_states;
    state_set* candidate_states = new state_set();
    for(state_set::iterator it = blue_states.begin();it != blue_states.end();++it){
        for(merged_APTA_iterator Ait = merged_APTA_iterator(*it); *Ait != 0; ++Ait){
            candidate_states->insert(*Ait);
        }
    }
    return *candidate_states;
}

state_set& state_merger::get_sink_states(){
    state_set states = blue_states;
    state_set* sink_states = new state_set();
    for(state_set::iterator it = blue_states.begin();it != blue_states.end();++it){
        if(sink_type(*it) != -1){
            for(merged_APTA_iterator Ait = merged_APTA_iterator(*it); *Ait != 0; ++Ait){
                sink_states->insert(*Ait);
            }
        }
    }
    return *sink_states;
}

int state_merger::get_final_apta_size(){
    int result = 0;
    for(merged_APTA_iterator Ait = merged_APTA_iterator(aut->root); *Ait != 0; ++Ait){
        result++;
    }
    return result;
}
/* END get special state sets, these are used by the SAT encoding */

/* BEGIN basic state merging routines, these can be accessed directly, for instance to compute a conflict graph
 * the search routines do not access these methods directly, but use the perform and test merge routines below */

/* standard merge, the process is interupted when an inconsistency is found */
void state_merger::pre_split(apta_node* left, apta_node* right){
    /*
    todot();
    std::ostringstream oss2;
    oss2 << "before_presplit.dot";
    ofstream output(oss2.str().c_str());
    output << dot_output;
    output.close();
    cerr << "---- pre splitting " << left->number << " " << right->number << endl;
    for(guard_map::iterator it = right->guards.begin(); it != right->guards.end(); ++it){
        cerr << it->first << " " << it->second << endl;
    }
    */
    
    for(split_list::reverse_iterator it = left->performed_splits.rbegin(); it != left->performed_splits.rend(); ++it){
        //apta_node* child = right->child(it->first);
        //if(child != 0) perform_split(child, it->first, it->second);
        perform_split(right, it->first, it->second);
    }
    
    
    /*
    for(guard_map::iterator it = right->guards.begin(); it != right->guards.end(); ++it){
        cerr << it->first << " " << it->second << endl;
    }

    cerr << "----" << endl;
    todot();
    std::ostringstream oss;
    oss << "after_presplit.dot";
    ofstream output2(oss.str().c_str());
    output2 << dot_output;
    output2.close();
     */
}

void state_merger::undo_pre_split(apta_node* left, apta_node* right){
    /*
    todot();
    std::ostringstream oss2;
    oss2 << "before_presplit_undo.dot";
    ofstream output(oss2.str().c_str());
    output << dot_output;
    output.close();

    cerr << "----" << endl;
    for(guard_map::iterator it = right->guards.begin(); it != right->guards.end(); ++it){
        cerr << it->first << " " << it->second << endl;
    }
     * */
    
    while(!right->performed_splits.empty()){
        split_list::iterator it = right->performed_splits.begin();
        //apta_node* child = right->child(it->first);
        //if(child != 0)
        undo_perform_split(right, it->first, it->second);
    }
    
    
    /*
    for(guard_map::iterator it = right->guards.begin(); it != right->guards.end(); ++it){
        cerr << it->first << " " << it->second << endl;
    }
    cerr << "---- pre splitting undone" << left->number << " " << right->number << endl;
            
    todot();
    std::ostringstream oss;
    oss << "after_presplit_undo.dot";
    ofstream output2(oss.str().c_str());
    output2 << dot_output;
    output2.close();
    for(merged_APTA_iterator Ait = merged_APTA_iterator(aut->root); *Ait != 0; ++Ait){
        apta_node* n = *Ait;
    }
     */
}

/* old tryout:
 bool found = false;
 for(guard_map::iterator it = left->guards.begin(); it != left->guards.end(); ++it){
 apta_guard* g = it->second;
 if(!g->min_attribute_values.empty()){ found = true; break; }
 if(!g->max_attribute_values.empty()){ found = true; break; }
 }
 if(!found) return;
 
 guard_map old_guards = guard_map(right->guards);
 right->guards = guard_map();
 for(guard_map::iterator it = left->guards.begin();it != left->guards.end(); ++it){
 apta_guard* new_guard = new apta_guard(it->second);
 int symbol = it->first;
 new_guard->target = old_guards[symbol];
 if(old_guards[symbol] != 0 && !new_guard->min_attribute_values.empty())
 new_guard->target = new apta_node();
 right->guards.insert(pair<int,apta_guard*>(it->first,new_guard));
 }
 
 for(guard_map::iterator it = old_guards.begin(); it != old_guards.end(); ++it){
 apta_guard* old_guard = it->second;
 if(old_guard->target != 0){
 apta_node* node = old_guard->target;
 for(tail_iterator it2 = tail_iterator(node); *it2 != 0; ++it2){
 tail* old_tail = *it2;
 apta_node* new_node = right->child(old_tail);
 if(node != new_node){
 tail* new_tail = t->split();
 new_node->add_tail(new_tail);
 new_node->size++;
 old_node->size--;
 }
 }
 }
 }
 */

bool state_merger::merge(apta_node* left, apta_node* right, int depth, bool evaluate, bool perform, bool test){
    bool result = true;
    if(left == 0 || right == 0) return true;

    if(test) {
        // test early stopping conditions
        if (KTAIL != -1 && depth > KTAIL) return true;
        if (KSTATE != -1 && (left->size < KSTATE || right->size < KSTATE)) return true;

        if (evaluate && eval->consistent(this, left, right) == false) return false;

        if (left->red && RED_FIXED) {
            for (guard_map::iterator it = right->guards.begin(); it != right->guards.end(); ++it) {
                apta_node *child = it->second->target;
                if (child == 0) continue;
                if (left->child(it->first) == 0) return false;
                //&& !it->second->target->data->sink_consistent(it->second->target, 0)) return false;
            }
        }
    }
    
    if(perform){
        pre_split(left, right);
    }
    if(evaluate){
        eval->update_score(this, left, right);
    }
    if(perform){
        right->merge_with(left);
        left->data->update(right->data);
    }
    if(evaluate) {
        eval->update_score_after(this, left, right);
    }

    for(guard_map::iterator it = right->guards.begin(); it != right->guards.end(); ++it){
        if(it->second->target == 0) continue;
        int i = it->first;
        apta_guard* right_guard = it->second;
        apta_guard* left_guard = left->guard(i, it->second);
        
        if(left_guard == 0 || left_guard->target == 0){
            if(perform){
                left->set_child(i, right_guard->target);
            }
        } else {
            apta_node* right_child = right_guard->target;
            apta_node* left_child = left_guard->target;
            
            apta_node* child = left_child->find();
            apta_node* other_child = right_child->find();
            
            if(child != other_child){
                if(perform){
                    other_child->set_undo(i, right);
                }
                result = merge(child, other_child, depth + 1, evaluate, perform, test);
                if(test && result == false) break;
            }
        }
    }
    return result;
}
bool state_merger::merge(apta_node* left, apta_node* right) {
    return merge(left, right, 0, true, true, true);
}

/* forced merge, continued even when inconsistent */
void state_merger::merge_force(apta_node* left, apta_node* right){
    merge(left, right, 0, false, true, false);
    /*
    if(left == 0 || right == 0) return;
    
    pre_split(left, right);
    
    right->merge_with(left);
    left->data->update(right->data);
    //right->representative = left;
    //left->size += right->size;
    
    for(guard_map::iterator it = right->guards.begin();it != right->guards.end(); ++it){
        if(it->second->target == 0) continue;
        int i = it->first;
        apta_guard* right_guard = it->second;
        apta_guard* left_guard = left->guard(i, it->second);
        
        if(left_guard == 0 || left_guard->target == 0){
            left->set_child(i, right_guard->target);
        } else {
            apta_node* right_child = right_guard->target;
            apta_node* left_child = left_guard->target;
            
            apta_node* child = left_child->find();
            apta_node* other_child = right_child->find();
            
            if(child != other_child){
                other_child->set_undo(i, right);
                merge_force(child, other_child);
            }
        }
    }
     */
}

/* testing merge, no actual merge is performed, only consistency and score are computed */
bool state_merger::merge_test(apta_node* left, apta_node* right){
    return merge(left, right, 0, true, false, true);
    /*bool result = true;
    if(left == 0 || right == 0) return true;
    if(eval->consistent(this, left, right) == false) return false;
    
    if(left->red && RED_FIXED){
        for(guard_map::iterator it = right->guards.begin();it != right->guards.end(); ++it)
            if(left->child(it->first) == 0) return false;
            //&& !it->second->target->data->sink_consistent(it->second->target, 0)) return false;
    }
    eval->update_score(this, left, right);
    eval->update_score_after(this, left, right);
    for(guard_map::iterator it = right->guards.begin();it != right->guards.end(); ++it){
        if(it->second->target == 0) continue;
        int i = it->first;
        apta_guard* right_guard = it->second;
        apta_guard* left_guard = left->guard(i, it->second);
        
        if(left_guard != 0 && left_guard->target != 0){
            apta_node* right_child = right_guard->target;
            apta_node* left_child = left_guard->target;
            
            apta_node* child = left_child->find();
            apta_node* other_child = right_child->find();
            
            if(child != other_child){
                result = merge_test(child, other_child, depth + 1);
                if(result == false) break;
            }
        }
    }
    return result;*/
}

/* undo merge, works for both forcing and standard merging, not needed for testing merge */
void state_merger::undo_merge(apta_node* left, apta_node* right){
    
    if(left == 0 || right == 0) return;
    
    if(right->representative != left) return;
    
    for(guard_map::reverse_iterator it = right->guards.rbegin();it != right->guards.rend(); ++it){
        if(it->second->target == 0) continue;
        int i = it->first;
        apta_guard* right_guard = it->second;
        apta_guard* left_guard = left->guard(i, it->second);
        if(left_guard == 0) continue;
        
        apta_node* right_child = right_guard->target;
        apta_node* left_child = left_guard->target;
        if(left_child == right_child){
            left_guard->target = 0;//->set_child(i,0);
        } else if(left_child != 0){
            //apta_node* other_child = right_child->find_until(right, i);
            apta_node* child = right_child->representative;
            if(child != right_child){
                undo_merge(child, right_child);
            }
            right_child->set_undo(i,0);
        }
    }
    
    right->undo_merge_with(left);
    left->data->undo(right->data);
    //left->size -= right->size;
    //right->representative = 0;
    //eval->undo_update(this, left, right);
    undo_pre_split(left, right);
}

bool state_merger::split_single(apta_node* new_node, apta_node* old_node, tail* t, bool evaluate){
    //if(old_node->size == 0) return true;
    
    int symbol = inputdata::get_symbol(t);
    if(symbol == -1) return true;
    //cerr << "++  " << inputdata::alphabet[symbol] << endl;
    apta_node* old_child = old_node->child(symbol)->find();
    apta_node* new_child = new_node->child(symbol);
    if(new_child == 0){
        new_child = new apta_node(aut);
        new_node->set_child(symbol,new_child);
        new_child->source = new_node;
        new_child->label  = symbol;
        //new_child->number = node_number++;
        new_child->depth = new_node->depth + 1;
    }
    //cerr << "childs " << old_child << " " << new_child << endl;
    
    if(t->future() == 0) return true;
    
    if(evaluate) eval->split_update_score_before(this, old_child, new_child, t->future());
    new_child->size++;
    old_child->size--;
    //old_child->data->split_update_single(new_child->data, t->future());
    
    old_child->data->del_tail(t->future());
    tail* new_tail = t->future()->split();
    new_child->add_tail(new_tail);
    new_child->data->add_tail(new_tail);    
    split_single(new_child, old_child, new_tail, evaluate);
    
    if(evaluate) eval->split_update_score_after(this, old_child, new_child, t->future());
    return true;
};

void state_merger::undo_split_single(apta_node* new_node, apta_node* old_node){
    //if(new_node == old_node) return;
    
    for(guard_map::iterator it = new_node->guards.begin(); it != new_node->guards.end(); ++it){
        if(it->second->target == 0) continue;
        int i = it->first;
        apta_guard* new_guard = it->second;
        apta_guard* old_guard = old_node->guard(i, it->second);
        
        apta_node* new_child = new_guard->target->find();
        apta_node* old_child = old_guard->target->find();
        
        undo_split_single(new_child, old_child);
        
        delete new_child;
    }
    
    old_node->size = old_node->size + new_node->size;
    for(tail* t = new_node->tails_head; t != 0; t = t->next()){
        t->split_from->undo_split();
    }
    
    old_node->data->split_undo(new_node->data);
};

/* new_node and old_node are already split
 * this function performs the split recursively on their children */
bool state_merger::split(apta_node* new_node, apta_node* old_node){
    if(old_node->size == 0 && false){
        /*
        new_node->guards = guard_map(old_node->guards);
        for(guard_map::iterator it = new_node->guards.begin(); it != new_node->guards.end(); ++it){
            if(it->second->target != 0) it->second->target->source = new_node;
        }
        old_node->guards.clear();
         */
        old_node->data->split_undo(new_node->data);
        old_node->size = old_node->size + new_node->size;
        for(tail* t = new_node->tails_head; t != 0; t = t->next()){
            t->split_from->undo_split();
        }
        delete new_node;
    } else {
        for (tail *t = new_node->tails_head; t != 0; t = t->next()) {
            if (t->get_index() == -1) continue;
            if (t->future() == 0) continue;

            int symbol = inputdata::get_symbol(t);
            apta_node *old_child = old_node->child(symbol)->find();
            apta_node *new_child = new_node->child(symbol);
            if (new_child == 0) {
                new_child = new apta_node(aut);
                new_node->set_child(symbol, new_child);
                new_child->source = new_node;
                new_child->label = symbol;
                new_child->depth = new_node->depth + 1;
            }
            new_child->size++;

            //old_child->size--;
            //old_child->data->split_update_single(new_child->data, t->future());

            tail *new_tail = t->future()->split();
            new_child->add_tail(new_tail);
            new_child->data->add_tail(new_tail);

            /*
            new_child->size++;
            //new_child->data->read_to(t);
            //new_child->data->read_from(t->future());
            tail* new_tail = t->future()->split();
            new_child->add_tail(new_tail);
            */
        }
        for (guard_map::iterator it = new_node->guards.begin(); it != new_node->guards.end(); ++it) {
            if (it->second->target == 0) continue;

            apta_guard *new_guard = it->second;
            apta_guard *old_guard = old_node->guard(it->first);

            apta_node* new_child = new_guard->target;
            apta_node* old_child = old_guard->target->find();

            if(new_child->size == old_child->size){
                for (tail *t = new_child->tails_head; t != 0; t = t->next()) {
                    t->split_from->undo_split();
                }
                delete new_child;
                new_guard->target = old_child;
                old_guard->target = 0;

                if(old_child->original_source == 0) old_child->original_source = old_child->source;
                old_child->source = new_node;
            } else {
                old_child->size -= new_child->size;
                old_child->data->split_update(new_child->data);
                split(new_child, old_child);
            }
        }
        return true;
    }
};

/* new_node and old_node are not yet unsplit
 * this function first performs the unsplit recursively on their children */
void state_merger::undo_split(apta_node* new_node, apta_node* old_node){
    if(old_node->size == 0 && false){
        old_node->guards = guard_map(new_node->guards);
        for(guard_map::iterator it = old_node->guards.begin(); it != old_node->guards.end(); ++it){
            if(it->second->target != 0){
                it->second->target->source = old_node;
                apta_node* oi = old_node->next_merged_node;
                while(oi != 0) {
                    for (guard_map::iterator it2 = oi->guards.begin(); it2 != oi->guards.end(); ++it2){
                        if(it2->second->target == it->second->target)
                            it->second->target->source = oi;
                    }
                    oi = oi->next_merged_node;
                }
            }
        }
    } else {
        for(guard_map::reverse_iterator it = new_node->guards.rbegin(); it != new_node->guards.rend(); ++it){
            if(it->second->target == 0) continue;
            int i = it->first;
            apta_guard* new_guard = it->second;
            apta_node *new_child = new_guard->target;
            apta_guard* old_guard = old_node->guard(i, it->second);

            if(old_guard->target == 0){
                old_guard->target = new_child;
                new_guard->target = 0;
                if(new_child->original_source->find() == old_node){
                    new_child->source = new_child->original_source;
                } else {
                    new_child->source = old_node;
                }
                continue;
            }
            apta_node *old_child = old_guard->target->find();

            undo_split(new_child, old_child);

            old_child->data->split_undo(new_child->data);
            old_child->size = old_child->size + new_child->size;
            for (tail *t = new_child->tails_head; t != 0; t = t->next()) {
                t->split_from->undo_split();
            }
            delete new_child;
        }
    }
        
    /*old_node->size = old_node->size + new_node->size;
    for(tail* t = new_node->tails_head; t != 0; t = t->next()){
        t->split_from->undo_split();
    }
    
    old_node->data->split_undo(new_node->data);*/
};

/* old version using split single
 * 
 * void state_merger::perform_split(apta_node* blue, tail* t, int attr){
    int symbol = inputdata::get_symbol(t);
    float val  = inputdata::get_value(t, attr);
    
    apta_node* red = blue->source->find();
    apta_guard* old_guard = red->guard(t);
    
    red->performed_splits.push_front(pair< tail*, int >(t, attr));
    
    if(old_guard != 0){
        apta_guard* new_guard = new apta_guard(old_guard);
        
        old_guard->min_attribute_values[attr] = val;
        new_guard->max_attribute_values[attr] = val;
        
        apta_node* new_child = new apta_node(aut);
        
        for(tail_iterator it = tail_iterator(blue); *it != 0; ++it){
            tail* t2 = *it;
            if(inputdata::get_value(t2->past_tail, attr) < val){
                //cerr << "++ " << symbol << " " << inputdata::get_value(t2->past_tail, attr) << endl;

                blue->data->split_update_single(new_child->data, t2);

                tail* new_tail = t2->split();
                new_child->add_tail(new_tail);
                new_child->size++;
                //new_child->data->read_from(new_tail);
                //new_child->data->read_to(new_tail->past_tail);
                blue->size--;

                split_single(new_child, blue, new_tail, false);
                //eval->split_update_score_after(this, blue, new_child, t);
            
            } else {
                //cerr << "-- " << symbol << " " << inputdata::get_value(t2->past_tail, attr) << endl;
            }
        }
        new_guard->target = new_child;
        new_child->source = red;
        red->guards.insert(pair<int, apta_guard*>(symbol, new_guard));
        //if(blue->source != red) red->guards.insert(pair<int, apta_guard*>(symbol, new_guard));
        new_child->label = symbol;
        
        //cerr << "attr " << attr << " " << inputdata::get_value(t, attr) << endl; 
        //cerr << "old: " << old_guard << " " << old_guard->target << " " << blue << endl;
        //cerr << "new: " << new_guard << " " << new_guard->target << " " << new_child << endl;
        
        //eval->update_score_after(this, blue, new_child);
        //split(new_child, blue);
        //eval->update_score(this, blue, new_child);
    } else {
        cerr << "empty guard?" << endl;
    }
}*/

void state_merger::perform_split(apta_node* red, tail* t, int attr){
    int symbol = inputdata::get_symbol(t);
    float val  = inputdata::get_value(t, attr);
    
    //apta_node* red = blue->source->find();
    apta_guard* old_guard = red->guard(t);
    if(old_guard != 0){
        red->performed_splits.push_front(pair< tail*, int >(t, attr));

        apta_node* blue = old_guard->target;

        apta_guard* new_guard = new apta_guard(old_guard);
        old_guard->min_attribute_values[attr] = val;
        new_guard->max_attribute_values[attr] = val;
        red->guards.insert(pair<int, apta_guard *>(symbol, new_guard));

        if(blue == 0) return;

        int split_count = 0;
        for(tail_iterator it = tail_iterator(blue); *it != 0; ++it) {
            tail *t2 = *it;
            if (inputdata::get_value(t2->past_tail, attr) < val) split_count++;
        }

        if(split_count == 0){
            return;
        }
        if(split_count == blue->size){
            new_guard->target = blue;
            old_guard->target = 0;
            /*if(new_guard->target->original_source == 0){
                new_guard->target->original_source = new_guard->target->source;
            }
            new_guard->target->source = blue;*/
        } else {
            apta_node *new_child = new apta_node(aut);

            for (tail_iterator it = tail_iterator(blue); *it != 0; ++it) {
                tail *t2 = *it;
                if (inputdata::get_value(t2->past_tail, attr) < val) {
                    //cerr << "++ " << symbol << " " << inputdata::get_value(t2->past_tail, attr) << endl;

                    //blue->data->split_update_single(new_child->data, t2);

                    tail *new_tail = t2->split();
                    new_child->add_tail(new_tail);
                    new_child->data->add_tail(new_tail);
                    new_child->size++;
                    //new_child->data->read_from(new_tail);
                    //new_child->data->read_to(new_tail->past_tail);
                    //blue->size--;

                    //split_single(new_child, blue, new_tail, false);
                    //eval->split_update_score_after(this, blue, new_child, t);

                    /*tail* new_tail = t2->split();
                    new_child->add_tail(new_tail);
                    new_child->size++;
                    new_child->data->read_from(new_tail);
                    new_child->data->read_to(new_tail->past_tail);
                    blue->size--;*/
                }
            }

            new_guard->target = new_child;
            new_child->source = red;
            new_child->label = symbol;
            //if(blue->source != red) red->guards.insert(pair<int, apta_guard*>(symbol, new_guard));


            //cerr << "attr " << attr << " " << inputdata::get_value(t, attr) << endl;
            //cerr << "old: " << old_guard << " " << old_guard->target << " " << blue << endl;
            //cerr << "new: " << new_guard << " " << new_guard->target << " " << new_child << endl;

            //eval->update_score_after(this, blue, new_child);
            blue->data->split_update(new_child->data);
            blue->size -= new_child->size;
            split(new_child, blue);
            //eval->update_score(this, blue, new_child);
        }
    } else {
        //cerr << "empty guard?" << endl;
    }
}

void state_merger::undo_perform_split(apta_node* red, tail* t, int attr){
    int symbol = inputdata::get_symbol(t);
    float val  = inputdata::get_value(t, attr);
    
    //apta_node* source = blue->source->find();
    apta_guard* old_guard = red->guard(t);
    
    guard_map::iterator new_it = red->guards.upper_bound(symbol);
    new_it--;
    apta_guard* new_guard = (*new_it).second;

    if(old_guard != new_guard){
        if(new_guard->target == 0) {
            // do nothing
        }
        else if(old_guard->target == 0){
                old_guard->target = new_guard->target;
                new_guard->target = 0;
                /*if(old_guard->target->original_source->find() == red) {
                old_guard->target->source = old_guard->target->original_source;
                } else {
                old_guard->target->source = red;
                }*/
        } else {
            apta_node* blue = old_guard->target;
            apta_node *new_child = new_guard->target;

            //cerr << "attr " << attr << " " << inputdata::get_value(t, attr) << endl;
            //cerr << "undo old " << old_guard << " " << old_guard->target << endl;
            //cerr << "undo new " << new_guard << " " << new_child << endl;

            //undo_split_single(new_child, blue);
            undo_split(new_child, blue);
            blue->size += new_child->size;
            blue->data->split_undo(new_child->data);
            for (tail *t = new_child->tails_head; t != 0; t = t->next()) {
                t->split_from->undo_split();
            }
            delete new_child;
        }

        if (new_guard->min_attribute_values.find(attr) != new_guard->min_attribute_values.end())
            old_guard->min_attribute_values[attr] = new_guard->min_attribute_values[attr];
        else
            old_guard->min_attribute_values.erase(attr);

        red->guards.erase(new_it);
        delete new_guard;

            /*
             if(blue->source->representative != 0){
             source = blue->source->find();
             old_guard = source->guard(t);

             guard_map::iterator new_it = source->guards.upper_bound(symbol);
             new_it--;
             new_guard = (*new_it).second;
             source->guards.erase(new_it);
             delete new_guard;
             }
             */
        red->performed_splits.erase(red->performed_splits.begin());
    }
}

/* END basic state merging routines */

/* update red blue sets after performing a merge, we keep these in memory instead of recomputing them */
void state_merger::update_red_blue(){
    state_set new_red;
    state_set new_blue;
    for(state_set::iterator it = red_states.begin(); it != red_states.end(); ++it){
        apta_node* red = (*it)->find();
        new_red.insert(red);
        red->red = true;
    }
    for(state_set::iterator it = new_red.begin(); it != new_red.end(); ++it){
        apta_node* red = *it;
        for(int i = 0; i < alphabet_size; ++i){
            apta_node* child = red->get_child(i);
            if(child != 0 && new_red.find(child) == new_red.end()) new_blue.insert(child);
        }
    }
    red_states = new_red;
    blue_states = new_blue;
    eval->update(this);
}

/* BEGIN merge functions called by state merging algorithms */

/* make a given blue state red, and its children blue */
void state_merger::extend(apta_node* blue){
    blue->red = true;
    return;
    blue_states.erase(blue);
    red_states.insert(blue);
    
    for(guard_map::iterator it = blue->guards.begin(); it != blue->guards.end(); ++it){
        apta_node* child = it->second->target;
        if(child == 0) continue;
        blue_states.insert(child);
    }
}

/* undo making a given blue state red */
void state_merger::undo_extend(apta_node* blue){
    blue->red = false;
    return;
    
    blue_states.insert(blue);
    red_states.erase(blue);
    
    for(guard_map::iterator it = blue->guards.begin(); it != blue->guards.end(); ++it){
        apta_node* child = it->second->target;
        blue_states.erase(child);
    }
}

/* make the first blue state (ordered by size) red, and its children blue */
apta_node* state_merger::extend_red(){
    
    apta_node* top_state = 0;
    for(blue_state_iterator it = blue_state_iterator(aut->root); *it != 0; ++it){
        if(!MERGE_SINKS_DSOLVE && (sink_type(*it) != -1)) continue;
        if(top_state == 0 || top_state->size < (*it)->size) top_state = *it;
    }
    if(top_state != 0) extend(top_state);
    return top_state;
    
    state_set blues = state_set();
    for(blue_state_iterator it = blue_state_iterator(aut->root); *it != 0; ++it) blues.insert(*it);
    
    for(state_set::iterator it = blues.begin(); it != blues.end(); ++it){
        //for(state_set::iterator it = blue_states.begin(); it != blue_states.end(); ++it){
        apta_node* blue = *it;
        bool found = false;
        
        if(!MERGE_SINKS_DSOLVE && (sink_type(blue) != -1)) continue;
        
        for(red_state_iterator it2 = red_state_iterator(aut->root); *it2 != 0; ++it2){
            //for(state_set::iterator it2 = red_states.begin(); it2 != red_states.end(); ++it2){
            apta_node* red = *it2;
            refinement* ref = test_merge(red, blue);
            if(ref != 0){ found = true; break; }
        }
        
        if(found == false) {
            extend(blue);
            return blue;
        }
    }
    return 0;
}

/* perform a merge, assumed to succeed, no testing for consistency or score computation */
void state_merger::perform_merge(apta_node* left, apta_node* right){
    merge_force(left, right);
    if(STORE_MERGES) left->size_store.push_back(pair<int,int>(num_merges,right->size));
    num_merges++;
    update_red_blue();
}

/* undo a merge, assumed to succeed, no testing for consistency or score computation */
void state_merger::undo_perform_merge(apta_node* left, apta_node* right){
    undo_merge(left, right);
    if(STORE_MERGES) left->size_store.pop_back();
    update_red_blue();
}

refinement* state_merger::get_stored_merge(apta_node* left, apta_node* right){
    if(!STORE_MERGES) return 0;

    score_map::iterator it = right->eval_store.find(left);
    if(it != right->eval_store.end()){
        int merge_time = it->second.first.first;
        int merge_size = it->second.first.second;
        score_pair merge_score = it->second.second;

        int size_same = 0;
        for(size_list::reverse_iterator it2 = left->size_store.rbegin(); it2 != left->size_store.rend(); ++it2){
            if(it2->first <= merge_time){
                size_same = it2->second;
                break;
            }
        }
        int size_change = left->size - size_same;

        if(STORE_MERGES_KEEP_CONFLICT && size_same == merge_size && merge_score.first == false)
            return new merge_refinement(merge_score.second, left, right);

        if(size_change < STORE_MERGES_SIZE_THRESHOLD)
            return new merge_refinement(merge_score.second, left, right);

        if((double)size_change/(double)left->size < STORE_MERGES_RATIO_THRESHOLD)
            return new merge_refinement(merge_score.second, left, right);
    }
    return 0;
}

void state_merger::store_merge(bool merge_consistent, double merge_score, apta_node* left, apta_node* right){
    if(STORE_MERGES){
        right->eval_store[left] = ts_pair(pair<int,int>(num_merges,left->size),score_pair(merge_consistent, merge_score));
    }
}

/* test a merge, behavior depending on input parameters
 * it performs a merge, computes its consistency and score, and undos the merge
 * returns a <consistency,score> pair */
refinement* state_merger::test_merge(apta_node* left, apta_node* right){
    eval->reset(this);

    double score_result = -1;
    bool   merge_result = false;

    if(MERGE_LOCAL){
        if(left->depth_distance(right) > MERGE_LOCAL_MAXDIST && left->num_distinct_sources() > MERGE_LOCAL_COLLECTOR_COUNT){
            return 0;
        }
    }
    
    if(eval->compute_before_merge) score_result = eval->compute_score(this, left, right);
    
    if(MERGE_WHEN_TESTING){
        merge_result = merge(left,right);
    } else {
        merge_result = merge_test(left,right);
    }
    
    if(merge_result && !eval->compute_before_merge) score_result = eval->compute_score(this, left, right);
    
    if((merge_result && eval->compute_consistency(this, left, right) == false)) merge_result = false;
    
    if(USE_LOWER_BOUND && score_result < LOWER_BOUND) merge_result = false;
    
    if(MERGE_WHEN_TESTING) undo_merge(left,right);
    
    if(merge_result == false) return 0;
    return new merge_refinement(score_result, left, right);
}

int testnr = 1;

refinement_set* state_merger::test_splits(apta_node* blue){
    refinement_set* result = new refinement_set();
    double score = 0;
    for(int attr = 0; attr < inputdata::num_attributes; ++attr){
        eval->reset_split(this, blue);
        multimap<float, tail*> sorted_tails;
        for(tail_iterator it = tail_iterator(blue); *it != 0; ++it){
            tail* t = *it;
            sorted_tails.insert(pair<float, tail*>(inputdata::get_value(t->past_tail, attr),t));
        }
        
        apta_node* new_node = new apta_node(aut);
        
        float prev_val = (*(sorted_tails.begin())).first;
        //cerr << "splits: attr:" << attr << endl;
        for(multimap<float, tail*>::iterator it = sorted_tails.begin(); it != sorted_tails.end(); ++it){
            //cerr << "testing split " << it->first << endl;
            //cerr << it->first << " " << prev_val << endl;
            if(it->first > prev_val){
                score = eval->split_compute_score(this, blue, new_node);
                if(eval->split_compute_consistency(this, blue, new_node)){
                    result->insert(new split_refinement(score, blue->source->find(), it->second->past_tail, attr));
                }
                //cerr << score << " ";
                prev_val = it->first;
            }
            tail* t = it->second;

            eval->split_update_score_before(this, blue, new_node, t);
            //blue->data->split_update_single(new_node->data, t);
            
            blue->size--;
            blue->data->del_tail(t);
            tail* new_tail = t->split();
            new_node->size++;
            new_node->add_tail(new_tail);
            new_node->data->add_tail(new_tail);
            //cerr << "++" << endl;
            //cerr << "++  " << inputdata::alphabet[blue->label] << " <= " << prev_val << endl;
            split_single(new_node, blue, new_tail, true);
            eval->split_update_score_after(this, blue, new_node, t);
        }
        //cerr << endl;
        
        //cerr << "undoing split " << endl;
        undo_split_single(new_node, blue);
        
        delete new_node;
        sorted_tails.clear();
    }
    
    return result;
}

/* returns all possible refinements given the current sets of red and blue states
 * behavior depends on input parameters
 * the merge score is used as key in the returned refinement_set */
refinement_set* state_merger::get_possible_refinements(){
    refinement_set* result = new refinement_set();
    
    state_set blue_its = state_set();
    //bool found = false;
    
    for(blue_state_iterator it = blue_state_iterator(aut->root); *it != 0; ++it){
        if((*it)->size != 0) blue_its.insert(*it);
    }
    // DEBUG("checking for " << blue_its.size() << " blue states" <<endl);
    
    for(state_set::iterator it = blue_its.begin(); it != blue_its.end(); ++it){
        apta_node* blue = *it;
        bool found = false;
        if((sink_type(blue) != -1)) continue;
        
        // cerr << inputdata::num_attributes << endl;
        if(inputdata::num_attributes > 0){
            //cerr << "testing splits" << endl;
            refinement_set* refs = test_splits(blue);
            if(refs != 0){
                result->insert(refs->begin(), refs->end());
                found = true;
                delete refs;
            }
        }
        
        for(red_state_iterator it2 = red_state_iterator(aut->root); *it2 != 0; ++it2){
            apta_node* red = *it2;
            
            refinement* ref = test_merge(red,blue);
            
            if(ref != 0){
                result->insert(ref);
                found = true;
            }
        }
        
        if(MERGE_BLUE_BLUE){
            for(state_set::iterator it2 = blue_its.begin(); it2 != blue_its.end(); ++it2){
                apta_node* blue2 = *it2;
                
                if(blue == blue2) continue;
                
                if(!MERGE_SINKS_DSOLVE && (sink_type(blue2) != -1)) continue;
                
                refinement* ref = test_merge(blue2,blue);
                if(ref != 0){
                    result->insert(ref);
                    //mset->insert(pair<int, pair<apta_node*, apta_node*> >(score.second, pair<apta_node*, apta_node*>(blue2, blue)));
                    found = true;
                }
            }
        }
        
        // cerr << "found results: " << result->size() << endl;
        // why not result->size() == 0?
        if(found == false && sink_type(blue) == -1) {
            
            // cerr << "came out empty for merges" <<endl;
            
            if(EXTEND_ANY_RED){
                for(refinement_set::iterator it = result->begin(); it != result->end(); ++it) delete *it;
                result->clear();
                result->insert(new extend_refinement(blue));
                return result;
            }
            
        }
        
        //if(result->size() == 0)
        result->insert(new extend_refinement(blue));
        
        if(MERGE_MOST_VISITED) break;
    }
    return result;
}

// streaming version: relevant threshold count from Hoeffding bound and using epsilon to determine whether
// we indeed do have a biggest score
merge_map* state_merger::get_possible_merges(int count){
    merge_map* mset = new merge_map();
    
    for(red_state_iterator it = red_state_iterator(aut->root); *it != 0; ++it){
        //for(state_set::iterator it = red_states.begin(); it != red_states.end(); ++it){
        apta_node* red = *(it);
        
        if(red->size < count) continue;
        
        for(blue_state_iterator it2 = blue_state_iterator(aut->root); *it2 != 0; ++it2){
            //for(state_set::iterator it2 = blue_states.begin(); it2 != blue_states.end(); ++it2){
            apta_node* blue = *it2;
            if(blue->size < count) continue;
            
            if(!MERGE_SINKS_DSOLVE && (sink_type(blue) != -1)) continue;
            
            refinement* ref = test_merge(red,blue);
            if(ref != 0){
                mset->insert(pair<int, pair<apta_node*, apta_node*> >(ref->score, pair<apta_node*, apta_node*>(red, blue)));
            }
            
            if(MERGE_MOST_VISITED) break;
            blue->age = 0;
            
            if(MERGE_BLUE_BLUE){
                for(blue_state_iterator it3 = blue_state_iterator(aut->root); *it3 != 0; ++it3){
                    //for(state_set::iterator it2 = blue_states.begin(); it2 != blue_states.end(); ++it2){
                    apta_node* blue2 = *it3;
                    
                    if(blue->size < count) continue;
                    if(blue == blue2) continue;
                    
                    if(!MERGE_SINKS_DSOLVE && (sink_type(blue2) != -1)) continue;
                    
                    refinement* ref = test_merge(blue2,blue);
                    if(ref != 0){
                        mset->insert(pair<int, pair<apta_node*, apta_node*> >(ref->score, pair<apta_node*, apta_node*>(blue2, blue)));
                    }
                }
            }
        }
    }
    return mset;
}

/* returns the highest scoring merge given the current sets of red and blue states
 * behavior depends on input parameters
 * returns (0,0) if none exists (given the input parameters) */
refinement* state_merger::get_best_refinement(){
    refinement_set* rs = get_possible_refinements();
    refinement* r = 0;
    if(!rs->empty()){
        r = *rs->begin();
    }
    delete rs;
    return r;

    refinement* best = 0;
    
    state_set blue_its = state_set();
    for(blue_state_iterator it = blue_state_iterator(aut->root); *it != 0; ++it){
        if((*it)->size != 0) blue_its.insert(*it);
    }
    
    for(state_set::iterator it = blue_its.begin(); it != blue_its.end(); ++it){
        apta_node* blue = *it;
        bool found = false;
        
        if(!MERGE_SINKS_DSOLVE && (sink_type(blue) != -1)) continue;
        
        //cerr << inputdata::num_attributes << endl;
        if(inputdata::num_attributes > 0){
            //cerr << "testing splits" << endl;
            test_splits(blue);
        }
        
        for(red_state_iterator it2 = red_state_iterator(aut->root); *it2 != 0; ++it2){
            apta_node* red = *it2;
            
            refinement* ref = test_merge(red,blue);
            if(ref != 0 && (best == 0 || best->score < ref->score)){
                delete best;
                best = ref;
            } else {
                delete ref;
            }
        }
        
        if(MERGE_BLUE_BLUE){
            for(state_set::iterator it2 = blue_its.begin(); it2 != blue_its.end(); ++it2){
                apta_node* blue2 = *it2;
                
                if(blue == blue2) continue;
                
                refinement* ref = test_merge(blue2,blue);
                if(ref != 0 && (best == 0 || best->score < ref->score)){
                    delete best;
                    best = ref;
                } else {
                    delete ref;
                }
            }
        }
        
        if(found == false && sink_type(blue) == -1) {
            if(EXTEND_ANY_RED){
                delete best;
                return new extend_refinement(blue);
            }
            if(best == 0 || best->score < -1){
                delete best;
                best = new extend_refinement(blue);
            }
        }
        
        if(MERGE_MOST_VISITED) break;
    }
    return best;
}

merge_pair* state_merger::get_best_merge(int count){
    merge_pair* best = new merge_pair(0,0);
    double result = -1;
    
    for(blue_state_iterator it = blue_state_iterator(aut->root); *it != 0; ++it){
        //for(state_set::iterator it = blue_states.begin(); it != blue_states.end(); ++it){
        apta_node* blue = *(it);
        if(!MERGE_SINKS_DSOLVE && (sink_type(blue) != -1)) continue;
        
        for(red_state_iterator it2 = red_state_iterator(aut->root); *it2 != 0; ++it2){
            //for(state_set::iterator it2 = red_states.begin(); it2 != red_states.end(); ++it2){
            apta_node* red = *it2;
            
            refinement* ref = test_merge(red,blue);
            if(ref != 0 && ref->score > result){
                best->first = red;
                best->second = blue;
                result = ref->score;
            } else {
                delete ref;
            }
        }
        
        if(MERGE_BLUE_BLUE){
            for(blue_state_iterator it2 = blue_state_iterator(aut->root); *it2 != 0; ++it2){
                //for(state_set::iterator it2 = blue_states.begin(); it2 != blue_states.end(); ++it2){
                apta_node* blue2 = *it2;
                
                if(blue == blue2) continue;
                
                refinement* ref = test_merge(blue2,blue);
                if(ref != 0 && ref->score > result){
                    best->first = blue2;
                    best->second = blue;
                    result = ref->score;
                } else {
                    delete ref;
                }
            }
        }
        
        if(MERGE_MOST_VISITED) break;
    }
    return best;
}

/* END merge functions called by state merging algorithms */

/* input function 	        *
 * pass along to  eval fct      */

// streaming mode methods
void state_merger::init_apta(string data) {
    eval->init(data, this);
}

void state_merger::advance_apta(string data) {
    eval->add_sample(data, this);
}

// batch mode methods
void state_merger::read_apta(istream &input_stream){
    aut->read_file(input_stream);
}

void state_merger::read_apta(string dfa_file){
    ifstream input_stream(dfa_file);
    read_apta(input_stream);
    input_stream.close();
}

void state_merger::read_apta(vector<string> dfa_data){
    stringstream input_stream;
    for (int i = 0; i < dfa_data.size(); i++) {
        input_stream << dfa_data[i];
    }
    //eval->read_file(input_stream, this);
    read_apta(input_stream);
}

/* output functions */
void state_merger::todot(){
    stringstream dot_output_buf;
    aut->print_dot(dot_output_buf);
    dot_output = "// produced with flexfringe from git commit"  + string(gitversion) + '\n' + "// " + COMMAND + '\n'+ dot_output_buf.str();
}

void state_merger::tojson(){
    stringstream output_buf;
    aut->print_json(output_buf);
    json_output = output_buf.str(); // json does not support comments, maybe we need to introduce a meta section
}

void state_merger::tojsonsinks(){
    stringstream output_buf;
    aut->print_sinks_json(output_buf);
    json_output = output_buf.str(); // json does not support comments, maybe we need to introduce a meta section
}

void state_merger::print_json(FILE* output)
{
    fprintf(output, "%s", json_output.c_str());
}


void state_merger::print_dot(FILE* output)
{
    fprintf(output, "%s", dot_output.c_str());
}

int state_merger::sink_type(apta_node* node){
    return node->data->sink_type(node);
};

bool state_merger::sink_consistent(apta_node* node, int type){
    return node->data->sink_consistent(type);
};

int state_merger::num_sink_types(){
    return 0; //eval->num_sink_types();
}; // this got moved to eval data */

int state_merger::compute_global_score(){
    return get_final_apta_size();
};
