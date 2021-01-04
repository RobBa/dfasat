#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <stdio.h>
//#include <gsl/gsl_cdf.h>

#include "state_merger.h"
#include "evaluate.h"
#include "depth-driven.h"
#include "alergia.h"

#include "parameters.h"

int EVAL_TYPE = 1;

REGISTER_DEF_DATATYPE(alergia_data);
REGISTER_DEF_TYPE(alergia);

alergia_data::alergia_data() : count_data() {
	// LOG_S(INFO) << "Creating alergia_data class"; // example; warning, this is created in every APTA node
};

void alergia_data::initialize() {
    count_data::initialize();
    trans_counts.clear();
}

void alergia_data::read_from(int type, int index, int length, int symbol, string data){
    count_data::read_from(type, index, length, symbol, data);

    //if(!PERTYPE_DISTRIBUTIONS) type = 0;

    if(trans_counts.find(type) == trans_counts.end()){
        //trans_counts[type] = num_map();
        trans_counts[type][symbol] = 1;
    } else {
        if(trans_counts[type].find(symbol) == trans_counts[type].end()){
            trans_counts[type][symbol] = 1;
        } else {
            trans_counts[type][symbol]++;
        }
    }
};

void alergia_data::add_tail(tail* t){
    count_data::add_tail(t);

    if(t->get_index() == -1) return;

    int type = t->get_type();
    int symbol = t->get_symbol();
    if(trans_counts.find(type) == trans_counts.end()){
        //trans_counts[type] = num_map();
        trans_counts[type][symbol] = 1;
    } else {
        if(trans_counts[type].find(symbol) == trans_counts[type].end()){
            trans_counts[type][symbol] = 1;
        } else {
            trans_counts[type][symbol]++;
        }
    }
}
void alergia_data::del_tail(tail* t){
    count_data::del_tail(t);

    if(t->get_index() == -1) return;

    int type = t->get_type();
    //if(!PERTYPE_DISTRIBUTIONS) type = 1;
    int symbol = t->get_symbol();

    trans_counts[type][symbol]--;
}

void alergia_data::print_transition_label(iostream& output, int symbol, apta* apta_context){
    for(type_num_map::iterator it = trans_counts.begin(); it != trans_counts.end(); it++){
        output << it->second[symbol] << " ";
    }
};

void alergia_data::print_state_label(iostream& output, apta* aptacontext){
    count_data::print_state_label(output, aptacontext);
    output << "\n" << num_paths() << " " << num_final();
};

void alergia_data::update(evaluation_data* right){
    count_data::update(right);
    alergia_data* other = (alergia_data*)right;
    for(type_num_map::iterator it = other->trans_counts.begin(); it != other->trans_counts.end(); it++){
        num_map& that_nm = it->second;
        num_map& this_nm = trans_counts[it->first];
        for(num_map::iterator itnm = that_nm.begin(); itnm != that_nm.end(); ++itnm){
            this_nm[(*itnm).first] += (*itnm).second;
        }
    }
};

void alergia_data::undo(evaluation_data* right){
    count_data::undo(right);
    alergia_data* other = (alergia_data*)right;
    for(type_num_map::iterator it = other->trans_counts.begin(); it != other->trans_counts.end(); it++){
        num_map& that_nm = it->second;
        num_map& this_nm = trans_counts[it->first];
        for(num_map::iterator itnm = that_nm.begin(); itnm != that_nm.end(); ++itnm){
            this_nm[(*itnm).first] -= (*itnm).second;
        }
    }
};

inline void alergia::update_divider(double left_count, double right_count, double& left_divider, double& right_divider){
    if(left_count >= SYMBOL_COUNT && right_count >= SYMBOL_COUNT){
        left_divider += left_count + CORRECTION;
        right_divider += right_count + CORRECTION;
    }
}

inline void alergia::update_divider_pool(double left_pool, double right_pool, double& left_divider, double& right_divider){
    if(left_pool >= SYMBOL_COUNT || right_pool >= SYMBOL_COUNT){
        left_divider += left_pool + CORRECTION;
        right_divider += right_pool + CORRECTION;
    }
}

inline void alergia::update_left_pool(double left_count, double right_count, double& left_pool, double& right_pool){
    if(right_count < SYMBOL_COUNT){
        left_pool += left_count;
        right_pool += right_count;
    }
}

inline void alergia::update_right_pool(double left_count, double right_count, double& left_pool, double& right_pool){
    if(left_count < SYMBOL_COUNT){
        left_pool += left_count;
        right_pool += right_count;
    }
}

bool alergia::alergia_consistency(double right_count, double left_count, double right_total, double left_total){
    if(left_count >= SYMBOL_COUNT && right_count >= SYMBOL_COUNT){
        double bound = (1.0 / sqrt(left_total) + 1.0 / sqrt(right_total));
        bound = bound * sqrt(0.5 * log(2.0 / CHECK_PARAMETER));

        double gamma = (left_count / left_total) - (right_count / right_total);

        if(gamma > bound) return false;
        if(-gamma > bound) return false;
    }
    return true;
};

bool alergia::alergia_consistency_pool(double right_count, double left_count, double right_total, double left_total){
    if(left_count >= SYMBOL_COUNT || right_count >= SYMBOL_COUNT){
        double bound = (1.0 / sqrt(left_total) + 1.0 / sqrt(right_total));
        bound = bound * sqrt(0.5 * log(2.0 / CHECK_PARAMETER));

        double gamma = (left_count / left_total) - (right_count / right_total);

        if(gamma > bound) return false;
        if(-gamma > bound) return false;
    }
    return true;
};

bool alergia::data_consistent(alergia_data* l, alergia_data* r){
    /* we ignore low frequency states, decided by input parameter STATE_COUNT */
    if(FINAL_PROBABILITIES) if(r->num_paths() + r->num_final() < STATE_COUNT || l->num_paths() + l->num_final() < STATE_COUNT) return true;
    else if(r->num_paths() < STATE_COUNT || l->num_paths() < STATE_COUNT) return true;

    /* we treat type distributions as independent */
    for(type_num_map::iterator it = l->counts_begin(); it != l->counts_end(); it++){
        int type = it->first;
        num_map& nm = it->second;

        /* computing the dividers (denominator) */
        double left_divider = 0.0; double right_divider = 0.0;
        /* we pool low frequency counts (sum them up in a separate bin), decided by input parameter SYMBOL_COUNT
        * we create pools 1 and 2 separately for left and right low counts
        * in this way, we can detect differences in distributions even if all counts are low (i.e. [0,0,1,1] vs [1,1,0,0]) */
        double l1_pool = 0.0; double r1_pool = 0.0; double l2_pool = 0.0; double r2_pool = 0.0;

        int matching_right = 0;
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            int symbol = it2->first;
            double left_count = it2->second;
            if(left_count == 0) continue;

            double right_count = r->count(type,symbol);
            matching_right += right_count;

            update_divider(left_count, right_count, left_divider, right_divider);
            update_left_pool(left_count, right_count, l1_pool, r1_pool);
            update_right_pool(left_count, right_count, l2_pool, r2_pool);
        }
        r2_pool += r->num_paths(type) - matching_right;

        /* optionally add final probabilities (input parameter) */
        if(FINAL_PROBABILITIES){
            double left_count = l->num_final(type);
            double right_count = r->num_final(type);

            update_divider(left_count, right_count, left_divider, right_divider);
            update_left_pool(left_count, right_count, l1_pool, r1_pool);
            update_right_pool(left_count, right_count, l2_pool, r2_pool);
        }

        update_divider_pool(l1_pool, r1_pool, left_divider, right_divider);
        update_divider_pool(l2_pool, r2_pool, left_divider, right_divider);

        if(left_divider < STATE_COUNT || right_divider < STATE_COUNT) continue;

        /* we have calculated the dividers and pools */
        for(num_map::iterator it2 = nm.begin(); it2 != nm.end(); ++it2){
            int symbol = it2->first;
            double left_count = it2->second;
            if(left_count == 0) continue;

            double right_count = r->count(type,symbol);

            if(alergia_consistency(left_count, right_count, left_divider, right_divider) == false){
                inconsistency_found = true; return false;
            }
        }
        if(FINAL_PROBABILITIES){
            double left_count = l->num_final(type);
            double right_count = r->num_final(type);

            if(alergia_consistency(left_count, right_count, left_divider, right_divider) == false){
                inconsistency_found = true; return false;
            }
        }
        if(alergia_consistency_pool(l1_pool, r1_pool, left_divider, right_divider) == false){
            inconsistency_found = true; return false;
        }
        if(alergia_consistency_pool(l2_pool, r2_pool, left_divider, right_divider) == false){
            inconsistency_found = true; return false;
        }
    }
    return true;
};

/* ALERGIA, consistency based on Hoeffding bound, only uses positive (type=1) data, pools infrequent counts */
bool alergia::consistent(state_merger *merger, apta_node* left, apta_node* right){
    //if(count_driven::consistent(merger, left, right) == false){ inconsistency_found = true; return false; }
    //if(left->depth != right->depth) {inconsistency_found = true; return false;};
    alergia_data* l = (alergia_data*) left->data;
    alergia_data* r = (alergia_data*) right->data;
    
    return data_consistent(l, r);
};

/* When is an APTA node a sink state?
 * sink states are not considered merge candidates
 *
 * accepting sink = only accept, accept now, accept afterwards
 * rejecting sink = only reject, reject now, reject afterwards 
 * low count sink = frequency smaller than STATE_COUNT */
bool alergia_data::is_low_count_sink(){
    //cerr << num_paths() << " " << num_final() << "<" << SINK_COUNT << endl;
    return num_paths() + num_final() < SINK_COUNT;
    //if(EVAL_TYPE == 1) return pos_final() + neg_final() + pos_paths() + neg_paths() < SINK_COUNT;
    //return pos_paths() + neg_paths() < SINK_COUNT;
}

bool alergia_data::is_stream_sink(apta_node* node) {

   return node->size < STREAM_COUNT;
}

int alergia_data::sink_type(apta_node* node){
    if(!USE_SINKS) return -1;

    if (is_low_count_sink()) return 0;
    if(is_stream_sink(node)) return 2;
    return -1;
};

bool alergia_data::sink_consistent(int type){
    if(!USE_SINKS) return true;
    
    if(type == 0) return is_low_count_sink();
    //if(type == 1) return this->is_stream_sink();
    return true;
};

int alergia::num_sink_types(){
    if(!USE_SINKS) return 0;
    return 2;
};

/*void alergia::print_dot(iostream& output, state_merger* merger){
    apta* aut = merger->aut;
    state_set s  = merger->red_states;
    
    cerr << "size: " << s.size() << endl;
    
    output << "digraph DFA {\n";
    output << "\t" << aut->root->find()->number << " [label=\"root\" shape=box];\n";
    output << "\t\tI -> " << aut->root->find()->number << ";\n";
    for(state_set::iterator it = merger->red_states.begin(); it != merger->red_states.end(); ++it){
        apta_node* n = *it;
        alergia_data* l = (alergia_data*) n->data;

        output << "\t" << n->number << " [shape=ellipse label=\"[" << l->pos_final() << "]\\n[";

        int my_sum = 0;
        for(int i = 0; i < alphabet_size; ++i){
            my_sum += l->pos(i);
        }
        output << " " << my_sum;
        
        output << "]\"];\n";
        state_set childnodes;
        set<int> sinks;
        for(int i = 0; i < alphabet_size; ++i){
            apta_node* child = n->get_child(i);
            if(child == 0){
                // no output
            } else {
                 if(merger->sink_type(child) != -1){
                     sinks.insert(sink_type(child));
                 } else {
                     childnodes.insert(child);
                 }
            }
        }
        // should not need sinks
        for(set<int>::iterator it2 = sinks.begin(); it2 != sinks.end(); ++it2){
            int stype = *it2;
            output << "\tS" << n->number << "t" << stype << " [label=\"sink " << stype << "\" shape=box];\n";
            output << "\t\t" << n->number << " -> S" << n->number << "t" << stype << " [label=\"";
            for(int i = 0; i < alphabet_size; ++i){
                if(n->get_child(i) != 0 && sink_type(n->get_child(i)) == stype){
                    output << " " << aut->alph_str(i).c_str() << " [" << ((alergia_data*)n->data)->num_pos[i] << ":" << ((alergia_data*)n->data)->num_neg[i] << "]";

                }
            }
            output << "\"];\n";
        }
        // this is what i care about
        for(state_set::iterator it2 = childnodes.begin(); it2 != childnodes.end(); ++it2){
            apta_node* child = *it2;
            output << "\t\t" << n->number << " -> " << child->number << " [label=\"";
            for(int i = 0; i < alphabet_size; ++i){
                if(n->get_child(i) != 0 && n->get_child(i) == child){
                    output << " " << aut->alph_str(i).c_str() << " [" << ((alergia_data*)n->data)->num_pos[i] << ":" << ((alergia_data*)n->data)->num_neg[i] << "]";
                }
            }
            output << "\"];\n";
        }
    }

    // leak workaround
    state_set *state = &merger->get_candidate_states();
    
    s = *state; // merger->get_candidate_states();

    for(state_set::iterator it = s.begin(); it != s.end(); ++it){
        apta_node* n = *it;
        alergia_data* l = reinterpret_cast<alergia_data*>(n->data);
        if(l->pos_final() != 0){
            output << "\t" << n->number << " [shape=doublecircle style=dotted label=\"" << l->pos_final() << ":" << l->neg_final() << "\\n[" << l->pos_paths() << ":" << l->neg_paths() << "]\"];\n";
        } else if(l->neg_final() != 0){
            output << "\t" << n->number << " [shape=Mcircle style=dotted label=\"" << l->pos_final() << ":" << l->neg_final() << "\\n[" << l->pos_paths() << ":" << l->neg_paths() << "]\"];\n";
        } else {
            output << "\t" << n->number << " [shape=circle style=dotted label=\"0:0\\n[" << l->pos_paths() << ":" << l->neg_paths() << "]\"];\n";
        }
        for(child_map::iterator it2 = n->children.begin(); it2 != n->children.end(); ++it2){
            apta_node* child = it2->second;
            int symbol = it2->first;
            output << "\t\t" << n->number << " -> " << child->number << " [style=dotted label=\"" << aut->alph_str(symbol).c_str() << " [" << l->num_pos[symbol] << ":" << l->num_neg[symbol] << "]\"];\n";
        }
    }

    output << "}\n";

    delete state;
};*/
