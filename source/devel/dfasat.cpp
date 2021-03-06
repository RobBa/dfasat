#include "dfasat.h"
#include "random_greedy.h"
#include "evaluate.h"
//#include <malloc.h>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <set>
#include <vector>
#include <unistd.h>
#include <sys/wait.h>
#include <ctime>

bool MERGE_SINKS_PRESOLVE = 0;
int OFFSET = 1;
int EXTRA_STATES = 0;
bool TARGET_REJECTING = 0;
bool SYMMETRY_BREAKING = 0;
bool FORCING = 0;

int literal_counter = 1;
int clause_counter = 0;

bool computing_header = true;

// apta/input state i has color j
int **x;
// color i has a transition with label a to color j
int ***y;
// color i is an accepting state
int *z;

// color i has a transition with label a to sink j
int ***sy;
// state i is a sink or one of the parents of apta/input state i is a sink
int *sp;

// literals used for symmetry breaking
int **yt;
int **yp;
int **ya;

state_merger merger;
state_set red_states;
state_set non_red_states;
state_set sink_states;

FILE* sat_stream; // TODO: also one as member in state_merger

int dfa_size;
int sinks_size;
int num_states;
int new_states;
int new_init;

set<int> trueliterals;

int best_solution = -1;

void reset_literals(bool init){
    int v, i, j, a; // TODO: aren't there better names for those?

    literal_counter = 1;
    for(v = 0; v < num_states; ++v)
        for(i = 0; i < dfa_size; ++i)
            if(init || x[v][i] > 0) x[v][i] = literal_counter++;
    
    for(a = 0; a < alphabet_size; ++a)
        for(i = 0; i < dfa_size; ++i)
            for(j = 0; j < dfa_size; ++j)
                if(init || y[a][i][j] > 0) y[a][i][j] = literal_counter++;

    for(a = 0; a < alphabet_size; ++a)
        for(i = 0; i < dfa_size; ++i)
            for(j = 0; j < sinks_size; ++j)
                if(init || sy[a][i][j] > 0) sy[a][i][j] = literal_counter++;

    for(i = 0; i < num_states; ++i)
        if(init || sp[i] > 0) sp[i] = literal_counter++;
    
    for(i = 0; i < dfa_size; ++i)
        if(init || z[i] > 0) z[i] = literal_counter++;

    for(i = 0; i < dfa_size; ++i)
        for(j = 0; j < new_states; ++j)
            if(init || yt[i][j] > 0) yt[i][j] = literal_counter++;
    
    for(i = 0; i < dfa_size; ++i)
        for(j = 0; j < new_states; ++j)
            if(init || yp[i][j] > 0) yp[i][j] = literal_counter++;
    
    for(a = 0; a < alphabet_size; ++a)
        for(i = 0; i < new_states; ++i)
            if(init || ya[a][i] > 0) ya[a][i] = literal_counter++;
}

void create_literals(){
    int v, a, i;
    //X(STATE,COLOR)
    x = (int**) malloc( sizeof(int*) * num_states);
    for(v = 0; v < num_states; v++ )
        x[ v ] = (int*) malloc( sizeof(int) * dfa_size);
    
    //Y(LABEL,COLOR,COLOR)
    y = (int***) malloc( sizeof(int**) * alphabet_size);
    for(a = 0; a < alphabet_size; ++a){
        y[ a ] = (int**) malloc( sizeof(int*) * dfa_size);
        for(i = 0; i < dfa_size; ++i)
            y[ a ][ i ]  = (int*) malloc( sizeof(int) * dfa_size);
    }

    //SY(LABEL,COLOR,SINK)
    sy = (int***) malloc( sizeof(int**) * alphabet_size);
    for(a = 0; a < alphabet_size; ++a){
        sy[ a ] = (int**) malloc( sizeof(int*) * dfa_size);
        for(i = 0; i < dfa_size; ++i)
            sy[ a ][ i ]  = (int*) malloc( sizeof(int) * sinks_size);
    }

    //SP(STATE)
    sp = (int*) malloc( sizeof(int) * num_states);
    
    //Z(COLOR)
    z = (int*) malloc( sizeof(int) * dfa_size);
    
    //YT(COLOR,COLOR)
    yt = (int**) malloc( sizeof(int*) * dfa_size);
    for(i = 0; i < dfa_size; ++i)
        yt[ i ]  = (int*) malloc( sizeof(int) * new_states);
    
    //YP(COLOR,COLOR)
    yp = (int**) malloc( sizeof(int*) * dfa_size);
    for(i = 0; i < dfa_size; ++i)
        yp[ i ]  = (int*) malloc( sizeof(int) * new_states);
    
    //YA(LABEL,COLOR)
    ya = (int**) malloc( sizeof(int*) * alphabet_size);
    for(a = 0; a < alphabet_size; ++a)
        ya[ a ]  = (int*) malloc( sizeof(int) * new_states);
    
    // reset literal values
    reset_literals(true);
}

void delete_literals(){
    int v, a, i;
    for(v = 0; v < num_states; v++ )
        free(x[ v ]);
    free(x);
    for(a = 0; a < alphabet_size; ++a){
        for(i = 0; i < dfa_size; ++i)
            free(y[ a ][ i ]);
        free(y[ a ]);
    }
    free(y);
    for(a = 0; a < alphabet_size; ++a){
        for(i = 0; i < sinks_size; ++i)
            free(sy[ a ][ i ]);
        free(sy[ a ]);
    }
    free(sy);
    free(z);
    free(sp);
    for(i = 0; i < dfa_size; ++i)
        free(yt[ i ]);
    free(yt);
    for(i = 0; i < dfa_size; ++i)
        free(yp[ i ]);
    free(yp);
    for(a = 0; a < alphabet_size; ++a)
        free(ya[ a ]);
    free(ya);
}

/* Print clauses without eliminated literals, -2 = false, -1 = true */
int print_clause(bool v1, int l1, bool v2, int l2, bool v3, int l3, bool v4, int l4){
    if(v1 && l1 == -1) return 0;
    if(!v1 && l1 == -2) return 0;
    if(v2  && l2 == -1) return 0;
    if(!v2 && l2 == -2) return 0;
    if(v3 && l3 == -1) return 0;
    if(!v3 && l3 == -2) return 0;
    if(v4 && l4 == -1) return 0;
    if(!v4 && l4 == -2) return 0;
    
    if(computing_header) return 1;
    
    if(v1 == true  && l1 != -2) fprintf(sat_stream, "%i ", l1);
    if(v1 == false && l1 != -1) fprintf(sat_stream, "-%i ", l1);
    if(v2 == true  && l2 != -2) fprintf(sat_stream, "%i ", l2);
    if(v2 == false && l2 != -1) fprintf(sat_stream, "-%i ", l2);
    if(v3 == true  && l3 != -2) fprintf(sat_stream, "%i ", l3);
    if(v3 == false && l3 != -1) fprintf(sat_stream, "-%i ", l3);
    if(v4 == true  && l4 != -2) fprintf(sat_stream, "%i ", l4);
    if(v4 == false && l4 != -1) fprintf(sat_stream, "-%i ", l4);
    
    fprintf(sat_stream, " 0\n");
    return 1;
}

int print_clause(bool v1, int l1, bool v2, int l2, bool v3, int l3){
    if(v1 && l1 == -1) return 0; // TODO: use nullptr?
    if(!v1 && l1 == -2) return 0;
    if(v2  && l2 == -1) return 0;
    if(!v2 && l2 == -2) return 0;
    if(v3 && l3 == -1) return 0;
    if(!v3 && l3 == -2) return 0;
    
    if(computing_header) return 1;
    
    if(v1 == true  && l1 != -2) fprintf(sat_stream, "%i ", l1);
    if(v1 == false && l1 != -1) fprintf(sat_stream, "-%i ", l1);
    if(v2 == true  && l2 != -2) fprintf(sat_stream, "%i ", l2);
    if(v2 == false && l2 != -1) fprintf(sat_stream, "-%i ", l2);
    if(v3 == true  && l3 != -2) fprintf(sat_stream, "%i ", l3);
    if(v3 == false && l3 != -1) fprintf(sat_stream, "-%i ", l3);
    
    fprintf(sat_stream, " 0\n");
    return 1;
}

int print_clause(bool v1, int l1, bool v2, int l2){
    if(v1 && l1 == -1) return 0;
    if(!v1 && l1 == -2) return 0;
    if(v2  && l2 == -1) return 0;
    if(!v2 && l2 == -2) return 0;
    
    if(computing_header) return 1;
    
    if(v1 == true  && l1 != -2) fprintf(sat_stream, "%i ", l1);
    if(v1 == false && l1 != -1) fprintf(sat_stream, "-%i ", l1);
    if(v2 == true  && l2 != -2) fprintf(sat_stream, "%i ", l2);
    if(v2 == false && l2 != -1) fprintf(sat_stream, "-%i ", l2);
    
    fprintf(sat_stream, " 0\n");
    return 1;
}

bool always_true(int number, bool flag){
    if(number == -1 && flag == true)  return true;
    if(number == -2 && flag == false) return true;
    return false;
}

void print_lit(int number, bool flag){
    if(computing_header) return;
    if(number < 0) return;
    
    if(flag == true) fprintf(sat_stream, "%i ", number);
    else fprintf(sat_stream, "-%i ", number);
}

void print_clause_end(){
    if(computing_header) return;
    fprintf(sat_stream, " 0\n");
}

/* fix values for red states -2 = false, -1 = true */
void fix_red_values(){
    for(state_set::iterator it = red_states.begin();it != red_states.end();++it){
        apta_node* node = *it;
        
        for(int i = 0; i < dfa_size; ++i) x[node->satnumber][i] = -2;
        sp[node->satnumber] = -2;
        x[node->satnumber][node->colour] = -1;
        
        apta_node* source = *it;
        for(int label = 0; label < alphabet_size; ++label){
            apta_node* target = source->get_child(label);
            if(target != 0 && red_states.find(target) != red_states.end()){
                for(int i = 0; i < dfa_size; ++i) y[label][source->colour][i] = -2;
                for(int i = 0; i < sinks_size; ++i) sy[label][source->colour][i] = -2;
                y[label][source->colour][target->colour] = -1;
            } else if(MERGE_SINKS_PRESOLVE && target != 0 && sink_states.find(target) != sink_states.end()){
                for(int i = 0; i < dfa_size; ++i) y[label][source->colour][i] = -2;
                for(int i = 0; i < sinks_size; ++i)
                    if(sink_consistent(target, i) == false) sy[label][source->colour][i] = -2;
            } else if(TARGET_REJECTING && target == 0){
                for(int i = 0; i < dfa_size; ++i) y[label][source->colour][i] = -2;
                for(int i = 0; i < sinks_size; ++i) sy[label][source->colour][i] = -2;
            }
        }
        
        if(node->pos_final() != 0) z[node->colour] = -1;
        if(node->neg_final() != 0) z[node->colour] = -2;
    }
}

/* erase possible colors due to symmetry reduction
 should be compatible with BFS symmtry breaking, unchecked */
int set_symmetry(){
    int num = 0;
    int max_value = new_init;
    for(state_set::iterator it = red_states.begin(); it != red_states.end(); ++it){
        if(max_value + 1>= dfa_size)
            break;
        
        apta_node* node = *it;
        for(int a = 0; a < alphabet_size; ++a){
            if(max_value + 1>= dfa_size)
                break;
            
            apta_node* child = (*it)->get_child(a);
            if(child != 0 && red_states.find(child) == red_states.end()){
                if(MERGE_SINKS_PRESOLVE && (is_accepting_sink(child) || is_rejecting_sink(child)))
                    continue;
                
                for(int i = max_value + 1; i < dfa_size; ++i){
                    x[child->satnumber][i] = -2;
                }
                max_value++;
            }
        }
    }
    return num;
}

int print_symmetry(){
    int num = 0;
    for(int i = 0; i < dfa_size; ++i){
        for(int k = 0; k < new_states; ++k){
            for(int j = 0; j < i; ++j){
                for(int l = k + 1; l < new_states; l++){
                    num += print_clause(false, yp[i][k], false, yp[j][l]);
                }
            }
        }
    }
    for(int i = 0; i < dfa_size; ++i){
        for(int k = 0; k < new_states; ++k){
            for(int l = k + 1; l < new_states; l++){
                for(int a = 0; a < alphabet_size; ++a){
                    for(int b = 0; b < a; ++b){
                        num += print_clause(false, yp[i][k], false, yp[i][l], false, ya[a][k], false, ya[b][l]);
                    }
                }
            }
        }
    }
    return num;
}

/* eliminate literals for merges that conflict with the red states */
void erase_red_conflict_colours(){
    for(state_set::iterator it = red_states.begin(); it != red_states.end(); ++it){
        apta_node* left = *it;
        for(state_set::iterator it2 = non_red_states.begin(); it2 != non_red_states.end(); ++it2){
            apta_node* right = *it2;
            if(merger.test_merge(left,right) == -1) x[right->satnumber][left->colour] = -2;
            if(right->pos_paths() != 0 || right->pos_final() != 0) x[right->satnumber][0] = -2;
            if(right->neg_paths() != 0 || right->neg_final() != 0) x[right->satnumber][1] = -2;
        }
    }
}

/* print the at least one en at most one clauses for x */
int print_colours(){
    int num = 0;
    bool altr = false;
    // at least one
    for(state_set::iterator it = non_red_states.begin(); it != non_red_states.end(); ++it){
        apta_node* node = *it;
        altr = always_true(sp[node->satnumber], true);
        for(int k = 0; k < dfa_size; ++k){
            if(altr) break;
            altr = always_true(x[node->satnumber][k], true);
        }
        if(altr == false){
            for(int k = 0; k < dfa_size; ++k)
                print_lit(x[node->satnumber][k], true);
            print_lit(sp[node->satnumber], true);
            print_clause_end();
            num += 1;
        }
    }
    // at most one
    for(state_set::iterator it = non_red_states.begin(); it != non_red_states.end(); ++it){
        apta_node* node = *it;
        for(int a = 0; a < dfa_size; ++a)
            for(int b = a+1; b < dfa_size; ++b)
                num += print_clause(false, x[node->satnumber][a], false, x[node->satnumber][b]);
        for(int a = 0; a < dfa_size; ++a)
            num += print_clause(false, x[node->satnumber][a], false, sp[node->satnumber]);
    }
    return num;
}

/* print clauses restricting two unmergable states to have the same color *
 * excludes pairs of states that are covered by the z literals            */
int print_conflicts(){
    int num = 0;
    for(state_set::iterator it = non_red_states.begin(); it != non_red_states.end(); ++it){
        apta_node* left = *it;
        state_set::iterator it2 = it;
        ++it2;
        while(it2 != non_red_states.end()){
            apta_node* right = *it2;
            ++it2;
            if(left->pos_final() != 0 && right->neg_final() != 0) continue;
            if(left->neg_final() != 0 && right->pos_final() != 0) continue;
            
            if(merger.test_merge(left, right) == -1){
                for(int k = 0; k < dfa_size; ++k)
                    num += print_clause(false, x[left->satnumber][k], false, x[right->satnumber][k]);
            }
        }
    }
    return num;
}

/* print de clauses voor z literals */
int print_accept(){
    int num = 0;
    for(state_set::iterator it = non_red_states.begin(); it != non_red_states.end(); ++it){
        apta_node* node = *it;
        for(int k = 0; k < dfa_size; ++k){
            if(node->pos_final() != 0) num += print_clause(false, x[node->satnumber][k], true, z[k]);
            if(node->neg_final() != 0) num += print_clause(false, x[node->satnumber][k], false, z[k]);
        }
    }
    return num;
}

/* print de clauses voor y literals */
int print_transitions(){
    int num = 0;
    for(int a = 0; a < alphabet_size; ++a)
        for(int i = 0; i < dfa_size; ++i)
            for(int j = 0; j < dfa_size; ++j)
                for(int h = 0; h < j; ++h)
                    num += print_clause(false, y[a][i][h], false, y[a][i][j]);
    for(int a = 0; a < alphabet_size; ++a)
        for(int i = 0; i < dfa_size; ++i)
            for(int j = 0; j < sinks_size; ++j)
                for(int h = 0; h < dfa_size; ++h)
                    num += print_clause(false, y[a][i][h], false, sy[a][i][j]);
    for(int a = 0; a < alphabet_size; ++a)
        for(int i = 0; i < dfa_size; ++i)
            for(int j = 0; j < sinks_size; ++j)
                for(int h = 0; h < j; ++h)
                    num += print_clause(false, sy[a][i][h], false, sy[a][i][j]);
    return num;
}

/* print transitions for any label yt */
int print_t_transitions(){
    int num = 0;
    for(int i = 0; i < dfa_size; ++i)
        for(int j = 0; j < new_states; ++j)
            for(int a = 0; a < alphabet_size; ++a)
                num += print_clause(false, y[a][i][new_init+j], true, yt[i][j]);
    
    for(int i = 0; i < dfa_size; ++i){
        for(int j = 0; j < new_states; ++j){
            bool altr = false;
            for(int a = 0; a < alphabet_size; ++a)
                if(y[a][i][new_init+j] == -1) altr = true;
            if(!altr){
                if(!computing_header){
                    print_lit(yt[i][j], false);
                    for(int a = 0; a < alphabet_size; ++a){
                        print_lit(y[a][i][new_init+j], true);
                    }
                    print_clause_end();
                }
                num++;
            }
        }
    }
    
    return num;
}

/* print BFS tree transitions */
int print_p_transitions(){
    int num = 0;
    for(int i = 0; i < dfa_size; ++i){
        for(int j = 0; j < new_states; ++j){
            for(int k = 0; k < i; ++k){
                num += print_clause(false, yp[i][j], false, yt[k][j]);
            }
            num += print_clause(false, yp[i][j], true, yt[i][j]);
        }
    }
    for(int i = 0; i < new_states; ++i){
        bool altr = false;
        for(int j = 0; j < new_init+i; ++j)
            if(yp[j][i] == -1) altr = true;
        if(!altr){
            if(!computing_header){
                for(int j = 0; j < new_init+i; ++j){
                    print_lit(yp[j][i], true);
                }
                print_clause_end();
            }
            num++;
        }
    }
    return num;
}

/* print BFS tree labels */
int print_a_transitions(){
    int num = 0;
    for(int i = 0; i < new_states; ++i){
        for(int a = 0; a < alphabet_size; ++a){
            for(int j = 0; j < dfa_size; ++j){
                for(int b = 0; b < a; ++b){
                    num += print_clause(false, ya[a][i], false, yp[j][i], false, y[b][j][new_init+i]);
                }
                num += print_clause(false, ya[a][i], false, yp[j][i], true, y[a][j][new_init+i]);
            }
        }
    }
    for(int i = 0; i < new_states; ++i){
        bool altr = false;
        for(int a = 0; a < alphabet_size; ++a){
            if(ya[a][i] == -1) altr = true;
        }
        if(!altr){
            if(!computing_header){
                for(int a = 0; a < alphabet_size; ++a){
                    print_lit(ya[a][i], true);
                }
                print_clause_end();
            }
            num++;
        }
    }
    return num;
}

/* print de clauses voor y literals */
int print_forcing_transitions(){
    int num = 0;
    bool altr = false;
    for (int label = 0; label < alphabet_size; ++label) {
        state_set label_states;
        for (state_set::iterator it = red_states.begin(); it != red_states.end(); ++it) {
            apta_node* source = *it;
            if(source->get_child(label) != 0) label_states.insert(source);
        }
        for (state_set::iterator it = non_red_states.begin(); it != non_red_states.end(); ++it) {
            apta_node* source = *it;
            if(source->get_child(label) != 0) label_states.insert(source);
        }
        
        for(int i = 0; i < dfa_size; ++i){
            for(int j = 0; j < dfa_size; ++j){
                altr = always_true(y[label][i][j], false);
                for (state_set::iterator it = label_states.begin(); it != label_states.end(); ++it) {
                    if(altr) break;
                    apta_node* source = *it;
                    altr = always_true(x[source->satnumber][i], true);
                }
                if(altr == false){
                    for (state_set::iterator it = label_states.begin(); it != label_states.end(); ++it) {
                        apta_node* source = *it;
                        print_lit(x[source->satnumber][i], true);
                    }
                    print_lit(y[label][i][j], false);
                    print_clause_end();
                    num += 1;
                }
            }
        }
        
        for(int i = 0; i < dfa_size; ++i){
            for(int j = 0; j < sinks_size; ++j){
                altr = always_true(sy[label][i][j], false);
                for (state_set::iterator it = label_states.begin(); it != label_states.end(); ++it) {
                    if(altr) break;
                    apta_node* source = *it;
                    altr = always_true(x[source->satnumber][i], true);
                }
                if(altr == false){
                    for (state_set::iterator it = label_states.begin(); it != label_states.end(); ++it) {
                        apta_node* source = *it;
                        print_lit(x[source->satnumber][i], true);
                    }
                    print_lit(sy[label][i][j], false);
                    print_clause_end();
                    num += 1;
                }
            }
        }
    }
    return num;
}

/* print de determinization constraint */
int print_paths(){
    int num = 0;
    for (state_set::iterator it = red_states.begin(); it != red_states.end(); ++it) {
        apta_node* source = *it;
        for (int label = 0; label < alphabet_size; ++label) {
            apta_node* target = source->get_child(label);
            if (target != 0 && sink_states.find(target) == sink_states.end()) {
                for (int i = 0; i < dfa_size; ++i)
                    for (int j = 0; j < dfa_size; ++j)
                        num += print_clause(true, y[label][i][j], false, x[source->satnumber][i], false, x[target->satnumber][j]);
            }
        }
    }
    for (state_set::iterator it = non_red_states.begin(); it != non_red_states.end(); ++it) {
        apta_node* source = *it;
        for (int label = 0; label < alphabet_size; ++label) {
            apta_node* target = source->get_child(label);
            if (target != 0) {
                for (int i = 0; i < dfa_size; ++i)
                    for (int j = 0; j < dfa_size; ++j)
                        num += print_clause(true, y[label][i][j], false, x[source->satnumber][i], false, x[target->satnumber][j]);
            }
        }
    }
    return num;
}

/* print sink paths */
int print_sink_paths(){
    int num = 0;
    bool altr = false;
    for (state_set::iterator it = red_states.begin(); it != red_states.end(); ++it) {
        apta_node* source = *it;
        for (int label = 0; label < alphabet_size; ++label) {
            apta_node* target = source->get_child(label);
            if (target != 0 && sink_states.find(target) == sink_states.end()) {
                for (int i = 0; i < dfa_size; ++i)
                    for (int j = 0; j < sinks_size; ++j)
                        if(sink_consistent(target, j) == false)
                            num += print_clause(false, sy[label][i][j], false, x[source->satnumber][i]);
                
                for (int i = 0; i < dfa_size; ++i){
                    altr = always_true(x[source->satnumber][i], false);
                    if(!altr) altr = always_true(sp[target->satnumber], false);
                    for(int j = 0; j < sinks_size; ++j){
                        if(altr) break;
                        if(sink_consistent(target, j) == true) altr = always_true(sy[label][i][j], true);
                    }
                    
                    if(altr == false){
                        for(int j = 0; j < sinks_size; ++j)
                            if(sink_consistent(target, j) == true) print_lit(sy[label][i][j], true);
                        print_lit(x[source->satnumber][i], false);
                        print_lit(sp[target->satnumber], false);
                        print_clause_end();
                        num += 1;
                    }
                }
                num += print_clause(false, sp[source->satnumber], true, sp[target->satnumber]);
            }
        }
    }
    for (state_set::iterator it = non_red_states.begin(); it != non_red_states.end(); ++it) {
        apta_node* source = *it;
        for (int label = 0; label < alphabet_size; ++label) {
            apta_node* target = source->get_child(label);
            if (target != 0) {
                for (int i = 0; i < dfa_size; ++i)
                    for (int j = 0; j < sinks_size; ++j)
                        if(sink_consistent(target, j) == false)
                            num += print_clause(false, sy[label][i][j], false, x[source->satnumber][i]);
                for (int i = 0; i < dfa_size; ++i){
                    altr = always_true(x[source->satnumber][i], false);
                    if(!altr) altr = always_true(sp[target->satnumber], false);
                    for(int j = 0; j < sinks_size; ++j){
                        if(altr) break;
                        if(sink_consistent(target, j) == true) altr = always_true(sy[label][i][j], true);
                    }
                    
                    if(altr == false){
                        for(int j = 0; j < sinks_size; ++j)
                            if(sink_consistent(target, j) == true) print_lit(sy[label][i][j], true);
                        print_lit(x[source->satnumber][i], false);
                        print_lit(sp[target->satnumber], false);
                        print_clause_end();
                        num += 1;
                    }
                }
                num += print_clause(false, sp[source->satnumber], true, sp[target->satnumber]);
            }
        }
    }
    return num;
}

/* output result to dot */
void print_dot_output(const char* dot_output){
    FILE* output = fopen(dot_output, "w");
    apta* aut = merger.aut;
    int v,i,a,j;
    
    fprintf(output,"digraph DFA {\n");
    fprintf(output,"\t\tI -> %i;\n", aut->root->find()->satnumber);
    
    set<int>::iterator it = trueliterals.begin();
    for(v = 0; v < num_states; ++v)
        for(i = 0; i < dfa_size; ++i)
            if(x[v][i] == *it) ++it;
    
    for(a = 0; a < alphabet_size; ++a){
        for(i = 0; i < dfa_size; ++i){
            for(j = 0; j < dfa_size; ++j) {
                if(y[a][i][j] == *it){
                    if(j != 0)
                        fprintf(output,"\t\t%i -> %i [label=\"%i\"];\n", i, j, a);
                    ++it;
                }
                if(y[a][i][j] == -1){
                    if(j != 0)
                        fprintf(output,"\t\t%i -> %i [label=\"%i\"];\n", i, j, a);
                }
            }
        }
    }
    
    for(a = 0; a < alphabet_size; ++a){
        for(i = 0; i < dfa_size; ++i){
            for(j = 0; j < sinks_size; ++j){
                if(sy[a][i][j] == *it){
                    ++it;
                    fprintf(output,"\t\t %i -> s%i [label=\"%i\"];\n", i, j, a);
                }
            }
        }
    }

    for(j = 0; j < sinks_size; ++j){
        fprintf(output,"\ts%i [shape=box];\n", j);
    }
    
    for(i = 0; i < num_states; ++i){
        if(sp[i] == *it){
            //cerr << "sp " << i << endl;
            ++it;
        }
    }
    
    for(i = 0; i < dfa_size; ++i){
        if(z[i] == *it){
            ++it;
            fprintf(output,"\t%i [shape=doublecircle];\n", i);
        } else if(z[i] == -1){
            fprintf(output,"\t%i [shape=doublecircle];\n", i);
        } else {
            fprintf(output,"\t%i [shape=Mcircle];\n", i);
        }
    }
    
    fprintf(output,"}\n");
    fclose(output);
}

/* output result to aut, for later processing in i.e. ensembles */
void print_aut_output(const char* aut_output){
    FILE* output = fopen(aut_output, "w");
    apta* aut = merger.aut;
    int v,i,a,j;
    
    fprintf(output,"%i %i\n", dfa_size, alphabet_size);
    fprintf(output,"%i\n", aut->root->find()->satnumber);
    
    set<int>::iterator it = trueliterals.begin();
    for(v = 0; v < num_states; ++v)
        for(i = 0; i < dfa_size; ++i)
            if(x[v][i] == *it) ++it;
    
    for(a = 0; a < alphabet_size; ++a){
        for(i = 0; i < dfa_size; ++i){
            for(j = 0; j < dfa_size; ++j) {
                if(y[a][i][j] == *it){
                    fprintf(output,"t %i %i %i\n", i, a, j);
                    ++it;
                }
                if(y[a][i][j] == -1){
                    fprintf(output,"t %i %i %i\n", i, a, j);
                }
            }
        }
    }
    
    for(a = 0; a < alphabet_size; ++a){
        for(i = 0; i < dfa_size; ++i){
            for(j = 0; j < sinks_size; ++j){
                if(sy[a][i][j] == *it){
                    ++it;
                    fprintf(output,"t %i %i %i;\n", i, a, dfa_size+j);
                }
            }
        }
    }
    
    for(i = 0; i < num_states; ++i){
        if(sp[i] == *it){
            //cerr << "sp " << i << endl;
            ++it;
        }
    }
    
    for(i = 0; i < dfa_size; ++i){
        if(z[i] == *it){
            ++it;
            fprintf(output,"a %i 1\n", i);
        } else if(z[i] == -1){
            fprintf(output,"a %i 1\n", i);
        } else {
            fprintf(output,"a %i 0\n", i);
        }
    }

    for(i = 0; i < sinks_size; ++i){
        fprintf(output,"a %i s%i\n", i+dfa_size, i);
    }
    fclose(output);
}

/* the main routine:
 * run greedy state merging runs
 * convert result to satisfiability formula
 * run sat solver
 * translate result to a DFA
 * print result
 * */
int dfasat(state_merger &m, const char* sat_program, const char* dot_output, const char* aut_output){
    int i,j,l,v,a,h,k;
    merger = m;
    apta* the_apta = merger.aut;
    
    merge_list merges = random_greedy_bounded_run(&merger);
    
    std::ostringstream oss2;
    oss2 << "pre_" << dot_output;
    FILE* output = fopen(oss2.str().c_str(), "w");
    merger.todot(output);
    fclose(output);
    
    non_red_states = merger.get_candidate_states();
    red_states     = merger.red_states;
    sink_states    = merger.get_sink_states();
    
    if(best_solution != -1 && merger.red_states.size() >= best_solution + EXTRA_STATES){
        cerr << "Greedy preprocessing resulted in too many red states." << endl;
        while(!merges.empty()){
            merge_pair performed_merge = merges.front();
            merges.pop_front();
            merger.undo_merge(performed_merge.first, performed_merge.second);
        }
        return -1;
    }
    
    dfa_size = merger.red_states.size() + OFFSET;
    sinks_size = 0;
    
    if(MERGE_SINKS_PRESOLVE) sinks_size = num_sink_types();
    else non_red_states.insert(sink_states.begin(), sink_states.end());
    num_states = red_states.size() + non_red_states.size();
    
    if(best_solution != -1) dfa_size = min(dfa_size, best_solution);
    new_states = dfa_size - merger.red_states.size();
    new_init = merger.red_states.size();
    
    // assign a unique number to every state
    i = 0;
    for(state_set::iterator it = red_states.begin(); it != red_states.end(); ++it){
        apta_node* node = *it;
        node->satnumber = i;
        node->colour = i;
        i++;
    }
    for(state_set::iterator it = non_red_states.begin(); it != non_red_states.end(); ++it){
        apta_node* node = *it;
        node->satnumber = i;
        i++;
    }
    
    clause_counter = 0;
    literal_counter = 1;
    
    cerr << "creating literals..." << endl;
    create_literals();
    
    cerr << "number of states: " << the_apta->get_states().size() << endl;
    cerr << "number of red states: " << red_states.size() << endl;
    cerr << "number of non_red states: " << non_red_states.size() << endl;
    cerr << "number of sink states: " << sink_states.size() << endl;
    cerr << "dfa size: " << dfa_size << endl;
    cerr << "sink types: " << sinks_size << endl;
    cerr << "new states: " << new_states << endl;
    cerr << "new init: " << new_init << endl;

    fix_red_values();
    erase_red_conflict_colours();
    set_symmetry();
    
    // renumber literals to account for eliminated ones
    reset_literals(false);

    computing_header = true;
    
    clause_counter = 0;
    clause_counter += print_colours();
    clause_counter += print_conflicts();
    clause_counter += print_accept();
    clause_counter += print_transitions();
    if(SYMMETRY_BREAKING){
        clause_counter += print_t_transitions();
        clause_counter += print_p_transitions();
        clause_counter += print_a_transitions();
        clause_counter += print_symmetry();
    }
    if(FORCING){
        clause_counter += print_forcing_transitions();
    }
    clause_counter += print_paths();
    clause_counter += print_sink_paths();
    
    cerr << "header: p cnf " << literal_counter - 1 << " " << clause_counter << endl;
    computing_header = false;
    
    int pipetosat[2];
    int pipefromsat[2];
    if (pipe(pipetosat) < 0 || pipe(pipefromsat) < 0){
        cerr << "Unable to create pipe for SAT solver: " << strerror(errno) << endl;
        exit(1);
    }
    pid_t pid = fork();
    if (pid == 0){
        close(pipetosat[1]);
        dup2(pipetosat[0], STDIN_FILENO);
        close(pipetosat[0]);
        
        close(pipefromsat[0]);
        dup2(pipefromsat[1], STDOUT_FILENO);
        close(pipefromsat[1]);
        
        cerr << "starting SAT solver " << sat_program << endl;
        char* copy_sat = strdup(sat_program);
        char* pch = strtok (copy_sat," ");
        vector<char*> args;
        while (pch != NULL){
            args.push_back(strdup(pch));
            pch = strtok (NULL," ");
        }
        free(copy_sat);
        free(pch);
        args.push_back((char*)NULL);
        execvp(args[0], &args[0]);
        cerr << "finished SAT solver" << endl;
        for(int argi = 0; argi < args.size(); ++argi) free(args[argi]);
        int* status;
        WIFEXITED(status);
        wait(status);
    }
    else
    {
        close(pipetosat[0]);
        close(pipefromsat[1]);
        
        sat_stream = (FILE*) fdopen(pipetosat[1], "w");
        //sat_stream = (FILE*) fopen("test.out", "w");
        if (sat_stream == 0){
            cerr << "Cannot open pipe to SAT solver: " << strerror(errno) << endl;
            exit(1);
        }
        fprintf(sat_stream, "p cnf %i %i\n", literal_counter - 1, clause_counter);

        print_colours();
        print_conflicts();
        print_accept();
        print_transitions();
        if(SYMMETRY_BREAKING){
            print_symmetry();
            print_t_transitions();
            print_p_transitions();
            print_a_transitions();
        }
        if(FORCING){
            print_forcing_transitions();
        }
        print_paths();
        print_sink_paths();
        
        fclose(sat_stream);
        
        cerr << "sent problem to SAT solver" << endl;
        
        time_t begin_time = time(nullptr);
        
        trueliterals = set<int>();
        
        char line[500];
        sat_stream = fdopen ( pipefromsat[0], "r" );
        
        bool improved = false;
        while(fgets ( line, sizeof line, sat_stream ) != NULL){
            char* pch = strtok (line," ");
            if(strcmp(pch,"s") == 0){
                pch = strtok (NULL, " ");
                cerr << pch << endl;
                if(strcmp(pch,"SATISFIABLE\n")==0){
                    cerr << "new solution, size = " << dfa_size << endl;
                    if(best_solution ==-1 || best_solution > dfa_size){
                        cerr << "new best solution, size = " << dfa_size << endl;
                        best_solution = dfa_size;
                        improved = true;
                    }
                }
            }
            if(strcmp(pch,"v") == 0){
                pch = strtok (NULL, " ");
                while(pch != NULL){
                    int val = atoi(pch);
                    if(val > 0) trueliterals.insert(val);
                    pch = strtok (NULL, " ");
                }
            }
        }
        fclose(sat_stream);
        
        cerr << "solving took " << (time(nullptr) - begin_time) << " seconds" << endl;

        if(improved){
            print_dot_output(dot_output);
            print_aut_output(aut_output);
        }
        
        while(!merges.empty()){
            merge_pair performed_merge = merges.front();
            merges.pop_front();
            merger.undo_merge(performed_merge.first, performed_merge.second);
        }
    }
    
    delete_literals();
    return best_solution;
};
