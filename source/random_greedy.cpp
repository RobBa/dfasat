//#include <malloc.h>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <cstdlib>

#include "refinement.h"
#include "random_greedy.h"
#include "parameters.h"

// global variable
bool print_all_models = true;

// function in object oriented design
refinement_list* random_greedy_bounded_run(state_merger* merger){
    cerr << "starting greedy merging" << endl;
    int num = 1;
    refinement_list* all_refs = new refinement_list();
    merger->eval->initialize(merger);
    while( true ){
        merger->reset();
        while( true ){
            cout << " ";
            
            if(print_all_models){
            // output: dot legacy format
            merger->todot();
            std::ostringstream oss2;
            oss2 << "pre_" << num << ".dot";
            ofstream output1(oss2.str().c_str());
            output1 << merger->dot_output;
            output1.close();
 
            // output: new json format
            merger->tojson();
            std::ostringstream oss3;
            oss3 << "pre_" << num << ".json";
            ofstream output2(oss3.str().c_str());
            output2 << merger->json_output;
            output2.close();
            }

            // if(EXTEND_ANY_RED) while(merger->extend_red() != 0) cerr << "+ ";
            // leak here, too
            //merge_map* possible_merges = merger->get_possible_merges();
            
            refinement_set* refs = merger->get_possible_refinements();

            if(print_all_models) {
                merger->todot();
                std::ostringstream oss2;
                oss2 << "pre_refs" << num << ".dot";
                ofstream output1(oss2.str().c_str());
                output1 << merger->dot_output;
                output1.close();
            }

            /* cerr << endl;

            for(refinement_set::iterator it = refs->begin(); it != refs->end(); ++it){
                (*it)->print_short();
                cerr << " , ";
            }
            cerr << endl; */

            if(refs->empty()){
                cerr << "no more possible merges" << endl;
                delete refs;
                break;
            }
            if(merger->red_states.size() > CLIQUE_BOUND){
                cerr << "too many red states" << endl;
                for(refinement_set::iterator it = refs->begin(); it != refs->end(); ++it){
                    (*it)->erase();
                }
                delete refs;
                break;
            }
            // FIXME
            /*if(merger->get_final_apta_size() <= APTA_BOUND){
               cerr << "APTA too small" << endl;
               break;
            }*/

            refinement* best_ref = *refs->begin();
            /* if(GREEDY_METHOD == RANDOMG){
                merge_map randomized_merges;
                for(merge_map::reverse_iterator it = possible_merges->rbegin(); it != possible_merges->rend(); it++){
                    //if(it->first < LOWER_BOUND) break;
                    randomized_merges.insert(pair<int, merge_pair>(it->first * (rand() / (double)RAND_MAX), it->second));
                }
                top_score = (*randomized_merges.rbegin()).first;
                top_pair = (*randomized_merges.rbegin()).second;
            }*/

            best_ref->print_short();
            cout << " ";
            std::cout.flush();
            best_ref->doref(merger);
            all_refs->push_back(best_ref);
            
            if(print_all_models) {
                merger->todot();
                std::ostringstream oss3;
                oss3 << "post_refs" << num << ".dot";
                ofstream output2(oss3.str().c_str());
                output2 << merger->dot_output;
                output2.close();
            }

            for(refinement_set::iterator it = refs->begin(); it != refs->end(); ++it){
                if(*it != best_ref) (*it)->erase();
            }
            delete refs;
            num = num + 1;
        }
        cout << endl;
        int size =  merger->get_final_apta_size();
        int red_size = merger->red_states.size();
        cout << endl << "found intermediate solution with " << size << " and " << red_size << " red states" << endl;

        std::ostringstream ref_os;
        ref_os << "refinements.json";
        ofstream ref_output(ref_os.str().c_str());
        stringstream output_buf;
        refinement::print_refinement_list_json(output_buf, all_refs);
        ref_output << output_buf.str();
        ref_output.close();

        return all_refs;
    }
    return all_refs;
};


