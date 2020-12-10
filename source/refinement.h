/*
 *  RTI (real-time inference)
 *  Searcher.cpp, the header file for the search routines
 *  Currently, only a simple greedy (best-first) routine is implemented, search routines will be added later.
 *
 *  A refinement is either a point (merge), split, or color (adding of a new state) in the current real-time automaton, as in:
 *  Sicco Verwer and Mathijs de Weerdt and Cees Witteveen (2007),
 *  An algorithm for learning real-time automata,
 *  In Maarten van Someren and Sophia Katrenko and Pieter Adriaans (Eds.),
 *  Proceedings of the Sixteenth Annual Machine Learning Conference of Belgium and the Netherlands (Benelearn),
 *  pp. 128-135.
 *  
 *  Copyright 2009 - Sicco Verwer, jan-2009
 *  This program is released under the GNU General Public License
 *  Info online: http://www.gnu.org/licenses/quick-guide-gplv3.html
 *  Or in the file: licence.txt
 *  For information/questions contact: siccoverwer@gmail.com
 *
 *  I will try to keep this software updated and will also try to add new statistics or search methods.
 *  Also I will add comments to the source in the near future.
 *
 *  Feel free to adapt the code to your needs, please inform me of (potential) improvements.
 */

#ifndef _REFINEMENT_H_
#define _REFINEMENT_H_

using namespace std;

#include <list>
#include <queue>
#include <map>
#include <set>

class refinement;
class merge_refinement;
class split_refinement;
class extend_refinement;
struct score_compare;

typedef list<refinement*> refinement_list;
typedef set<refinement*, score_compare > refinement_set;

#include "apta.h"

/**
 * @brief Base class for refinements. Specialized to either
 * a point (merge), split, or color (adding of a new state).
 *
 */
class refinement{
public:
    double score;
	tail* red;

	apta_node* tempnode;
    apta_node* tempblue;

	virtual void print() const;
	virtual void print_short() const;
	virtual void doref(state_merger* m);
	virtual void undo(state_merger* m);
    virtual bool testref(state_merger* m);

    virtual void erase();

    virtual void print_json(iostream &output) const;

    static void print_refinement_list_json(iostream &output, refinement_list *list);
};

/**
 * @brief A merge-refinement assigns a score to the merge of
 * a left and right state. 
 *
 */
class merge_refinement : public refinement {
public:
	tail* blue;

	merge_refinement(state_merger* m, double s, apta_node* l, apta_node* r);
    void initialize(state_merger* m, double s, apta_node* l, apta_node* r);

	virtual inline void print() const;
	virtual inline void print_short() const;
	virtual inline void doref(state_merger* m);
	virtual inline void undo(state_merger* m);
    virtual inline bool testref(state_merger* m);

    virtual inline void erase();

    virtual void print_json(iostream &output) const;
};

 /**
 * @brief A extend-refinement makes a blue state red. The
 * score is the size (frequency) of the state in the APTA.
 *
 */
class extend_refinement : public refinement {
public:
    int size;

	extend_refinement(state_merger* m, apta_node* r);
    void initialize(state_merger* m, apta_node* r);

	virtual inline void print() const;
	virtual inline void print_short() const;
	virtual inline void doref(state_merger* m);
	virtual inline void undo(state_merger* m);
    virtual inline bool testref(state_merger* m);

    virtual inline void erase();

    virtual void print_json(iostream &output) const;
};

class split_refinement : public refinement {
public:
    tail* split_point;
    int attribute;

	split_refinement(state_merger* m, double s, apta_node* l, tail* t, int a);
	void initialize(state_merger* m, double s, apta_node* l, tail* t, int a);

	virtual inline void print() const;
	virtual inline void print_short() const;
	virtual inline void doref(state_merger* m);
	virtual inline void undo(state_merger* m);
    virtual inline bool testref(state_merger* m);

    virtual inline void erase();

    virtual void print_json(iostream &output) const;
};

 /**
 * @brief Compare function for refinements, based on scores.
 *
 */
struct score_compare {
    merge_refinement* mref;
    split_refinement* sref;

    inline bool operator()(refinement* r1, refinement* r2) const {
        if(typeid(r1) == typeid(sref) && typeid(r2) != typeid(sref)) return 1;
        if(typeid(r1) != typeid(sref) && typeid(r2) == typeid(sref)) return 0;
        if(typeid(r1) == typeid(mref) && typeid(r2) != typeid(mref)) return 1;
        if(typeid(r1) != typeid(mref) && typeid(r2) == typeid(mref)) return 0;
            if(r1->score == r2->score){
                return r1 > r2;
                /*
                if(r1->right->size == r2->right->size){
                    return r1 > r2;
                } else {
                    return r1->right->size > r2->right->size;
                }
                 */
            }
            return r1->score > r2->score;
    }
};

#endif /* _REFINEMENT_H_ */
