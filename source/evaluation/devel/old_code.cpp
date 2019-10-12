//
//  old_code.cpp
//  
//
//  Created by sicco on 05/04/16.
//
//

            else {
                if(positive){
                    node->num_pos[c] = node->pos(c) + 1;
                    node->pos_paths()++;
                } else {
                    node->num_neg[c] = node->neg(c) + 1;
                    node->neg_paths()++;
                }
            }
            //node->input_output
            if(occ >= 0)
                node->occs.push_front(occ);
                //node->child(c)->occs.push_front(occ);
        if(positive) node->pos_final()++;
        else node->neg_final()++;

    right->merge_point = left->conflicts.end();
    --(right->merge_point);
    left->conflicts.splice(left->conflicts.end(), right->conflicts);
    ++(right->merge_point);

    right->occ_merge_point = left->occs.end();
    --(right->occ_merge_point);
    left->occs.splice(left->occs.end(), right->occs);
    ++(right->occ_merge_point);

    for(num_map::iterator it = right->num_pos.begin();it != right->num_pos.end(); ++it){
        left->num_pos[(*it).first] = left->pos((*it).first) + (*it).second;
    }
    for(num_map::iterator it = right->num_neg.begin();it != right->num_neg.end(); ++it){
        left->num_neg[(*it).first] = left->neg((*it).first) + (*it).second;
    }

    left->pos_final() -= right->pos_final();
    left->neg_final() -= right->neg_final();
    left->pos_paths() -= right->pos_paths();
    left->neg_paths() -= right->neg_paths();

    //left->depth = left->old_depth;

    for(num_map::iterator it = right->num_pos.begin();it != right->num_pos.end(); ++it){
        left->num_pos[(*it).first] = left->pos((*it).first) - (*it).second;
    }
    for(num_map::iterator it = right->num_neg.begin();it != right->num_neg.end(); ++it){
        left->num_neg[(*it).first] = left->neg((*it).first) - (*it).second;
    }

