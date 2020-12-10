/*
 */
#include "utility/loguru.hpp"

#include "inputdata.h"


json inputdata::all_data;
vector<string> inputdata::alphabet;
map<string, int> inputdata::r_alphabet;
vector<string> inputdata::types;
map<string, int> inputdata::r_types;
int inputdata::num_symbol_attributes;
int inputdata::num_trace_attributes;
int inputdata::num_attributes;
vector<attribute> inputdata::trace_attributes;
vector<attribute> inputdata::symbol_attributes;
int inputdata::num_sequences;
int inputdata::alphabet_size;
int inputdata::node_number;

void tail::split(tail* t){
    t->split_from = this;
    t->future_tail = future_tail;
    t->past_tail = past_tail;
    //if(past_tail != 0) past_tail->future_tail = t;
    //if(future_tail != 0) future_tail->past_tail = t;
    split_to = t;
};

void tail::undo_split(){
    //delete split_to;
    split_to = 0;
    //if(past_tail != 0) past_tail->future_tail = this;
    //if(future_tail != 0) future_tail->past_tail = this;
};

tail* tail::next(){
    if(split_to == 0) return next_in_list;
    if(next_in_list == 0) return 0;
    return next_in_list->next();
};

tail* tail::splitted(){
    if(split_to == 0) return this;
    return split_to->splitted();
};

tail* tail::future(){
    if(future_tail == 0) return 0;
    if(split_to == 0) return future_tail->splitted();
    return split_to->future();
};

tail* tail::past(){
    if(past_tail == 0) return 0;
    if(split_to == 0) return past_tail->splitted();
    return split_to->past();
};

tail::tail(tail* ot){
    td = ot->td;

    past_tail = 0;
    future_tail = 0;
    next_in_list = 0;
    split_from = 0;
    split_to = 0;
};

void tail::initialize(tail* ot){
    td = ot->td;

    past_tail = 0;
    future_tail = 0;
    next_in_list = 0;
    split_from = 0;
    split_to = 0;
};

tail::tail(int seq, int i, tail* pt){
    td = new tail_data();

    td->sequence = seq;
    td->index = i;

    td->type = inputdata::get_type(seq);
    td->length = inputdata::get_length(seq);
    td->symbol = inputdata::get_symbol(seq,i);
    td->attr = new int[inputdata::num_symbol_attributes];
    for(int k = 0; k < inputdata::num_symbol_attributes; k++){
        td->attr[k] = inputdata::get_value(seq, i, k);
    }
    td->data = inputdata::get_data(seq, i);

    past_tail = pt;
    if(past_tail != 0){
        past_tail->future_tail = this;
        td->trace_attr = pt->td->trace_attr;
    } else {
        td->trace_attr = new int[inputdata::num_trace_attributes];
        for(int k = 0; k < inputdata::num_trace_attributes; k++){
            td->trace_attr[k] = inputdata::get_trace_value(seq, i, k);
        }
    }
    
    future_tail = 0;
    next_in_list = 0;
    split_from = 0;
    split_to = 0;
};

attribute::attribute(string input){
    discrete = false;
    splittable = false;
    distributionable = false;
    target = false;

    if(input.find("d") != std::string::npos) discrete = true;
    if(input.find("s") != std::string::npos) splittable = true;
    if(input.find("f") != std::string::npos) distributionable = true;
    if(input.find("t") != std::string::npos) target = true;
}

inputdata::inputdata() {
    node_number = 0;
}

void inputdata::read_json_file(istream &input_stream){

    all_data = json::parse(input_stream);

    //for each json line
    for (int i = 0; i < all_data.size(); ++i) {
        //get the objects of each sequence
        for (auto& el : all_data[i].items()){
            //the alphabet needs to be identified, thus keep track of unique symbols
            if (el.key() == "S"){
                for (int j = 0; j < el.value().size(); ++j) {
                    int temp_symbol = el.value()[j];
                    string symbol = std::to_string(temp_symbol);
                    //if symbol not in alphabet, add it
                    if (std::find(alphabet.begin(), alphabet.end(), symbol)==alphabet.end()){
                        alphabet.push_back(symbol);
                    }
                }
            }
        }
        if (alphabet.size()==0){
            LOG_S(ERROR) << "Json wrongly formatted / no symbol alphabet. Aborting...";
            exit(-1);
        }
    }

};

void inputdata::abbadingo_init(istream &input_stream){
    input_stream >> inputdata::num_sequences;

    string tuple;
    input_stream >> tuple;

    std::stringstream lineStream;
    lineStream.str(tuple);

    string alph;
    std::getline(lineStream,alph,':');
    inputdata::alphabet_size = stoi(alph);

    string trace_attr;
    std::getline(lineStream,trace_attr, ':');
    string symbol_attr;
    std::getline(lineStream,symbol_attr);
    if(symbol_attr.empty()){
        symbol_attr = trace_attr;
        trace_attr = "";
    }

    if(!trace_attr.empty()){
        lineStream.str(trace_attr);
        string attr;
        std::getline(lineStream,attr, ',');
        while(!attr.empty()){
            inputdata::trace_attributes.push_back(attribute(attr));
            std::getline(lineStream,attr, ',');
        }
    }

    if(!symbol_attr.empty()){
        lineStream.str(symbol_attr);
        string attr;
        std::getline(lineStream,attr, ',');
        while(!attr.empty()){
            inputdata::symbol_attributes.push_back(attribute(attr));
            std::getline(lineStream,attr, ',');
        }
    }
};

void inputdata::read_abbadingo_file(istream &input_stream){
    int num_sequences;
    input_stream >> num_sequences;
    
    string tuple;
    input_stream >> tuple;

    std::stringstream lineStream;
    lineStream.str(tuple);

    string alph;
    std::getline(lineStream,alph,':');
    inputdata::alphabet_size = stoi(alph);

    string trace_attr;
    std::getline(lineStream,trace_attr, ':');
    string symbol_attr;
    std::getline(lineStream,symbol_attr);
    if(symbol_attr.empty()){
        symbol_attr = trace_attr;
        trace_attr = "";
    }
    lineStream.clear();

    int num = 0;
    if(!trace_attr.empty()){
        lineStream.str(trace_attr);
        string attr;
        std::getline(lineStream,attr, ',');
        while(!std::getline(lineStream,attr, ',').eof()){
            num += 1;
            inputdata::trace_attributes.push_back(attribute(attr));
        }
        num += 1;
        inputdata::trace_attributes.push_back(attribute(attr));
        lineStream.clear();
    }
    inputdata::num_trace_attributes = num;

    num = 0;
    if(!symbol_attr.empty()){
        lineStream.str(symbol_attr);
        string attr;
        std::getline(lineStream,attr, ',');
        while(!std::getline(lineStream,attr, ',').eof()){
            num += 1;
            inputdata::symbol_attributes.push_back(attribute(attr));
        }
        num += 1;
        inputdata::symbol_attributes.push_back(attribute(attr));
        lineStream.clear();
    }
    inputdata::num_symbol_attributes = num;

    inputdata::num_attributes = inputdata::num_trace_attributes + inputdata::num_symbol_attributes;

    lineStream.clear();

    for(int line = 0; line < num_sequences; ++line){
        read_abbadingo_sequence(input_stream);
    }
};

void inputdata::read_abbadingo_sequence(istream &input_stream){
    string temp, temp_symbol, data, type_string, type_attr, symbol_string, symbol_attr, val;
    int length;
    std::stringstream l1, l2, l3;

    data = "";

    json sequence;

    input_stream >> temp;
    l1.str(temp);
    std::getline(l1,type_string,':');
    std::getline(l1,type_attr);
    l1.clear();

    if(r_types.find(type_string) == r_types.end()){
        r_types[type_string] = types.size();
        types.push_back(type_string);
    }

    input_stream >> length;

    vector<int> symbols(length);
    vector< float > t_values(num_trace_attributes);
    vector< vector<float> > s_values(num_symbol_attributes, vector<float>(length));
    vector< string > datas(length);

    if(num_trace_attributes > 0){
        l2.str(type_attr);
        for(int i = 0; i < num_trace_attributes; ++i){
            if(i < num_trace_attributes - 1) std::getline(l2,val,',');
            else std::getline(l2,val);
            t_values[i] = trace_attributes[i].get_value(val);
        }
        l2.clear();
    }

    sequence["T"] = r_types[type_string];
    sequence["L"] = length;

    for(int index = 0; index < length; ++index){
        input_stream >> temp;

        l1.str(temp);
        std::getline(l1,temp_symbol,'/');
        std::getline(l1,data);
        l1.clear();

        l2.str(temp_symbol);
        std::getline(l2,symbol_string,':');
        std::getline(l2,symbol_attr);
        l2.clear();

        //cerr << temp << " - " << temp_symbol << endl;

        if(r_alphabet.find(symbol_string) == r_alphabet.end()){
            r_alphabet[symbol_string] = alphabet.size();
            alphabet.push_back(symbol_string);
        }
        
        symbols[index] = r_alphabet[symbol_string];
        datas[index] = data;

        if(num_symbol_attributes > 0){
            l3.str(symbol_attr);
            for(int i = 0; i < num_symbol_attributes; ++i){
                if(i < num_symbol_attributes - 1) std::getline(l3,val,',');
                else std::getline(l3,val);
                s_values[i][index] = symbol_attributes[i].get_value(val);
            }
            l3.clear();
        }

        /*string val;
        if(num_attributes != 0){
            std::stringstream l3;
            l3.str(vals);
            for(int i = 0; i < num_attributes-1; ++i){
                std::getline(l3,val,',');
                //cerr << val << endl;
                values[i][index] = stof(val);
            }
            std::getline(l3,val);
            //cerr << val << endl;
            values[num_attributes-1][index] = stof(val);
            //cerr << val << endl;
        }
         */
    }

    //for(int i = 0; i < length; ++i)
    //    cerr << symbols[i] << endl;

    sequence["S"] = symbols;
    sequence["VT"] = t_values;
    for(int i = 0; i < num_symbol_attributes; ++i){
        sequence["V" + to_string(i)] = s_values[i];
    }
    sequence["D"] = datas;
    
    all_data.push_back(sequence);
};

void inputdata::add_data_to_apta(apta* the_apta){
    for(int i = 0; i < all_data.size(); ++i){
        add_sequence_to_apta(the_apta, i);
    } 
};

void inputdata::add_sequence_to_apta(apta* the_apta, int seq_nr){
    json sequence = all_data[seq_nr];
    
    int depth = 0;
    apta_node* node = the_apta->root;
    node->label = -1;
    tail* ot = 0;

    for(int index = 0; index < sequence["L"]; index++){
        depth++;
        tail* nt = new tail(seq_nr, index, ot);
        int symbol = sequence["S"][index];
        if(node->child(symbol) == 0){
            apta_node* next_node = new apta_node(the_apta);
            node->set_child(symbol, next_node);
            next_node->source = node;
            next_node->label  = symbol;
            next_node->depth  = depth;
            next_node->number = ++(this->node_number);
        }
        node->size = node->size + 1;
        node->add_tail(nt);
        node->data->add_tail(nt);
        //node->data->read_from(seq_nr, index);
        //node->data->read_from(nt);
        apta_node* node2 = node->child(symbol);
        //node2->data->read_to(seq_nr, index);
        //node2->data->read_to(nt);
        node = node2;
        ot = nt;
    }
    
    tail* nt = new tail(seq_nr, -1, ot);
    //node->data->read_to(nt);
    node->type = sequence["T"];
    node->size = node->size + 1;
    node->final = node->final + 1;
    node->add_tail(nt);
    node->data->add_tail(nt);
};

const string inputdata::to_json_str() const{
    ostringstream ostr;
    ostr << all_data;
    return ostr.str();
};

const string tail::to_string(){
    ostringstream ostr;
    tail* t = this;
    while(t->past() != 0) t = t->past_tail;

    ostr << "[ ";
    while(t != this->future_tail){
        ostr << "\"" << inputdata::alphabet[inputdata::get_symbol(t)];
        if(inputdata::num_symbol_attributes > 0){
            ostr << ":";
            for(int i = 0; i < inputdata::num_symbol_attributes; i++){
                ostr << inputdata::get_value(t, i);
                if(i + 1 < inputdata::num_symbol_attributes)
                    ostr << ",";
            }
        }
        if(inputdata::get_data(t) != "") {
            ostr << "/" << inputdata::get_data(t);
        }
        t = t->future_tail;
        if(t != this->future_tail){
            ostr << "\" , ";
        }
    }
    ostr << "\" ]";
    return ostr.str();
};

tail::~tail(){
    if(split_from == 0){
        delete td->attr;
        delete td;
    }
};

/*
const string inputdata::to_abbadingo_str() const{
    ostringstream ostr;
	ostr << num_sequences  << " " << alph_size << ":" << num_attributes << "\n";
	for(int line = 0; line < num_sequences; ++line) {
	    sequence* seq = sequences[line];
	    ostr << seq->type << " " << seq->length << " ";
	    for(int index = 0; index < seq->length; ++index){
            ostr << inputdata::alphabet[seq->symbols[index]];
            if(inputdata::num_symbol_attributes != 0){
                ostr << inputdata::alphabet[seq->symbols[index]] << ":";
                for(int val = 0; val < inputdata::num_symbol_attributes-1; ++val){
                    ostr << seq->values[index][val] << ",";
                }
                ostr << seq->values[index][inputdata::num_symbol_attributes-1] << " ";
            }
            if(!seq->data[index].empty()){
                ostr << "/" << seq->data;
            }
            ostr << " ";
        }
        ostr << "\n";
    }
    return ostr.str();
};

*/

