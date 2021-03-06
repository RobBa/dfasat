#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include <string>
#include <vector>

using namespace std;

extern int alphabet_size;
extern bool MERGE_SINKS;
extern int STATE_COUNT;
extern int SYMBOL_COUNT;
extern int SINK_COUNT;
extern float CORRECTION;
extern float CHECK_PARAMETER;
extern bool USE_SINKS;
extern bool USE_LOWER_BOUND;
extern float LOWER_BOUND;
extern int alphabet_size;
extern int GREEDY_METHOD;
extern int APTA_BOUND;
extern int CLIQUE_BOUND;
extern bool EXTEND_ANY_RED;
extern bool MERGE_SINKS_PRESOLVE;
extern int OFFSET;
extern int EXTRA_STATES;
extern bool TARGET_REJECTING;
extern bool SYMMETRY_BREAKING;
extern bool FORCING;
extern bool MERGE_SINKS_DSOLVE;
extern string eval_string;
extern string OUTPUT;
extern bool MERGE_MOST_VISITED;
extern bool MERGE_BLUE_BLUE;
extern bool RED_FIXED;
extern bool MERGE_WHEN_TESTING;
extern bool DEPTH_FIRST;
extern int RANGE;
extern int STORE_MERGES;
extern int STORE_MERGES_KEEP_CONFLICT;
extern int STORE_MERGES_SIZE_THRESHOLD;
extern double STORE_MERGES_RATIO_THRESHOLD;
extern string COMMAND;
extern int STREAM_COUNT;
extern bool EXCEPTION4OVERLAP;
extern string EVALPAR;
extern bool FINAL_PROBABILITIES;
extern int MERGE_LOCAL;
extern int MERGE_LOCAL_COLLECTOR_COUNT;
extern int KTAIL;
extern int KSTATE;
extern bool MARKOVIAN_MODEL;
extern bool MERGE_SINKS_WITH_CORE;

extern bool PERTYPE_DISTRIBUTIONS;
extern bool TYPE_DISTRIBUTIONS;
extern bool TYPE_CONSISTENT;

extern bool MERGE_ROOT;
extern bool PRINT_WHITE;
extern bool PRINT_BLUE;

extern string DFA_FILE;

class parameters{
public:
    string command;
    string dfa_file; // TODO: name misleading. Shall we replace it by something better?
    vector<string> dfa_data;
    string dot_file;
    string sat_program;
    string hName;
    string hData;
    string output;
    int runs;
    int sinkson;
    int seed;
    int sataptabound;
    int satdfabound;
    float lower_bound;
    int satextra;
    int mergesinks;
    int satmergesinks;
    int method;
    int extend;
    int heuristic;
    int symbol_count;
    int state_count;
    int sink_count;
    float correction;
    float extrapar;
    int satplus;
    int satfinalred;
    int symmetry;
    int forcing;
    int blueblue;
    int finalred;
    int largestblue;
    int testmerge;
    int shallowfirst;
    string mode;
    string evalpar;
    int batchsize;
    float delta;
    float epsilon;
    int debugging;
    parameters();
};

void init_with_params(parameters*);
#endif
