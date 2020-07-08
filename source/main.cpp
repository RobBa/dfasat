/// @file
/// \brief Command line parameter processing and entry point.

#include <cstdlib>
#include "random_greedy.h"
#include "state_merger.h"
#include "evaluate.h"
#include "dfasat.h"
#include <sstream>
#include <iostream>
#include "evaluation_factory.h"
#include <string>
#include "searcher.h"
#include "stream.h"
#include <vector>
#include "interactive.h"
#include "common.h"
#include "inputdata.h"


#include "utility/loguru.hpp"
#include "utility/CLI11.hpp"

#include "parameters.h"

// this file is generaterated by collector.sh
// during make, so contributors don't have to
// touch anything than their own files
#include "evaluators.h"

using namespace std;

string COMMAND;
bool debugging_enabled = false;

/*
 * Input parameters, see 'man popt'
 */

// map_type* BaseFactory::map = NULL;


/*! \brief Set global vars to values in parameter object.
 *         
 * @param[in] param pointer to initialized parameter object.
 */
void init_with_params(parameters* param) {

    srand(param->seed);

    APTA_BOUND = param->sataptabound;
    CLIQUE_BOUND = param->satdfabound;

    STATE_COUNT = param->state_count;
    SINK_COUNT = param->sink_count;
    SYMBOL_COUNT = param->symbol_count;
    CHECK_PARAMETER = param->extrapar;
    CORRECTION = param->correction;

    LOWER_BOUND = param->lower_bound;
    USE_LOWER_BOUND = true;
    OFFSET = param->satextra;
    USE_SINKS = param->sinkson;
    MERGE_SINKS_PRESOLVE = param->mergesinks;
    MERGE_SINKS_DSOLVE = param->satmergesinks;
    EXTEND_ANY_RED = param->extend;

    SYMMETRY_BREAKING = param->symmetry;
    FORCING = param->forcing;

    EXTRA_STATES = param->satplus;
    TARGET_REJECTING = param->satfinalred;

    /* new since master-dev merge */
    MERGE_MOST_VISITED = param->largestblue;
    MERGE_BLUE_BLUE = param->blueblue;
    RED_FIXED = param->finalred;
    MERGE_WHEN_TESTING = !param->testmerge;
    DEPTH_FIRST = param->shallowfirst;

    COMMAND = param->command;

    if(param->debugging > 0) debugging_enabled = true;

    if(param->method == 1) GREEDY_METHOD = RANDOMG;
    if(param->method == 2) GREEDY_METHOD = NORMALG;

    OUTPUT = param->output;
    OUTPUT = "both";

    eval_string = param->hData;

    try {
       (DerivedRegister<evaluation_function>::getMap())->at(param->hName);
       std::cout << "valid: " << param->hName << std::endl;
       if(param->hName == "overlap4logs") {
           EXCEPTION4OVERLAP = true;
           LOG_S(ERROR) << "overlap4logs backcompatilibity";
       }

    } catch(const std::out_of_range& oor ) {
        std::cout << "invalid: " << param->hName << std::endl;
        LOG_S(ERROR) << "invalid: " << param->hName;
    }
}


void run(parameters* param) {

    state_merger merger;

    init_with_params(param);

    evaluation_function *eval;

    for(auto myit = DerivedRegister<evaluation_function>::getMap()->begin(); myit != DerivedRegister<evaluation_function>::getMap()->end() && debugging_enabled; myit++   ) {
       cout << myit->first << " " << myit->second << endl;
    }

    try {
       eval = (DerivedRegister<evaluation_function>::getMap())->at(param->hName)();
       std::cout << "Using heuristic " << param->hName << std::endl;
       LOG_S(INFO) <<  "Using heuristic " << param->hName;
    } catch(const std::out_of_range& oor ) {
        LOG_S(WARNING) << "No named heuristic found, defaulting back on -h flag";
        std::cerr << "No named heuristic found, defaulting back on -h flag" << std::endl;
    }

    ifstream input_stream(param->dfa_file);

    if(!input_stream) {
        LOG_S(ERROR) << "Input file not found, aborting";
        std::cerr << "Input file not found, aborting" << std::endl;
        exit(-1);
    }

    inputdata id;
    string file_name = param->dfa_file;
    // based on the file extension run the corresponding parser
    if (file_name.substr(file_name.find_last_of(".") + 1) == "json") {
        LOG_S(INFO) << "Parse json input file";
        id.read_json_file(input_stream);
    } else {
        LOG_S(INFO) << "Parse dot input file";
        id.read_abbadingo_file(input_stream);
    }
    input_stream.close();

    apta* the_apta = new apta();

    eval->set_params(param->evalpar);

    cout << "Creating apta " <<  "using evaluation class " << eval_string << endl;

    
    eval->initialize_globals();

    cerr << "done parsing" << endl;
    
    merger = state_merger(eval,the_apta);
    the_apta->context = &merger;

    //cerr << id.to_json_str();
    
    if(param->mode == "batch") {
        cout << "batch mode selected" << endl;

        //merger.read_apta(input_stream); TODO: error checking and logging
        LOG_S(INFO) << "Start reading data in batch mode";
        id.add_data_to_apta(the_apta);
        the_apta-> alp = inputdata::alphabet;
        LOG_S(INFO) << "Finished reading data in batch mode";

        cout << "reading data finished, processing:" << endl;
        // run the state merger
        int solution = -1;

        std::ostringstream oss3;
        oss3 << param->dot_file  << "init_dfa";

        if(OUTPUT == "dot" || OUTPUT == "both") {
            oss3 << ".dot";
            FILE* output = fopen(oss3.str().c_str(), "w");
            merger.todot();
            merger.print_dot(output);
            fclose(output);
        }

        if(OUTPUT == "json" || OUTPUT == "both") {
            oss3<< ".json";
            FILE* output = fopen(oss3.str().c_str(), "w");
            merger.tojson();
            merger.print_json(output);
            fclose(output);
        }

        //bestfirst(&merger);

       for(int i = 0; i < param->runs; ++i){
          std::ostringstream oss;
          oss << param->dot_file << "dfa" << (i+1) << ".aut";
          std::ostringstream oss2;
          oss2 << param->dot_file << "dfa" << (i+1);

          LOG_S(INFO) << "Starting state-merging in batch mode";
          cout << "dfasat running";

          solution = dfasat(merger, param->sat_program, oss2.str().c_str(), oss.str().c_str());
          if(solution != -1)
             CLIQUE_BOUND = min(CLIQUE_BOUND, solution - OFFSET + EXTRA_STATES);
         }

    } else if(param->mode == "stream") {
       cout << "stream mode selected" << endl;
       LOG_S(INFO) << "Stream mode selected, starting run";

       stream_mode(&merger, param, input_stream);
    } else if(param->mode == "inter") {

       LOG_S(INFO) << "Reading data in interactive mode";

       // process data
       id.add_data_to_apta(the_apta);
       the_apta-> alp = id.alphabet;
       //merger.read_apta(input_stream);
       //input_stream.close();

       LOG_S(INFO) << "Finished reading input data, starting state-merging";
       cout << "Reading data finished, processing:" << endl;
 
       // run interactive loop
       interactive(&merger, param); 
    } else {
       LOG_S(ERROR) << "Unkown mode of operation selected, please chose batch, stream, or inter. Provided was " << param->mode;
       cerr << R"(unknown mode of operation selected, valid options are "batch", "stream", and "inter", while the parameter provided was )" << param->mode << endl;
       exit(1);
    }

    std::ostringstream oss2;
    oss2 << param->dot_file << "final";
    std::ostringstream oss3;
    oss3 << param->dot_file << "sinksfinal";

    OUTPUT = "both";

    if(OUTPUT == "dot" || OUTPUT == "both") {
      oss2 << ".dot";
      ofstream output(oss2.str().c_str());
      merger.todot();
      output << merger.dot_output;
      output.close();
    }

    if(OUTPUT == "json" || OUTPUT == "both") {
        oss2 << ".json";
        ofstream output(oss2.str().c_str());
        merger.tojson();
        output << merger.json_output;
        output.close();

        oss3 << ".json";
        ofstream output2(oss3.str().c_str());
        merger.tojsonsinks();
        output2 << merger.json_output;
        output2.close();
    }

} // end run

#ifndef UNIT_TESTING
int main(int argc, char *argv[]){

    loguru::g_stderr_verbosity = loguru::Verbosity_OFF;
    loguru::init(argc, argv);
    loguru::add_file("flexfringe.log", loguru::Append, loguru::Verbosity_MAX);
    LOG_S(INFO) << "Starting flexfringe run";

    auto* param = new parameters();

    for(int i = 0; i < argc; i++) {
      param->command += string(argv[i]) + string(" ");
    }

    cout << "welcome, running git commit " << gitversion <<  " with: "<< param->command << endl;

    // CLI11-based parameter parsing to replace libpopt
    // https://cliutils.gitlab.io/CLI11Tutorial/
    CLI::App app{"flexFringe"
                 "Copyright 2020 Sicco Verwer, Delft University of Technology"
                 "with contributions from Christian Hammerschmidt, Delft University of Technology,\
                  University of Luxembourg"
                 "with contributions from APTA Technologies B.V."
                 "based on"
                 "DFASAT with random greedy preprocessing"
                 "Copyright 2015 Sicco Verwer and Marijn Heule, Delft University of Technology."
    };

    // remove -h short-form from help because for backward-compatibility --heuristic needs nit
    app.set_help_flag("--help", "Print this help message and exit");

    // read parameters from ini file if desired
    std::string default_file_name = "flexfringe.ini";
    app.set_config("-i,--ini", default_file_name, "Read an ini file", false);

    // { "debug", 'V', POPT_ARG_INT, &(param->debugging), 'V', "Debug mode and verbosity evel", "integer" },
    param->debugging = 0;
    app.add_option("-V,--debug", param->debugging, "A help string");

    // { "output-dir", 'o', POPT_ARG_STRING, &(dot_file), 'o', "Relative path for output files with trailing /, default: \"./\".", "string" },
    param->output = "./";
    app.add_option("-o,--output-dir", param->dot_file, "Relative path for output files with trailing /, default: \"./\".");

    //{ "output", 'T', POPT_ARG_STRING, &(output), 'T', "Switch between output in dot, json (default), or both formats" , "string" },
    string output = "both";
    app.add_option("-T,--output", output, "Switch between output in dot, json, or both (default) formats.");

    // { "heuristic-name", 'h', POPT_ARG_STRING, &(hName), 'h', "Name of the merge heurstic to use; default count_driven. Use any heuristic in the evaluation directory. It is often beneficial to write your own, as heuristics are very application specific.", "string" },
    app.add_option("-h,--heuristic-name", param->hName, "Name of the merge heurstic to use; default count_driven. Use any heuristic in the evaluation directory. It is often beneficial to write your own, as heuristics are very application specific.")->required();

    // { "data-name", 'd', POPT_ARG_STRING, &(hData), 'd', "Name of the merge data class to use; default count_data. Use any heuristic in the evaluation directory.", "string" },
    app.add_option("-d,--data-name", param->hData, "Name of the merge data class to use; default count_data. Use any heuristic in the evaluation directory.")->required();

    // { "mode", 'M', POPT_ARG_STRING, &(mode), 'M', "batch or stream depending on the mode of operation", "string" },
    param->mode = "batch";
    app.add_option("-M,--mode", param->mode, "batch (default), interactive, or stream depending on the mode of operation.");

    // { "evalpar", 'X', POPT_ARG_STRING, &(evalpar), 'X', "string of key-value pairs", "string" }
    param->evalpar = "mode=default";
    app.add_option("-X,--evalpar", param->evalpar, "string of key-value pairs for evaluation functions");

    // { "method", 'm', POPT_ARG_INT, &(param->method), 'm', "Method to use when merging states, default value 1 is random greedy (used in Stamina winner), 2 is one standard (non-random) greedy.", "integer" },
    app.add_option("-m,--method", param->method, "Method to use when merging states, default value 1 is random greedy (used in Stamina winner), 2 is one standard (non-random) greedy.");

    // { "seed", 's', POPT_ARG_INT, &(param->seed), 's', "Seed for random merge heuristic; default=12345678", "integer" },
    app.add_option("-s,--seed", param->seed, "Seed for random merge heuristic; default=12345678");

    // { "runs", 'n', POPT_ARG_INT, &(param->runs), 'n', "Number of random greedy runs/iterations; default=1. Advice: when using random greedy, a higher value is recommended (100 was used in Stamina winner).\n\nSettings that modify the red-blue state-merging framework:", "integer" },
    app.add_option("-n,--runs", param->runs, "Number of random greedy runs/iterations; default=1. Advice: when using random greedy, a higher value is recommended (100 was used in Stamina winner).");

    //  { "extend", 'x', POPT_ARG_INT, &(param->extend), 'x', "When set to 1, any merge candidate (blue) that cannot be merged with any target (red) is immediately changed into a (red) target; default=1. If set to 0, a merge candidate is only changed into a target when no more merges are possible. Advice: unclear which strategy is best, when using statistical (or count-based) consistency checks, keep in mind that merge consistency between states may change due to other performed merges. This will especially influence low frequency states. When there are a lot of those, we therefore recommend setting x=0.", "integer" },
    app.add_option("-x,--extend", param->extend, "When set to 1, any merge candidate (blue) that cannot be merged with any target (red) is immediately changed into a (red) target; default=1. If set to 0, a merge candidate is only changed into a target when no more merges are possible. Advice: unclear which strategy is best, when using statistical (or count-based) consistency checks, keep in mind that merge consistency between states may change due to other performed merges. This will especially influence low frequency states. When there are a lot of those, we therefore recommend setting x=0.");

    // { "shallowfirst", 'w', POPT_ARG_INT, &(param->shallowfirst), 'w', "When set to 1, the ordering of the nodes is changed from most frequent first (default) to most shallow (smallest depth) first; default=0. Advice: use depth-first when learning from characteristic samples.", "integer" },
    app.add_option("-w,--shallowfirst", param->shallowfirst, "When set to 1, the ordering of the nodes is changed from most frequent first (default) to most shallow (smallest depth) first; default=0. Advice: use depth-first when learning from characteristic samples.");

    // { "largestblue", 'a', POPT_ARG_INT, &(param->largestblue), 'a', "When set to 1, the algorithm only tries to merge the most frequent (or most shallow if w=1) candidate (blue) states with any target (red) state, instead of all candidates; default=0. Advice: this reduces run-time significantly but comes with a potential decrease in merge quality.", "integer" },
    app.add_option("-a,--largestblue", param->largestblue, "When set to 1, the algorithm only tries to merge the most frequent (or most shallow if w=1) candidate (blue) states with any target (red) state, instead of all candidates; default=0. Advice: this reduces run-time significantly but comes with a potential decrease in merge quality");

    // { "blueblue", 'b', POPT_ARG_INT, &(param->blueblue), 'b', "When set to 1, the algorithm tries to merge candidate (blue) states with candiate (blue) states in addition to candidate (blue) target (red) merges; default=0. Advice: this adds run-time to the merging process in exchange for potential improvement in merge quality.", "integer" },
    app.add_option("-b,--blueblue", param->blueblue, "When set to 1, the algorithm tries to merge candidate (blue) states with candiate (blue) states in addition to candidate (blue) target (red) merges; default=0. Advice: this adds run-time to the merging process in exchange for potential improvement in merge quality.");

    // { "finalred", 'f', POPT_ARG_INT, &(param->finalred), 'f', "When set to 1, merges that add new transitions to red states are considered inconsistent. Merges with red states will also not modify any of the counts used in evaluation functions. Once a red state has been learned, it is considered final and unmodifiable; default=0. Advice: setting this to 1 frequently results in easier to vizualize and more insightful models.\n\nSettings that influece the use of sinks:", "integer" },
    app.add_option("-f,--finalred", param->finalred, "When set to 1, merges that add new transitions to red states are considered inconsistent. Merges with red states will also not modify any of the counts used in evaluation functions. Once a red state has been learned, it is considered final and unmodifiable; default=0. Advice: setting this to 1 frequently results in easier to vizualize and more insightful models.");

    // { "sinkson", 'I', POPT_ARG_INT, &(param->sinkson), 'I', "Set to 1 to use sink states; default=1. Advice: leads to much more concise and easier to vizualize models, but can cost predictive performance depending on the sink definitions.", "integer"},
    app.add_option("-I,--sinkson", param->sinkson, "Set to 1 to use sink states; default=1. Advice: leads to much more concise and easier to vizualize models, but can cost predictive performance depending on the sink definitions.");

    // { "sinkcount", 'Q', POPT_ARG_INT, &(param->sink_count), 'Q', "The maximum number of occurrences of a state for it to be a low count sink (see evaluation functions); default=10.", "integer" },
    app.add_option("-Q,--sinkcount", param->sink_count, "The maximum number of occurrences of a state for it to be a low count sink (see evaluation functions); default=10.");

    //   { "mergesinks", 'J', POPT_ARG_INT, &(param->mergesinks), 'J', "Sink nodes are candidates for merging during the greedy runs (setting 0 or 1); default=0. Advice: merging sinks typically only makes the learned model worse. Keep in mind that sinks can become non-sinks due to other merges that influcence the occurrence counts.", "integer" },
    app.add_option("-J,--mergesinks", param->mergesinks, "The maximum number of occurrences of a state for it to be a low count sink (see evaluation functions); default=10.");

    // { "satmergesinks", 'K', POPT_ARG_INT, &(param->satmergesinks), 'K', "Merge all sink nodes of the same type before sending the problem to the SAT solver (setting 0 or 1); default=1. Advice: radically improves runtime, only set to 0 when sinks of the same type can be different states in the final model.\n\nSettings that influence merge evaluations:", "integer" },
    app.add_option("-K,--satmergesinks", param->satmergesinks, "Merge all sink nodes of the same type before sending the problem to the SAT solver (setting 0 or 1); default=1. Advice: radically improves runtime, only set to 0 when sinks of the same type can be different states in the final model.");

    // { "testmerge", 't', POPT_ARG_INT, &(param->testmerge), 't', "When set to 1, merge tries in order to compute the evaluation scores do not actually perform the merges themselves. Thus the consistency and score evaluation for states in merges that add recursive loops are uninfluenced by earlier merges; default=0. Advice: setting this to 1 reduces run-time and can be useful when learning models using statistical evaluation functions, but can lead to inconsistencies when learning from labeled data.", "integer" },
    app.add_option("-t,--testmerge", param->testmerge, "When set to 1, merge tries in order to compute the evaluation scores do not actually perform the merges themselves. Thus the consistency and score evaluation for states in merges that add recursive loops are uninfluenced by earlier merges; default=0. Advice: setting this to 1 reduces run-time and can be useful when learning models using statistical evaluation functions, but can lead to inconsistencies when learning from labeled data.");

    // { "lowerbound", 'l', POPT_ARG_FLOAT, &(param->lower_bound), 'l', "Minimum value of the heuristic function, smaller values are treated as inconsistent, also used as the paramater value in any statistical tests; default=-1. Advice: state merging is forced to perform the merge with best heuristic value, it can sometimes be better to color a state red rather then performing a bad merge. This is achieved using a positive lower bound value. Models learned with positive lower bound are frequently more interpretable.", "float" },
    app.add_option("-l,--lowerbound", param->lower_bound, "Minimum value of the heuristic function, smaller values are treated as inconsistent, also used as the paramater value in any statistical tests; default=-1. Advice: state merging is forced to perform the merge with best heuristic value, it can sometimes be better to color a state red rather then performing a bad merge. This is achieved using a positive lower bound value. Models learned with positive lower bound are frequently more interpretable");

    // { "state_count", 'q', POPT_ARG_INT, &(param->state_count), 'q', "The minimum number of positive occurrences of a state for it to be included in overlap/statistical checks (see evaluation functions); default=25. Advice: low frequency states can have an undesired influence on statistical tests, set to at least 10. Note that different evaluation functions can use this parameter in different ways.", "integer" },
    app.add_option("-q,--state_count", param->state_count, "The minimum number of positive occurrences of a state for it to be included in overlap/statistical checks (see evaluation functions); default=25. Advice: low frequency states can have an undesired influence on statistical tests, set to at least 10. Note that different evaluation functions can use this parameter in different ways.");

    //   { "symbol_count", 'y', POPT_ARG_INT, &(param->symbol_count), 'y', "The minimum number of positive occurrences of a symbol/transition for it to be included in overlap/statistical checks, symbols with less occurrences are binned together; default=10. Advice: low frequency transitions can have an undesired influence on statistical tests, set to at least 4. Note that different evaluation functions can use this parameter in different ways.", "integer" },
    app.add_option("-y,--symbol_count", param->symbol_count, "When set to 1, merge tries in order to compute the evaluation scores do not actually perform the merges themselves. Thus the consistency and score evaluation for states in merges that add recursive loops are uninfluenced by earlier merges; default=0. Advice: setting this to 1 reduces run-time and can be useful when learning models using statistical evaluation functions, but can lead to inconsistencies when learning from labeled data.");

    // { "epsilon", 'e', POPT_ARG_FLOAT, &(param->epsilon), 'e', "epsilon for Hoeffding condition", "float" },
    app.add_option("-e,--epsilon", param->epsilon, "epsilon for Hoeffding condition");

    // { "delta", 'D', POPT_ARG_FLOAT, &(param->delta), 'D', "delta for Hoeffding condition", "float" },
    app.add_option("-D,--delta", param->delta, "delta for Hoeffding condition");

    // { "batchsize", 'B', POPT_ARG_INT, &(param->batchsize), 'B', "batchsize for streaming", "int" },
    app.add_option("-B,--batch_size", param->batchsize, "batchsize for streaming");

    // { "correction", 'c', POPT_ARG_FLOAT, &(param->correction), 'c', "Value of a Laplace correction (smoothing) added to all symbol counts when computing statistical tests (in ALERGIA, LIKELIHOODRATIO, AIC, and KULLBACK-LEIBLER); default=0.0. Advice: unclear whether smoothing is needed for the different tests, more smoothing typically leads to smaller models.", "float" },
    app.add_option("-c,--correction", param->correction, "Value of a Laplace correction (smoothing) added to all symbol counts when computing statistical tests (in ALERGIA, LIKELIHOODRATIO, AIC, and KULLBACK-LEIBLER); default=0.0. Advice: unclear whether smoothing is needed for the different tests, more smoothing typically leads to smaller models.");

    // { "extrapar", 'p', POPT_ARG_FLOAT, &(param->extrapar), 'p', "Extra parameter used during statistical tests, the significance level for the likelihood ratio test, the alpha value for ALERGIA; default=0.5. Advice: look up the statistical test performed, this parameter is not always the same as a p-value.\n\nSettings influencing the SAT solving procedures:", "float" },
    app.add_option("-p,--extrapar", param->extrapar, "Extra parameter used during statistical tests, the significance level for the likelihood ratio test, the alpha value for ALERGIA; default=0.5. Advice: look up the statistical test performed, this parameter is not always the same as a p-value.");

    // { "sataptabound", 'A', POPT_ARG_INT, &(param->sataptabound), 'A', "Maximum number of remaining states in the partially learned DFA before starting the SAT search process. The higher this value, the larger the problem sent to the SAT solver; default=2000. Advice: try sending problem instances that are as large as possible, since larger instances take more time, test what time is acceptable for you.", "integer" },
    app.add_option("-A,--sataptabound", param->sataptabound, "Maximum number of remaining states in the partially learned DFA before starting the SAT search process. The higher this value, the larger the problem sent to the SAT solver; default=2000. Advice: try sending problem instances that are as large as possible, since larger instances take more time, test what time is acceptable for you.");

    // { "satdfabound", 'D', POPT_ARG_INT, &(param->satdfabound), 'D', "Maximum size of the partially learned DFA before starting the SAT search process; default=50. Advice: when the merging process performs bad merges, it can blow-up the size of the learned model. This test ensres that models that are too large do not get solved.", "integer" },
    // this clashed with short option -D for --delta
    app.add_option("-U,--satdfabound", param->satdfabound, "Maximum size of the partially learned DFA before starting the SAT search process; default=50. Advice: when the merging process performs bad merges, it can blow-up the size of the learned model. This test ensres that models that are too large do not get solved.");

    // { "satextra", 'E', POPT_ARG_INT, &(param->satextra), 'E', "DFASAT runs a SAT solver to find a solution of size at most the size of the partially learned DFA + E; default=5. Advice: larger values greatly increases run-time. Setting it to 0 is frequently sufficient (when the merge heuristic works well).", "integer" },
    app.add_option("-E,--satextra", param->satextra, "DFASAT runs a SAT solver to find a solution of size at most the size of the partially learned DFA + E; default=5. Advice: larger values greatly increases run-time. Setting it to 0 is frequently sufficient (when the merge heuristic works well).");

    // { "satplus", 'P', POPT_ARG_INT, &(param->satplus), 'P', "With every iteration, DFASAT tries to find solutions of size at most the best solution found + P, default=0. Advice: current setting only searches for better solutions. If a few extra states is OK, set it higher.", "int" },
    app.add_option("-P,--satplus", param->satplus, "With every iteration, DFASAT tries to find solutions of size at most the best solution found + P, default=0. Advice: current setting only searches for better solutions. If a few extra states is OK, set it higher.");

    // { "satfinalred", 'F', POPT_ARG_INT, &(param->satfinalred), 'F', "Make all transitions from red states without any occurrences force to have 0 occurrences (similar to targeting a rejecting sink), (setting 0 or 1) before sending the problem to the SAT solver; default=0. Advice: the same as finalred but for the SAT solver. Setting it to 1 greatly improves solving speed.", "integer" },
    app.add_option("-F,--satfinalred", param->satfinalred, "Make all transitions from red states without any occurrences force to have 0 occurrences (similar to targeting a rejecting sink), (setting 0 or 1) before sending the problem to the SAT solver; default=0. Advice: the same as finalred but for the SAT solver. Setting it to 1 greatly improves solving speed.");

    // { "satsymmetry", 'Y', POPT_ARG_INT, &(param->symmetry), 'Y', "Add symmetry breaking predicates to the SAT encoding (setting 0 or 1), based on Ulyantsev et al. BFS symmetry breaking; default=1. Advice: in our experience this only improves solving speed.", "integer" },
    app.add_option("-Y,--symmetry", param->symmetry, "dd symmetry breaking predicates to the SAT encoding (setting 0 or 1), based on Ulyantsev et al. BFS symmetry breaking; default=1. Advice: in our experience this only improves solving speed.");

    // { "satonlyinputs", 'O', POPT_ARG_INT, &(param->forcing), 'O', "Add predicates to the SAT encoding that force transitions in the learned DFA to be used by input examples (setting 0 or 1); default=0. Advice: leads to non-complete models. When the data is sparse, this should be set to 1. It does make the instance larger and can have a negative effect on the solving time.", "integer" },
    app.add_option("-O,--forcing", param->forcing, "Add predicates to the SAT encoding that force transitions in the learned DFA to be used by input examples (setting 0 or 1); default=0. Advice: leads to non-complete models. When the data is sparse, this should be set to 1. It does make the instance larger and can have a negative effect on the solving time.");

    // input file name
    app.add_option("tracefile", param->dfa_file, "Name of the input file containing the traces, either in Abadingo or JSON format.")->required();

    CLI11_PARSE(app, argc, argv);


    run(param);

    delete param;

    LOG_S(INFO) << "Ending flexfringe run normally";

    return 0;
}
#endif
