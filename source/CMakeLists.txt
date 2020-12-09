include_directories("${CMAKE_CURRENT_SOURCE_DIR}/evaluation"
                    "${CMAKE_CURRENT_SOURCE_DIR}/utility"
                    )

add_subdirectory(evaluation)
add_subdirectory(utility)

## create the gitversion.cpp
execute_process(COMMAND "git" "rev-parse" "HEAD"
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                OUTPUT_VARIABLE Gitversion)
string(REPLACE "\n" "" Gitversion ${Gitversion})
# TODO: check version and append either version or "No commit info available\";"
file(WRITE "gitversion.cpp" "const char *gitversion = \"${Gitversion}\";")

## create the evaluators.h file
file(GLOB Files "${CMAKE_CURRENT_SOURCE_DIR}/evaluation/*.h")
file(WRITE "evaluators.h" "#ifndef __ALL_HEADERS__ \n#define __ALL_HEADERS__ \n\n")
foreach(Filepath ${Files})
    get_filename_component(Filename ${Filepath} NAME)
    file( APPEND "evaluators.h" "#include \"${Filename}\"\n")
endforeach()
file(APPEND "evaluators.h" "\n#endif")

add_library(Source STATIC
            apta.cpp
            apta.h
            common.h
            conflict_graph.cpp
            conflict_graph.h
            dfasat.cpp
            dfasat.h
            evaluation_factory.cpp
            evaluation_factory.h
            evaluators.h
            gitversion.cpp
            inputdata.cpp
            inputdata.h
            interactive.cpp
            interactive.h
            parameters.cpp
            parameters.h
            random_greedy.cpp
            random_greedy.h
            refinement.cpp
            refinement.h
            searcher.cpp
            searcher.h
            state_merger.cpp
            state_merger.h
            stream.cpp
            stream.h
            mem_store.cpp
            mem_store.h
            )

#target_include_directories(Source PUBLIC 
#                            "${CMAKE_CURRENT_SOURCE_DIR}/evaluation"
#                            "${CMAKE_CURRENT_SOURCE_DIR}/utility"
#                            )
