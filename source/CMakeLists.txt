include_directories("${CMAKE_CURRENT_SOURCE_DIR}/evaluation"
                    "${CMAKE_CURRENT_SOURCE_DIR}/utility"
                    )

add_subdirectory(evaluation)
add_subdirectory(utility)

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
