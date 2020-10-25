cmake_minimum_required(VERSION 3.15)
project(flexfringe)

set(CMAKE_CXX_STANDARD 14)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(source)
include_directories(source/evaluation)

add_compile_definitions(LOGURU_WITH_STREAMS=1)

add_executable(flexfringe
        source/utility/loguru.cpp
        source/evaluation/aike.cpp
        source/evaluation/aike.h
        source/evaluation/alergia.cpp
        source/evaluation/alergia.h
        source/evaluation/alergia94.cpp
        source/evaluation/alergia94.h
        source/evaluation/conflict-driven-edsm.cpp
        source/evaluation/conflict-driven-edsm.h
        source/evaluation/conflict-driven-overlap.cpp
        source/evaluation/conflict-driven-overlap.h
        source/evaluation/conflict-driven.cpp
        source/evaluation/conflict-driven.h
        source/evaluation/depth-driven.cpp
        source/evaluation/depth-driven.h
        source/evaluation/evaluate.cpp
        source/evaluation/evaluate.h
        source/evaluation/evidence-driven.cpp
        source/evaluation/evidence-driven.h
        source/evaluation/fixed_depth_regression.cpp
        source/evaluation/fixed_depth_regression.h
        source/evaluation/full-overlap-driven.cpp
        source/evaluation/full-overlap-driven.h
        source/evaluation/kldistance.cpp
        source/evaluation/kldistance.h
        source/evaluation/likelihood.cpp
        source/evaluation/likelihood.h
        source/evaluation/mealy.cpp
        source/evaluation/mealy.h
        source/evaluation/mse-error.cpp
        source/evaluation/mse-error.h
        source/evaluation/num_count.cpp
        source/evaluation/num_count.h
        source/evaluation/overlap-driven.cpp
        source/evaluation/overlap-driven.h
        source/evaluation/overlap4logs.cpp
        source/evaluation/overlap4logs.h
        source/evaluation/process-mining.cpp
        source/evaluation/process-mining.h
        source/evaluation/rti+.cpp
        source/evaluation/rti+.h
        source/evaluation/tsp.cpp
        source/evaluation/tsp.h
        source/apta.cpp
        source/apta.h
        source/common.h
        source/conflict_graph.cpp
        source/conflict_graph.h
        source/dfasat.cpp
        source/dfasat.h
        source/evaluation_factory.cpp
        source/evaluation_factory.h
        source/evaluators.h
        source/gitversion.cpp
        source/inputdata.cpp
        source/inputdata.h
        source/interactive.cpp
        source/interactive.h
        source/utility/json.hpp
        source/main.cpp
        source/parameters.cpp
        source/parameters.h
        source/random_greedy.cpp
        source/random_greedy.h
        source/refinement.cpp
        source/refinement.h
        source/searcher.cpp
        source/searcher.h
        source/state_merger.cpp
        source/state_merger.h
        source/stream.cpp
        source/stream.h
        source/mem_store.cpp
        source/mem_store.h)

#find_package(GSL REQUIRED)
#target_link_libraries(flexfringe GSL::gsl GSL::gslcblas)

#target_link_libraries(flexfringe m)

find_package(Threads)
target_link_libraries(flexfringe ${CMAKE_THREAD_LIBS_INIT}) # For pthreads

if(NOT WIN32)
    target_link_libraries(flexfringe dl) # For ldl
endif()
