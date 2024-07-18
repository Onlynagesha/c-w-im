cmake_minimum_required(VERSION 3.20)
project(C_w_IM)

set(CMAKE_CXX_STANDARD 23)
# Required by reflect-cpp for huge structs
set(CMAKE_CXX_FLAGS -fconstexpr-depth=4096)

find_package(argparse CONFIG REQUIRED)
find_package(Boost COMPONENTS unit_test_framework REQUIRED) # Boost.assert
find_package(fmt CONFIG REQUIRED)
find_package(magic_enum CONFIG REQUIRED)
find_package(nameof CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED) # Dependency of NWGraph library

set(THIRD_PARTY_INCLUDE_DIRS 
    third_party/nwgraph/include
    third_party/reflect_cpp/include
    third_party/yalantinglibs/include)

set(COMMON_INCLUDE_DIRS
    ${Boost_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${THIRD_PARTY_INCLUDE_DIRS})

set(COMMON_LINK_LIBRARIES
    argparse::argparse
    fmt::fmt
    magic_enum::magic_enum
    nameof::nameof
    nlohmann_json::nlohmann_json
    TBB::tbb
    TBB::tbbmalloc)

# Enables customized Boost.assert fail handler (see utils/boost_assert.h for details)
set(COMMON_COMPILE_DEFINITIONS
    BOOST_ENABLE_ASSERT_DEBUG_HANDLER)
# Enables assertion for GCC libstdc++
if (CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
    list(APPEND COMMON_COMPILE_DEFINITIONS 
        _GLIBCXX_ASSERTIONS
        ENABLES_MYLOG_TRACE 
        DEBUG_FIXED_RANDOM_SEED)
endif ()

set(COMMON_COMPILE_OPTIONS)
# AddressSanitizer checking for debug mode
if (CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
    list(APPEND COMMON_COMPILE_OPTIONS -fsanitize=address -fno-omit-frame-pointer)
    list(APPEND COMMON_LINK_LIBRARIES -fsanitize=address -fno-omit-frame-pointer) # Required as link libraries
endif ()

# ---- Components ----

set(COMPONENT_OBJECTS
    yyjson-components
    io-components
    experiment-components
    contrast-components)

add_library(yyjson-components OBJECT
    third_party/reflect_cpp/src/yyjson.c)

add_library(io-components OBJECT
    argparse.cpp # Requires yyjson components
    dump.cpp
    graph_types.cpp)
target_link_libraries(io-components PRIVATE yyjson-components)

add_library(experiment-components OBJECT 
    experiments.cpp 
    experiments_internal_states.cpp # Requires yyjson components
    coarsening.cpp
    coarsening_dump.cpp
    wim.cpp
    wim_dump.cpp)
target_link_libraries(experiment-components PRIVATE yyjson-components)

add_library(contrast-components OBJECT
    contrast/imrank.cpp
    contrast/max_degree.cpp
    contrast/pagerank.cpp
    contrast/r_robust_scc.cpp)

# ---- Experiment Executables ----

set(EXPERIMENT_EXECUTABLES
    wim-experiment
    wbim-experiment
    wim-coarsening-experiment
    wbim-coarsening-experiment
    wim-contrast-experiment)

add_executable(wim-experiment
    exec/experiments/wim_experiment.cpp)

add_executable(wbim-experiment
    exec/experiments/wbim_experiment.cpp)

add_executable(wim-coarsening-experiment
    exec/experiments/wim_coarsening_experiment.cpp)

add_executable(wbim-coarsening-experiment
    exec/experiments/wbim_coarsening_experiment.cpp)

add_executable(wim-contrast-experiment
    exec/experiments/wim_contrast_experiment.cpp)
target_link_libraries(wim-contrast-experiment PRIVATE contrast-components)

foreach (target ${EXPERIMENT_EXECUTABLES})
    target_link_libraries(${target} PRIVATE yyjson-components io-components experiment-components)
endforeach()

# ---- Dataset IO Executables ----

set(DATASET_IO_EXECUTABLES
    create-wim-dataset
    create-wbim-dataset
    to-matrix-market)

add_executable(create-wim-dataset
    exec/datasets/create_wim_dataset.cpp
    create_dataset.cpp)

add_executable(create-wbim-dataset
    exec/datasets/create_wbim_dataset.cpp
    create_dataset.cpp)

add_executable(to-matrix-market
    exec/datasets/to_matrix_market.cpp
    to_matrix_market.cpp)

foreach (target ${DATASET_IO_EXECUTABLES})
    target_link_libraries(${target} PRIVATE yyjson-components io-components)
endforeach()

# ---- Tests ----

set(TEST_CASE_EXECUTABLES 
    test-graph-connectivity
    test-imrank
    test-max-degree
    test-pagerank
    test-wbim-activation-probability)

add_executable(test-graph-connectivity
    tests/graph_connectivity.cpp)

add_executable(test-imrank
    tests/imrank.cpp)

add_executable(test-max-degree
    tests/max_degree.cpp)

add_executable(test-pagerank
    tests/pagerank.cpp)

add_executable(test-wbim-activation-probability
    tests/wbim_activation_probability.cpp)
target_link_libraries(test-wbim-activation-probability PRIVATE 
    yyjson-components io-components experiment-components contrast-components)

foreach (target test-imrank test-max-degree test-pagerank)
    target_link_libraries(${target} PRIVATE contrast-components)
endforeach()

foreach (target ${TEST_CASE_EXECUTABLES})
    target_link_libraries(${target} PRIVATE ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    target_compile_definitions(${target} PRIVATE ENABLES_MYLOG_TRACE)
endforeach()

# ---- Playground (for Testing) ----

set(PLAYGROUND_TESTS 
    test-wim-coarsening 
    test-wbim-coarsening 
    test-prr-sketching)

add_executable(test-wim-coarsening
    playground/test_wim_coarsening.cpp)

add_executable(test-wbim-coarsening
    playground/test_wbim_coarsening.cpp)

add_executable(test-prr-sketching
    playground/test_prr_sketching.cpp)

foreach (target ${PLAYGROUND_TESTS})
    target_link_libraries(${target} PRIVATE yyjson-components io-components experiment-components)
endforeach()

# ---- Common ----

foreach(target ${COMPONENT_OBJECTS}
               ${EXPERIMENT_EXECUTABLES}
               ${DATASET_IO_EXECUTABLES}
               ${TEST_CASE_EXECUTABLES}
               ${PLAYGROUND_TESTS})
    target_compile_options(${target} PRIVATE ${COMMON_COMPILE_OPTIONS})
    target_compile_definitions(${target} PRIVATE ${COMMON_COMPILE_DEFINITIONS})
    target_include_directories(${target} PRIVATE ${COMMON_INCLUDE_DIRS})
    target_link_libraries(${target} PRIVATE ${COMMON_LINK_LIBRARIES})
endforeach()
