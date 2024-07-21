cmake_minimum_required(VERSION 3.20)
project(C_w_IM)

set(CMAKE_CXX_STANDARD 23)
# Required by reflect-cpp for huge structs
set(CMAKE_CXX_FLAGS -fconstexpr-depth=4096)

find_package(argparse CONFIG REQUIRED)
find_package(Boost CONFIG COMPONENTS unit_test_framework REQUIRED) # Boost.assert
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

add_library(dump-components SHARED
    third_party/reflect_cpp/src/yyjson.c
    argparse.cpp
    coarsening_dump.cpp
    dump.cpp
    wim_dump.cpp)

set(CONTRAST_SOURCE_FILES
    contrast/imrank.cpp
    contrast/max_degree.cpp
    contrast/pagerank.cpp
    contrast/r_robust_scc.cpp
    wim.cpp) # Required by r-Robust SCC algorithm

set(EXPERIMENT_COMMON_SOURCE_FILES
    graph_types.cpp # For graph IO
    experiments/internal_states.cpp
    experiments/components.cpp
    coarsening.cpp # Requires by experiments/components.cpp
    contrast/max_degree.cpp) # Required by generate_wbim_seeds in experiments/components.cpp

# ---- Experiment Executables ----

set(EXPERIMENT_EXECUTABLES
    wim-experiment
    wbim-experiment
    wim-coarsening-experiment
    wbim-coarsening-experiment
    wim-contrast-experiment)

add_executable(wim-experiment
    exec/experiments/wim_experiment.cpp
    ${EXPERIMENT_COMMON_SOURCE_FILES}
    experiments/impl/wim.cpp
    wim.cpp
    )

add_executable(wbim-experiment
    exec/experiments/wbim_experiment.cpp
    ${EXPERIMENT_COMMON_SOURCE_FILES}
    experiments/impl/wim.cpp
    wim.cpp)

add_executable(wim-coarsening-experiment
    exec/experiments/wim_coarsening_experiment.cpp
    ${EXPERIMENT_COMMON_SOURCE_FILES}
    experiments/impl/wim.cpp
    experiments/impl/wim_coarsening.cpp
    coarsening.cpp
    wim.cpp)

add_executable(wbim-coarsening-experiment
    exec/experiments/wbim_coarsening_experiment.cpp
    ${EXPERIMENT_COMMON_SOURCE_FILES}
    experiments/impl/wim.cpp
    experiments/impl/wim_coarsening.cpp
    coarsening.cpp
    wim.cpp)

add_executable(wim-contrast-experiment
    exec/experiments/wim_contrast_experiment.cpp
    ${EXPERIMENT_COMMON_SOURCE_FILES}
    ${CONTRAST_SOURCE_FILES}
    experiments/impl/wim.cpp
    experiments/impl/wim_coarsening.cpp
    experiments/impl/contrast.cpp
    coarsening.cpp)

foreach (target ${EXPERIMENT_EXECUTABLES})
    target_link_libraries(${target} PRIVATE dump-components)
endforeach()

# ---- Dataset IO Executables ----

set(DATASET_IO_EXECUTABLES
    create-wim-dataset
    create-wbim-dataset
    to-matrix-market)

add_executable(create-wim-dataset
    exec/datasets/create_wim_dataset.cpp
    create_dataset.cpp
    graph_types.cpp)

add_executable(create-wbim-dataset
    exec/datasets/create_wbim_dataset.cpp
    create_dataset.cpp
    graph_types.cpp)

add_executable(to-matrix-market
    exec/datasets/to_matrix_market.cpp
    to_matrix_market.cpp
    graph_types.cpp)

foreach (target ${DATASET_IO_EXECUTABLES})
    target_link_libraries(${target} PRIVATE dump-components)
endforeach()

# ---- Tests ----

set(TEST_CASE_EXECUTABLES 
    test-graph-connectivity
    test-imrank
    test-max-degree
    test-pagerank
    test-radix-sort
    test-wbim-activation-probability)

add_executable(test-graph-connectivity
    tests/graph_connectivity.cpp)

add_executable(test-imrank
    tests/imrank.cpp
    contrast/max_degree.cpp # Required by contrast/imrank.cpp
    contrast/imrank.cpp)

add_executable(test-max-degree
    tests/max_degree.cpp
    contrast/max_degree.cpp)

add_executable(test-pagerank
    tests/pagerank.cpp
    contrast/pagerank.cpp)

add_executable(test-radix-sort
    tests/radix_sort.cpp)

add_executable(test-wbim-activation-probability
    tests/wbim_activation_probability.cpp
    wim.cpp)
target_link_libraries(test-wbim-activation-probability PRIVATE dump-components)

foreach (target ${TEST_CASE_EXECUTABLES})
    target_link_libraries(${target} PRIVATE Boost::unit_test_framework)
    target_compile_definitions(${target} PRIVATE ENABLES_MYLOG_TRACE)
endforeach()

# ---- Playground (for Testing) ----

set(PLAYGROUND_TESTS 
    test-wim-coarsening 
    test-wbim-coarsening 
    test-prr-sketching)

add_executable(test-wim-coarsening
    playground/test_wim_coarsening.cpp
    coarsening.cpp
    wim.cpp)

add_executable(test-wbim-coarsening
    playground/test_wbim_coarsening.cpp
    coarsening.cpp
    wim.cpp)

add_executable(test-prr-sketching
    playground/test_prr_sketching.cpp
    wim.cpp)

foreach (target ${PLAYGROUND_TESTS})
    target_link_libraries(${target} PRIVATE dump-components)
endforeach()

# ---- Common ----

foreach(target dump-components
               ${EXPERIMENT_EXECUTABLES}
               ${DATASET_IO_EXECUTABLES}
               ${TEST_CASE_EXECUTABLES}
               ${PLAYGROUND_TESTS})
    target_compile_options(${target} PRIVATE ${COMMON_COMPILE_OPTIONS})
    target_compile_definitions(${target} PRIVATE ${COMMON_COMPILE_DEFINITIONS})
    target_include_directories(${target} PRIVATE ${COMMON_INCLUDE_DIRS})
    target_link_libraries(${target} PRIVATE ${COMMON_LINK_LIBRARIES})
endforeach()
