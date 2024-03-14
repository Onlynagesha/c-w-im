cmake_minimum_required(VERSION 3.20)
project(C_w_IM)

set(CMAKE_CXX_STANDARD 23)
# Required by reflect-cpp for huge structs
set(CMAKE_CXX_FLAGS -fconstexpr-depth=4096)

find_package(argparse CONFIG REQUIRED)
find_package(Boost REQUIRED) # Boost.assert
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
    list(APPEND COMMON_COMPILE_DEFINITIONS _GLIBCXX_ASSERTIONS)
endif ()

set(COMMON_COMPILE_OPTIONS)
# AddressSanitizer checking for debug mode
if (CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
    list(APPEND COMMON_COMPILE_OPTIONS -fsanitize=address -fno-omit-frame-pointer)
    list(APPEND COMMON_LINK_LIBRARIES -fsanitize=address -fno-omit-frame-pointer) # Required as link libraries
endif ()

# ---- Components ----

set(COMPONENT_OBJECTS
    io-components
    experiment-components)

add_library(io-components OBJECT
    third_party/reflect_cpp/src/yyjson.c
    argparse.cpp
    dump.cpp
    graph_types.cpp # Graph IO
    )

add_library(experiment-components OBJECT
    experiments.cpp
    coarsening.cpp
    coarsening_dump.cpp
    contrast_algorithms.cpp
    wim.cpp
)

# ---- Executables ----

set(EXPERIMENT_EXECUTABLES
    wim-experiment
    wim-coarsening-experiment
    wim-contrast-experiment)

add_executable(wim-experiment
    main_wim_experiment.cpp)

add_executable(wim-coarsening-experiment
    main_wim_coarsening_experiment.cpp)

add_executable(wim-contrast-experiment
    main_wim_contrast_experiment.cpp)

foreach (target ${EXPERIMENT_EXECUTABLES})
    target_link_libraries(${target} PRIVATE io-components experiment-components)
endforeach()

add_executable(create-dataset
    main_create_dataset.cpp
    create_dataset.cpp)

add_executable(to-matrix-market
    main_to_matrix_market.cpp
    to_matrix_market.cpp)

foreach (target create-dataset to-matrix-market)
    target_link_libraries(${target} PRIVATE io-components)
endforeach()

add_executable(test-coarsening
    playground/test_coarsening.cpp
    coarsening.cpp
    coarsening_dump.cpp
    wim.cpp)

add_executable(test-ranking-algorithms
    playground/test_ranking_algorithms.cpp
    contrast_algorithms.cpp)

add_executable(test-seed-detection
    playground/test_seed_detection.cpp
    wim.cpp)

foreach (target test-coarsening test-ranking-algorithms test-seed-detection)
    target_link_libraries(${target} PRIVATE io-components)
    target_compile_definitions(${target} PRIVATE ENABLES_MYLOG_TRACE)
endforeach()

foreach(target ${COMPONENT_OBJECTS}
               ${EXPERIMENT_EXECUTABLES}
               create-dataset to-matrix-market
               test-coarsening test-ranking-algorithms test-seed-detection)
    target_compile_options(${target} PRIVATE ${COMMON_COMPILE_OPTIONS})
    target_compile_definitions(${target} PRIVATE ${COMMON_COMPILE_DEFINITIONS})
    target_include_directories(${target} PRIVATE ${COMMON_INCLUDE_DIRS})
    target_link_libraries(${target} PRIVATE ${COMMON_LINK_LIBRARIES})
endforeach()
