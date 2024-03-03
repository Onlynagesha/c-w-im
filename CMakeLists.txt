cmake_minimum_required(VERSION 3.20)
project(C_w_IM)

set(CMAKE_CXX_STANDARD 23)

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
set(THIRD_PARTY_SOURCE_FILES 
    third_party/reflect_cpp/src/yyjson.c)

# Enables customized Boost.assert fail handler (see utils/boost_assert.h for details)
set(COMMON_COMPILE_DEFINITIONS
    BOOST_ENABLE_ASSERT_DEBUG_HANDLER)
# Enables assertion for GCC libstdc++
if (CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
    list(APPEND COMMON_COMPILE_DEFINITIONS _GLIBCXX_ASSERTIONS)
endif ()

set(COMMON_INCLUDE_DIRS
    ${Boost_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${THIRD_PARTY_INCLUDE_DIRS})

set(COMMON_SOURCE_FILES
    ${THIRD_PARTY_SOURCE_FILES}
    argparse.cpp
    dump.cpp
    graph_types.cpp
    )
set(COMMON_LINK_LIBRARIES
    argparse::argparse
    fmt::fmt
    magic_enum::magic_enum
    nameof::nameof
    nlohmann_json::nlohmann_json
    TBB::tbb
    TBB::tbbmalloc)

add_executable(wim-experiment
    main_wim_experiment.cpp
    coarsening.cpp
    coarsening_dump.cpp
    experiments.cpp
    wim.cpp
    ${COMMON_SOURCE_FILES})

add_executable(wim-coarsening-experiment
    main_wim_coarsening_experiment.cpp
    coarsening.cpp
    coarsening_dump.cpp
    experiments.cpp
    wim.cpp
    ${COMMON_SOURCE_FILES})

add_executable(create-dataset
    main_create_dataset.cpp
    create_dataset.cpp
    ${COMMON_SOURCE_FILES})

add_executable(to-matrix-market
    main_to_matrix_market.cpp
    to_matrix_market.cpp
    ${COMMON_SOURCE_FILES})

add_executable(playground
    playground/coarsening.cpp
    coarsening.cpp
    coarsening_dump.cpp
    wim.cpp
    ${COMMON_SOURCE_FILES})
target_compile_definitions(playground PRIVATE ENABLES_MYLOG_TRACE)

foreach(target wim-experiment wim-coarsening-experiment create-dataset to-matrix-market playground)
    target_compile_definitions(${target} PRIVATE ${COMMON_COMPILE_DEFINITIONS})
    target_include_directories(${target} PRIVATE ${COMMON_INCLUDE_DIRS})
    target_link_libraries(${target} PRIVATE ${COMMON_LINK_LIBRARIES})

    # AddressSanitizer checking for debug mode
    if (CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
        target_compile_options(${target} PRIVATE -fsanitize=address -fno-omit-frame-pointer)
        target_link_libraries(${target} PRIVATE -fsanitize=address -fno-omit-frame-pointer)
    endif ()
endforeach()
