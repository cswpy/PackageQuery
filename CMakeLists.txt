cmake_minimum_required(VERSION 3.12)

option(CXX "enable C++ compilation" ON)
if (CXX)
  enable_language(CXX)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILE_FEATURES cxx_std_17)

project(PackageQuery VERSION 1.0.0)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(EXECUTABLE_NAME "main")

# MAIN
#set(MAIN_NAME "main_simplex")
#set(MAIN_NAME "main_mc_matrix")
#set(MAIN_NAME "main_lattice_solver")
#set(MAIN_NAME "main_ccfits")
set(MAIN_NAME "main_reducer")

# TEST 
#set(MAIN_NAME "test_accuracy_lattice_walk")
#set(MAIN_NAME "test_speed_lattice_walk")
#set(MAIN_NAME "test_speed_simplex")
#set(MAIN_NAME "test_lattice_solver")
#set(MAIN_NAME "test_ppq")

message("PROJECT NAME: " ${PROJECT_NAME})
message("SOURCE FOLDER: " ${CMAKE_SOURCE_DIR})
message("BUILD FOLDER: " ${CMAKE_BINARY_DIR})
message("CMAKE MODULE FOLDER: " ${CMAKE_MODULE_PATH})
message("EXECUTABLE NAME: " ${MAIN_NAME})

include_directories("${PROJECT_SOURCE_DIR}")
include_directories(
	${PROJECT_SOURCE_DIR}/include
	${PROJECT_SOURCE_DIR}/src
)

file(GLOB all_SRCS
  "${CMAKE_SOURCE_DIR}/runnable/${MAIN_NAME}.cpp"
  "${PROJECT_SOURCE_DIR}/include/*.h"
  "${PROJECT_SOURCE_DIR}/include/*.hpp"
  "${PROJECT_SOURCE_DIR}/src/*.cpp"
  "${PROJECT_SOURCE_DIR}/src/*.c"
)

#OpenMP Deprecated
# find_package(OpenMP)
# if (OPENMP_FOUND)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Conan
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

# Gurobi
find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIR})

# Cplex
find_package(Cplex)
include_directories(${CPLEX_INCLUDE_DIRS})

add_executable(${EXECUTABLE_NAME} ${all_SRCS})

# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

if(CXX)
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
endif()

target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${CONAN_LIBS} ${GUROBI_LIBRARY} ${CMAKE_DL_LIBS})
target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${CPLEX_LIBRARIES})