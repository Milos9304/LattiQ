cmake_minimum_required (VERSION 3.17 FATAL_ERROR)

project (LattiQ
	LANGUAGES CXX
	DESCRIPTION "Experimental framework for lattice problems on variational quantum algorithms using FastVQA"
	VERSION 0.1.0
)

set(CMAKE_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")
set(DISTRIBUTED OFF)
message(
	"DISTRIBUTED=${DISTRIBUTED}"
)

find_package(FastVQA)

find_package(MPI)
if(MPI_FOUND AND DISTRIBUTED)
	set(DISTRIBUTED_OPTS_FILES ${PROJECT_SOURCE_DIR}/src/lattiq_distributed.cpp)
else()
	set(DISTRIBUTED_OPTS_FILES ${PROJECT_SOURCE_DIR}/src/lattiq.cpp)
endif()

file(GLOB SRC_ORIG
        "${DISTRIBUTED_OPTS_FILES}"
        ${PROJECT_SOURCE_DIR}/src/lattice/lattice.cpp
        ${PROJECT_SOURCE_DIR}/src/lattice/lattice_utils.cpp
        ${PROJECT_SOURCE_DIR}/src/run_paper_experiment.cpp
        ${PROJECT_SOURCE_DIR}/src/symbolic_manipulation.cpp
)
       
add_subdirectory(src)
       
add_executable(${PROJECT_NAME} ${SRC_ORIG})
include_directories(
	${PROJECT_SOURCE_DIR}/src
)
target_link_libraries(${PROJECT_NAME} fastVQA::FastVQA fplll gmpxx gmp) #lattice