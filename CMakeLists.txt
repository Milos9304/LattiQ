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
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

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

file(GLOB SRC_ORIG_AVERAGEDSVP
	${PROJECT_SOURCE_DIR}/src/lattice/lattice.cpp
    ${PROJECT_SOURCE_DIR}/src/lattice/lattice_utils.cpp
	${PROJECT_SOURCE_DIR}/src/averaged_svp.cpp
	${PROJECT_SOURCE_DIR}/src/averaged_svp_utils.cpp
	${PROJECT_SOURCE_DIR}/src/experiment_runner.cpp
)
       
add_subdirectory(src)
       
add_executable(${PROJECT_NAME} ${SRC_ORIG})
add_executable(averagedSvp ${SRC_ORIG_AVERAGEDSVP} ${PROJECT_SOURCE_DIR}/src/experiment_runner.cpp)


include_directories(
	${PROJECT_SOURCE_DIR}/src
)
target_link_libraries(${PROJECT_NAME} fastVQA::FastVQA Eigen3::Eigen fplll gmpxx gmp) #lattice
target_compile_definitions(${PROJECT_NAME}
    PRIVATE
    QuEST_PREC=4
)

target_link_libraries(averagedSvp fastVQA::FastVQA Eigen3::Eigen fplll gmpxx gmp)
target_compile_definitions(averagedSvp
    PRIVATE
    QuEST_PREC=4
)
