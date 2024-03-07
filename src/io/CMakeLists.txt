set(LIBRARY_NAME io)

file(GLOB SRC *.cpp)
set(SRC_ORIG ${SRC})

add_library(${LIBRARY_NAME} SHARED ${SRC})
find_package(Boost 1.55 COMPONENTS serialization system)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${Boost_INCLUDE_DIR})
include_directories(../sqlite/include)

target_link_libraries(${LIBRARY_NAME} PUBLIC
	${Boost_LIBRARIES}
	Eigen3::Eigen
	)
