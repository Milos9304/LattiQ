#set(LIBRARY_NAME src)

#file(GLOB SRC *.cpp)
			
#set(SRC_ORIG ${SRC})

add_subdirectory(io)
#add_subdirectory(accelerator)
#add_subdirectory(optimizer)
#add_subdirectory(circuit)
#add_subdirectory(ansatz)
#add_subdirectory(hamiltonian)
#add_subdirectory(qaoa)
#add_subdirectory(vqe)
#add_subdirectory(lattice)

#add_library(${LIBRARY_NAME} SHARED ${SRC})
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../bin")
#add_executable(fastVQA ${SRC_ORIG})
#target_link_libraries(${LIBRARY_NAME} PRIVATE lattice)#io vqe qaoa hamiltonian accelerator optimizer lattice fplll gmp gmpxx QuEST)
