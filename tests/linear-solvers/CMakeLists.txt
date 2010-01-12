project(solver)

include(CMake.vars OPTIONAL)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})
include_directories(${hermes3d_SOURCE_DIR})

if(REAL)
	add_executable(${PROJECT_NAME}
		main.cc
		${HERMES_COMMON_DIR}/trace.cc
	)
	target_link_libraries(${PROJECT_NAME} ${TESTS_LIBRARIES})
	target_link_libraries(${PROJECT_NAME} ${HERMES_REAL_BIN})
endif(REAL)
if(COMPLEX)
	add_executable(${PROJECT_NAME}-cplx
		main.cc
		${HERMES_COMMON_DIR}/trace.cc
	)
	set_target_properties(${PROJECT_NAME}-cplx PROPERTIES COMPILE_FLAGS -DCOMPLEX)
	target_link_libraries(${PROJECT_NAME}-cplx ${TESTS_LIBRARIES})
	target_link_libraries(${PROJECT_NAME}-cplx ${HERMES_CPLX_BIN})
endif(COMPLEX)

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Tests

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})
if(REAL)
if(WITH_PETSC)
	add_test(petsc-solver-1 sh -c "${BIN} petsc ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
#	add_test(petsc-solver-2 sh -c "${BIN} petsc ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-singular")
endif(WITH_PETSC)

if(WITH_UMFPACK)
	add_test(umfpack-solver-1 sh -c "${BIN} umfpack ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(umfpack-solver-2 sh -c "${BIN} umfpack ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
	add_test(umfpack-solver-3 sh -c "${BIN} umfpack ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
	add_test(umfpack-solver-4 sh -c "${BIN} umfpack ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-singular")
	
	add_test(umfpack-solver-b-1 sh -c "${BIN} umfpack-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(umfpack-solver-b-2 sh -c "${BIN} umfpack-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
	add_test(umfpack-solver-b-3 sh -c "${BIN} umfpack-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
	add_test(umfpack-solver-b-4 sh -c "${BIN} umfpack-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-singular")
endif(WITH_UMFPACK)

if(WITH_PARDISO)
	add_test(pardiso-solver-1 sh -c "${BIN} pardiso ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(pardiso-solver-2 sh -c "${BIN} pardiso ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
	add_test(pardiso-solver-3 sh -c "${BIN} pardiso ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
	# NOTE: PARDISO is not able to recognize singular system	
#	add_test(pardiso-solver-4 sh -c "${BIN} pardiso ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-singular")
	
	add_test(pardiso-solver-b-1 sh -c "${BIN} pardiso-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(pardiso-solver-b-2 sh -c "${BIN} pardiso-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
	add_test(pardiso-solver-b-3 sh -c "${BIN} pardiso-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
	# NOTE: PARDISO is not able to recognize singular system	
#	add_test(pardiso-solver-b-4 sh -c "${BIN} pardiso-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-singular")
endif(WITH_PARDISO)

if(WITH_TRILINOS)
	# TODO: if we do not have aztecoo, do not run these tests
	add_test(aztecoo-solver-1 sh -c "${BIN} aztecoo ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(aztecoo-solver-2 sh -c "${BIN} aztecoo ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
	add_test(aztecoo-solver-3 sh -c "${BIN} aztecoo ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")

	add_test(aztecoo-solver-b-1 sh -c "${BIN} aztecoo-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(aztecoo-solver-b-2 sh -c "${BIN} aztecoo-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
	add_test(aztecoo-solver-b-3 sh -c "${BIN} aztecoo-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")

	# TODO: if we do not have amesos, do not run these tests
	add_test(amesos-solver-1 sh -c "${BIN} amesos ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(amesos-solver-2 sh -c "${BIN} amesos ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
	add_test(amesos-solver-3 sh -c "${BIN} amesos ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")

	add_test(amesos-solver-b-1 sh -c "${BIN} amesos-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(amesos-solver-b-2 sh -c "${BIN} amesos-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
	add_test(amesos-solver-b-3 sh -c "${BIN} amesos-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
endif(WITH_TRILINOS)

if(WITH_MUMPS)
	add_test(mumps-solver-1 sh -c "${BIN} mumps ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(mumps-solver-2 sh -c "${BIN} mumps ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
	add_test(mumps-solver-3 sh -c "${BIN} mumps ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
	add_test(mumps-solver-4 sh -c "${BIN} mumps ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-singular")
	
	add_test(mumps-solver-b-1 sh -c "${BIN} mumps-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(mumps-solver-b-2 sh -c "${BIN} mumps-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
	add_test(mumps-solver-b-3 sh -c "${BIN} mumps-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
	add_test(mumps-solver-b-4 sh -c "${BIN} mumps-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-singular")
endif(WITH_MUMPS)

endif(REAL)

if(COMPLEX)
set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-cplx)

if(WITH_UMFPACK)
	add_test(umfpack-solver-cplx-1 sh -c "${BIN} umfpack | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
	add_test(umfpack-solver-cplx-b-1 sh -c "${BIN} umfpack-block | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
endif(WITH_UMFPACK)

if(WITH_TRILINOS)
	# TODO: if we do not have aztecoo, do not run these tests
	add_test(aztecoo-solver-cplx-1 sh -c "${BIN} aztecoo | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
	add_test(aztecoo-solver-cplx-b-1 sh -c "${BIN} aztecoo-block | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
endif(WITH_TRILINOS)

if(WITH_MUMPS)
	add_test(mumps-solver-cplx-1 sh -c "${BIN} mumps | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
	add_test(mumps-solver-cplx-b-1 sh -c "${BIN} mumps-block | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
endif(WITH_MUMPS)

endif(COMPLEX)
