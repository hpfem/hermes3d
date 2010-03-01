project(hermes3d)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})		# for config.h

# source files for the hermes3d library
set(SRC
	hermes3d.cc
	filter.cc
	forms.cc
	function.cc
	mesh.cc
	feproblem.cc
	linproblem.cc
	loader/exodusii.cc
	loader/mesh3d.cc
	loader/hdf5.cc
	matrix.cc
	norm.cc
	output/gmsh.cc
	output/vtk.cc
	output/graph.cc
	precond/ifpack.cc
	precond/ml.cc
	qsort.cc
	quadcheb.cc
	quadstd.cc
	refdomain.cc
	refmap.cc
	shapefn.cc
	shapeset.cc
	shapeset/lobatto.cc
	shapeset/h1lobattotetra.cc
	shapeset/h1lobattotetradx.cc
	shapeset/h1lobattotetrady.cc
	shapeset/h1lobattotetradz.cc
	shapeset/h1lobattohex.cc
	shapeset/hcurllobattohex.cc
	shapeset/refmapss.cc
	solution.cc
	space.cc
	space/h1.cc
	space/hcurl.cc
	transform.cc
	traverse.cc
	weakform.cc
	solver/pardiso.cc
	solver/petsc.cc
	solver/umfpack.cc
	solver/epetra.cc
	solver/amesos.cc
	solver/aztecoo.cc
	solver/nox.cc
	solver/mumps.cc
	adapt/h1.cc
	adapt/proj.cc
	adapt/projipol.cc
	adapt/h1proj.cc
	adapt/h1projipol.cc
	${HERMES_COMMON_DIR}/callstack.cc
	${HERMES_COMMON_DIR}/trace.cc
	${HERMES_COMMON_DIR}/error.cc
	${HERMES_COMMON_DIR}/utils.cc
	${HERMES_COMMON_DIR}/timer.cc
)

macro(LIBRARY_SETTINGS BIN)
	target_link_libraries(${BIN} ${JUDY_LIBRARY})
	target_link_libraries(${BIN} ${BLAS_LIBRARIES})
	if(WITH_PETSC)
		include_directories(${PETSC_INCLUDE_DIR})
		target_link_libraries(${BIN} ${PETSC_LIBRARIES})
		target_link_libraries(${BIN} ${PETSC_LIBS})
		target_link_libraries(${BIN} ${LAPACK_LIBRARIES})
	endif(WITH_PETSC)

	if(WITH_UMFPACK)
		include_directories(${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})
		target_link_libraries(${BIN} ${UMFPACK_LIBRARY} ${AMD_LIBRARY})
		target_link_libraries(${BIN} ${LAPACK_LIBRARIES})
	endif(WITH_UMFPACK)

	if(WITH_PARDISO)
		target_link_libraries(${BIN} ${PARDISO_LIBRARY})
		target_link_libraries(${BIN} ${LAPACK_LIBRARIES})
	endif(WITH_PARDISO)

	if(WITH_MUMPS)
		include_directories(${MUMPS_INCLUDE_DIR})
		target_link_libraries(${BIN} ${MUMPS_LIBRARIES})
		target_link_libraries(${BIN} ${LAPACK_LIBRARIES})
	endif(WITH_MUMPS)

	if(WITH_TRILINOS)
		include_directories(${TRILINOS_INCLUDE_DIR})
		target_link_libraries(${BIN} ${TRILINOS_LIBRARIES})
	endif(WITH_TRILINOS)

	if(WITH_HDF5)
	endif(WITH_HDF5)

	if(WITH_EXODUSII)
		include_directories(${EXODUSII_INCLUDE_DIR})
		target_link_libraries(${BIN} ${EXODUSII_LIBRARIES})
	endif(WITH_EXODUSII)

	target_link_libraries(${BIN} ${ADDITIONAL_LIBS})
endmacro(LIBRARY_SETTINGS)

# real version
if(REAL)
	add_library(${HERMES_REAL_BIN} SHARED ${SRC})
	LIBRARY_SETTINGS(${HERMES_REAL_BIN})
endif(REAL)

# complex version
if(COMPLEX)
	add_library(${HERMES_CPLX_BIN} SHARED ${SRC})
	set_target_properties(${HERMES_CPLX_BIN} PROPERTIES COMPILE_FLAGS -DCOMPLEX)
	LIBRARY_SETTINGS(${HERMES_CPLX_BIN})
endif(COMPLEX)

#
configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/h3dconfig.h
)


#
# INSTALL
#

# install header files
file(GLOB INC           "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
file(GLOB INC_COMMON    "${HERMES_COMMON_DIR}/*.h")
file(GLOB INC_SPACE     "${CMAKE_CURRENT_SOURCE_DIR}/space/*.h")
file(GLOB INC_SHAPESET  "${CMAKE_CURRENT_SOURCE_DIR}/shapeset/*.h")
file(GLOB INC_LOADER    "${CMAKE_CURRENT_SOURCE_DIR}/loader/*.h")
file(GLOB INC_SOLVER    "${CMAKE_CURRENT_SOURCE_DIR}/solver/*.h")
file(GLOB INC_PRECOND   "${CMAKE_CURRENT_SOURCE_DIR}/precond/*.h")
file(GLOB INC_OUTPUT    "${CMAKE_CURRENT_SOURCE_DIR}/output/*.h")
file(GLOB INC_INTEGRALS "${CMAKE_CURRENT_SOURCE_DIR}/integrals/*.h")
file(GLOB INC_ADAPT     "${CMAKE_CURRENT_SOURCE_DIR}/adapt/*.h")

install(FILES ${INC}           DESTINATION include/)
install(FILES ${INC_COMMON}    DESTINATION include/common)
install(FILES ${INC_SPACE}     DESTINATION include/space)
install(FILES ${INC_SHAPESET}  DESTINATION include/shapeset)
install(FILES ${INC_LOADER}    DESTINATION include/loader)
install(FILES ${INC_SOLVER}    DESTINATION include/solver)
install(FILES ${INC_PRECOND}   DESTINATION include/precond)
install(FILES ${INC_OUTPUT}    DESTINATION include/output)
install(FILES ${INC_INTEGRALS} DESTINATION include/integrals)
install(FILES ${INC_ADAPT}     DESTINATION include/adapt)

# install library
if(REAL)
	install(TARGETS ${HERMES_REAL_BIN} LIBRARY DESTINATION lib)
endif(REAL)

if(COMPLEX)
	install(TARGETS ${HERMES_CPLX_BIN} LIBRARY DESTINATION lib)
endif(COMPLEX)
