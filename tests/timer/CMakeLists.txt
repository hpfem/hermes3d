project(timer)

include(CMake.vars OPTIONAL)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${hermes3d_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})

add_executable(${PROJECT_NAME}
	main.cc
	${HERMES_COMMON_DIR}/trace.cc
	${HERMES_COMMON_DIR}/timer.cc
)

target_link_libraries(${PROJECT_NAME} ${TESTS_LIBRARIES}) 

if(REAL)
	target_link_libraries(${PROJECT_NAME} ${HERMES_REAL_BIN})
elseif(COMPLEX)
	target_link_libraries(${PROJECT_NAME} ${HERMES_CPLX_BIN})
endif(REAL)

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Tests

add_test(timer-1 ${PROJECT_NAME})

