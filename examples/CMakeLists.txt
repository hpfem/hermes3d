#
# CMake for examples
#

include (../tests/common.cmake)

# HEX
if(WITH_HEX)
    # H1
	if(REAL)
		add_subdirectory(elastostatics)
		add_subdirectory(fichera)
		add_subdirectory(heat-conduction)
		add_subdirectory(heat-cubit)
		add_subdirectory(layer)
		add_subdirectory(singpert-aniso)
	endif(REAL)
	if(COMPLEX)
		add_subdirectory(bessel)
	endif(COMPLEX)
endif(WITH_HEX)
