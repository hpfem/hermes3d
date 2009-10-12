#
# CMake for test the calculations
#

# TETRA
if(WITH_TETRA)
	if(REAL)
		add_subdirectory(tetra-h1)
		add_subdirectory(tetra-h1-dirichlet)
		add_subdirectory(tetra-h1-neumann)
		add_subdirectory(tetra-h1-newton)
	endif(REAL)
endif(WITH_TETRA)

# HEX
if(WITH_HEX)
    # H1
	if(REAL)
		add_subdirectory(hex-h1)
		add_subdirectory(hex-h1-1face)
		add_subdirectory(hex-h1-dirichlet)
		add_subdirectory(hex-h1-neumann)
		add_subdirectory(hex-h1-newton)
		add_subdirectory(hex-h1-unsym)
		# systems of equations
		add_subdirectory(hex-h1-sys)
		add_subdirectory(hex-h1-sys-dirichlet)
		add_subdirectory(hex-h1-sys-coupled)
		# non linear
		if(HAVE_NOX)
			add_subdirectory(hex-h1-nonlin)
			add_subdirectory(hex-h1-matfree)
		endif(HAVE_NOX)
	endif(REAL)
	# Hcurl
	if(REAL)
		add_subdirectory(hex-hcurl)		
	endif(REAL)
	if(COMPLEX)
		add_subdirectory(hex-hcurl-cplx)
		add_subdirectory(hex-hcurl-imped)
	endif(COMPLEX)
endif(WITH_HEX)
