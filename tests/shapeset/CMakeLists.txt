#
# Testing shapesets
#

if(DEV_TESTS)
	if(REAL)
		# test cases with REAL version od Hermes3D
		add_subdirectory(hex-h1)
		#add_subdirectory(hex-h1-lobatto-ced)
		add_subdirectory(hex-hcurl)
	endif(REAL)
endif(DEV_TESTS)