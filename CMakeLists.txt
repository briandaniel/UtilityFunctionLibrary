
# This cmake script builds the library of some general functions used in other programs

# Standard cmake version requirement
cmake_minimum_required(VERSION 2.6)

# Inform user
MESSAGE("Configuring utility function library CPP.")

# Project (fortran must be enabled for OpenCMISS compatibility)
PROJECT( utilityFunctionLibrary LANGUAGES C CXX Fortran )

# Add all the local source files from the folder "src" to the makefile
file(GLOB_RECURSE UtilityFunctionCPP_SRC_LIBRARY
	"*.hpp"
	"*.cpp"
)

#---- external packages ----#
# 1. MPI
find_package( MPI REQUIRED  ) 
include_directories(${MPI_INCLUDE_PATH})


# Add the execubale using the source files in (src)
add_library(utilityFunctionLibrary SHARED
	${UtilityFunctionCPP_SRC_LIBRARY}
)

# Set version
set_target_properties(utilityFunctionLibrary PROPERTIES SOVERSION 1)

# Link the local and external libraries to the executable
target_link_libraries( 
	utilityFunctionLibrary
	mpi
	mpi_cxx
	${MPI_LIBRARIES}
	${VTK_LIBRARIES}
)

# Set the install location
install(TARGETS utilityFunctionLibrary 
			  LIBRARY DESTINATION lib
)
   
# Install the header files
INSTALL ( DIRECTORY ${CMAKE_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
  FILES_MATCHING PATTERN "*.hpp"
)

# Install the header files
INSTALL ( DIRECTORY ${CMAKE_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
  FILES_MATCHING PATTERN "*.h"
  )




































