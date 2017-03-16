CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(AlphaMate)
ENABLE_LANGUAGE(Fortran)

OPTION(RUN_TESTS "Run and compile tests" ON)
OPTION(USE_HDF5 "Use hdf" OFF)

SET(CLUSTER 0 CACHE STRING "cluster size")
# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../AlphaCMakeFiles/cmake/modules")
SET(CMAKE_MACOSX_RPATH 0)

INCLUDE(${CMAKE_MODULE_PATH}/setVersion.cmake)

ADD_DEFINITIONS(-DCLUSTER=${CLUSTER})

# Define the executable name
SET(ALPHAMATEEXE AlphaMate)
MESSAGE("Prog name ${ALPHAMATEEXE}")
# Define the library name
SET(AHLIB AlphaHouse)
SET(ALPHARELATELIB AlphaRelate)

# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(OBJ ${CMAKE_SOURCE_DIR}/objs)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)
SET(TESTS ${CMAKE_SOURCE_DIR}/tests)
SET(SRCALPHAMATE ${SRC}/)

# Defined library source

# Uncomment if it is required that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
  MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)

IF(${RUN_TESTS})
  SET(testDeps ${SRC}/AlphaMateModule.f90)
ENDIF(${RUN_TESTS})

# set(_files "test.pf")

INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

INCLUDE(${CMAKE_MODULE_PATH}/findMKL.cmake)
# Setup the LAPACK libraries.  This also takes care of peculiarities, such as
# the fact the searching for MKL requires a C compiler, and that the results
# are not stored in the cache.
#INCLUDE(${CMAKE_MODULE_PATH}/SetUpLAPACK.cmake)
IF(${USE_HDF5})
  INCLUDE(${CMAKE_MODULE_PATH}/findHDF5.cmake)
ENDIF(${USE_HDF5})

# Doxygen
FIND_PACKAGE(Doxygen)
INCLUDE(${CMAKE_MODULE_PATH}/findDoxygen.cmake)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)
SET(SRCAH ../AlphaHouse/src/)

# Have the .mod files placed in the OBJ folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${OBJ})

ADD_SUBDIRECTORY(${SRCAH} ${OBJ}/AH)

# TODO add tests to correct location
ADD_SUBDIRECTORY(${SRCALPHAMATE} ${BIN})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_MODULE_PATH}/distclean.cmake
)

IF(${RUN_TESTS})
  INCLUDE(${CMAKE_MODULE_PATH}/findPFUnit.cmake)
  TARGET_LINK_LIBRARIES(pftest_alltests ${ALPHARELATELIB})
ENDIF(${RUN_TESTS})

