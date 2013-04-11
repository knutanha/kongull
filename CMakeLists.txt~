PROJECT(poisson)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# comment if you only going to write fortran code
ENABLE_LANGUAGE(C)

# needed no matter since the fst is written in fortran
ENABLE_LANGUAGE(Fortran)

#Mine kommandoer:
OPTION(USE_OPENMP "Use OpenMP?" ON)
OPTION(USE_MPI "Use MPI?" ON)

IF(USE_OPENMP)
  FIND_PACKAGE(OpenMP)
ENDIF(USE_OPENMP)

IF(USE_MPI)
	FIND_PACKAGE(MPI)
ENDIF(USE_MPI)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

IF(OPENMP_FOUND)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_OPENMP=1")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ENDIF(OPENMP_FOUND)

IF(MPI_FOUND)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_MPI=1")
  SET(INCLUDES ${INCLUDES} ${MPI_INCLUDE_PATH})
  SET(DEPLIBS ${DEPLIBS} ${MPI_C_LIBRARIES})
  MESSAGE("NÅ SKAL MPI VÆRE MED!")
ENDIF(MPI_FOUND)


# Send the include directories to the compiler
INCLUDE_DIRECTORIES(${INCLUDES})


#ADD_EXECUTABLE(poisson poisson.c fst.f)
#SET_TARGET_PROPERTIES(poisson PROPERTIES LINKER_LANGUAGE C)

#ADD_EXECUTABLE(poisson_f poisson.f fst.f)

#ADD_EXECUTABLE(mat matrisesjekk.c)

#ADD_EXECUTABLE(alle alltoallcheck.c)


# Jeg har kalt fila jeg kjører fra for poisson_mp.c
ADD_EXECUTABLE(testlol poisson_paral.c fst.f)
SET_TARGET_PROPERTIES(testlol PROPERTIES LINKER_LANGUAGE C)
ADD_EXECUTABLE(poisson_omp poisson_OMP.c fst.f)
SET_TARGET_PROPERTIES(poisson_omp PROPERTIES LINKER_LANGUAGE C)
ADD_EXECUTABLE(poisson_time poisson_time.c fst.f)
SET_TARGET_PROPERTIES(poisson_time PROPERTIES LINKER_LANGUAGE C)
ADD_EXECUTABLE(poisson_conv convergence.c fst.f)
SET_TARGET_PROPERTIES(poisson_conv PROPERTIES LINKER_LANGUAGE C)
ADD_EXECUTABLE(conv_time conv_time.c fst.f)
SET_TARGET_PROPERTIES(conv_time PROPERTIES LINKER_LANGUAGE C)


SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)
