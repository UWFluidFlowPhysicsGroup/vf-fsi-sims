#Variable declaration to be used later
#variables declared as "SET(NAME value)" and called as "${NAME}"
SET(TARGET "main")
#setting name of target cpp file
SET(TARGET_SRC
    ${TARGET}.cpp
)
enable_language(Fortran)

#setting minimum versions of programs from OpenIFEM CMakeLists
SET(CMAKE_CXX_STANDARD 17)
SET(DEAL_II_VERSION_REQUIRED 9.3.0)
#TODO find proper file version from package info, not needed
#SET(OPENIFEM_VERSION_REQUIRED 1.0.0)

CMAKE_MINIMUM_REQUIRED(VERSION 3.5.0)

#find deal.II package, mostly imported from examples
#deal.II_DIR is at /home/nhewko/local/lib/cmake/deal.II

# Might need to change hints for server simulations,
# server dealii found at $EBROOTDEALII
find_package(deal.II ${DEAL_II_VERSION_REQUIRED}
    # HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../../../ $ENV{DEAL_II_DIR}
)

IF(NOT ${deal.II_FOUND})
    MESSAGE(FATAL_ERROR"\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

# TODO add checks for dependencies?
if(NOT DEAL_II_WITH_UMFPACK) # keep in one line
  message(FATAL_ERROR "
Error! This tutorial requires a deal.II library that was configured with the following options:
    DEAL_II_WITH_UMFPACK = ON
However, the deal.II library found at ${DEAL_II_PATH} was configured with these options
    DEAL_II_WITH_UMFPACK = ${DEAL_II_WITH_UMFPACK}
which conflict with the requirements."
    )
endif()

if(NOT (DEAL_II_WITH_PETSC OR DEAL_II_WITH_TRILINOS) OR NOT DEAL_II_WITH_P4EST OR DEAL_II_PETSC_WITH_COMPLEX) # keep in one line
  message(FATAL_ERROR "
Error! This tutorial requires a deal.II library that was configured with the following options:
    DEAL_II_WITH_PETSC = ON
    DEAL_II_PETSC_WITH_COMPLEX = OFF
    DEAL_II_WITH_P4EST = ON
or
    DEAL_II_WITH_TRILINOS = ON
    DEAL_II_WITH_P4EST = ON
However, the deal.II library found at ${DEAL_II_PATH} was configured with these options
    DEAL_II_WITH_PETSC = ${DEAL_II_WITH_PETSC}
    DEAL_II_PETSC_WITH_COMPLEX = ${DEAL_II_PETSC_WITH_COMPLEX}
    DEAL_II_WITH_P4EST = ${DEAL_II_WITH_P4EST}
    DEAL_II_WITH_TRILINOS = ${DEAL_II_WITH_TRILINOS}
which conflict with the requirements.
One or both of the aforementioned combinations of prerequisites are not met by your installation, but at least one is required for this tutorial step."
    )
endif()

if(NOT DEAL_II_WITH_METIS)
  message(FATAL_ERROR "
Error! OpenIFEM requires a deal.ii library that was configured with the following options:
    DEAL_II_WITH_METIS = ON
However, the deal.II library found at ${DEAL_II_PATH} was configured with these options
    DEAL_II_WITH_METIS = ${DEAL_II_WITH_METIS}
which conflict with the requirements.
The prerequisite are not met by your installation, but this is required for OpenIFEM."
   )
endif()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET} VERSION 1.0.0)
DEAL_II_INVOKE_AUTOPILOT()

# these 2 lines are very important in getting OpenIFEM to work properly
# targets header files in the include folder
target_include_directories(${TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../OpenIFEM/OpenIFEM-source/OpenIFEM/include)

# deal_ii_setup_target(${TARGET})
# targets libopenifem.so file file from build and imports std::experimental::filesystem
#target_link_libraries(${TARGET} "/home/nhewko/OpenIFEM/OpenIFEM-build/lib/libopenifem.so" stdc++fs)
target_link_libraries(${TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/../OpenIFEM/OpenIFEM-build/lib/libopenifem.so stdc++fs)

# use step 7 here
# https://p4est.github.io/api/p4est-latest/installing_p4est.html

find_package(MPI REQUIRED)
# message(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(${TARGET} ${MPI_LIBRARIES})

target_include_directories(${TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../p4est_build/local/include)
target_link_libraries(${TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/../p4est_build/local/lib/libp4est.so stdc++fs)
# target_include_directories(${TARGET} PUBLIC "/home/nhewko/p4est_build/local/include")
# target_link_libraries(${TARGET} "/home/nhewko/p4est_build/local/lib/libp4est.so" stdc++fs)

# include(CTest)
# enable_testing()
