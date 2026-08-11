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

CMAKE_MINIMUM_REQUIRED(VERSION 3.5.0)

find_package(deal.II ${DEAL_II_VERSION_REQUIRED}
    HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../../../ $ENV{DEAL_II_DIR} $EBROOTDEALII
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

# targets libopenifem.so file file from build and imports std::experimental::filesystem
target_link_libraries(${TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/../OpenIFEM/OpenIFEM-build/lib/libopenifem.so stdc++fs)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(${TARGET} ${MPI_LIBRARIES})

# Create list of all parameter and shell files found in source directory
# FILES = list of all files in source directory (with directory path)
# FILE = currently selected file (either .prm or .sh) file with directory (iterating through FILES)
# FILE_NAME = only file name, removed directory location, which is needed for copying to build directory
file(GLOB FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.prm ${CMAKE_CURRENT_SOURCE_DIR}/*.sh
)

# Colours being used for file management/updating from source directory
if(NOT WIN32)
    string(ASCII 27 Esc)
    set(ColourReset "${Esc}[m")
    set(Green       "${Esc}[32m")
    set(Yellow         "${Esc}[33m")
endif()

# Iterate through each file found in build directory
foreach(FILE ${FILES})
    get_filename_component(FILE_NAME ${FILE} NAME)
    # Copy the current file from source directory to build if it does not already exist
    # Files are not overwritten if currently existing to avoid deleting local changes per simulation 
    if(NOT(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${FILE_NAME}))
        file(COPY ${FILE_NAME}
                DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
        )
        MESSAGE(STATUS "${Yellow}COPIED " ${FILE_NAME}
                    "${ColourReset} from " ${CMAKE_CURRENT_SOURCE_DIR})
    else()
        MESSAGE(STATUS "${Green}KEPT " ${FILE_NAME}
                    "${ColourReset} from " ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
endforeach()


# Create list of all mesh files from the source directory
# Same process as other files, just for meshes located in the /meshes folder
# MESHES = list of all meshes in source directory (with directory path)
# MESH = currently selected mesh with directory (iterating through MESHES)
# MESH_NAME = only mesh name, removed directory location, which is needed for copying to build directory
file(GLOB MESHES
    ${CMAKE_CURRENT_SOURCE_DIR}/meshes/*.msh
)

foreach(MESH ${MESHES})
    get_filename_component(MESH_NAME ${MESH} NAME)
    if(NOT(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/meshes/${MESH_NAME}))
        file(COPY ${MESH}
                DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
        MESSAGE(STATUS "${Yellow}COPIED " ${MESH_NAME}
                    "${ColourReset} from " ${CMAKE_CURRENT_SOURCE_DIR} "/meshes")
    else()
        MESSAGE(STATUS "${Green}KEPT " ${MESH_NAME}
                    "${ColourReset} from " ${CMAKE_CURRENT_SOURCE_DIR} "/meshes")
    endif()
endforeach()

# include(CTest)
# enable_testing()