#Variable declaration to be used later
#variables declared as "SET(NAME value)" and called as "${NAME}"
SET(TARGET "main")
#setting name of target cpp file
SET(TARGET_SRC
    ${TARGET}.cpp
)

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

# message(STATUS $ENV{DEAL_II_DIR})

IF(NOT ${deal.II_FOUND})
    MESSAGE(FATAL_ERROR"\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET} VERSION 1.0.0)
DEAL_II_INVOKE_AUTOPILOT()

# these 2 lines are very important in getting OpenIFEM to work properly
# targets header files in the include folder
target_include_directories(${TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../OpenIFEM/OpenIFEM-source/OpenIFEM/include)
# target_include_directories(${TARGET} PUBLIC "/home/nhewko/OpenIFEM/OpenIFEM-source/OpenIFEM/include")

# deal_ii_setup_target(${TARGET})
# targets libopenifem.so file file from build and imports std::experimental::filesystem
target_link_libraries(${TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/../OpenIFEM/OpenIFEM-build/lib/libopenifem.so stdc++fs)
# target_link_libraries(${TARGET} "/home/nhewko/OpenIFEM/OpenIFEM-build/lib/libopenifem.so" stdc++fs)

# include(CTest)
# enable_testing()
