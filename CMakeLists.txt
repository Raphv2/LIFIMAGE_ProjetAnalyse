# Global Project Configuration
cmake_minimum_required(VERSION 3.22)
project(Projet-Analyse)

# Find packages


find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

# Compiler options
set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_FLAGS_DEBUG "-ggdb -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")

# Set sources and includes
set(SOURCES
        # Classes
        # Other Sources
        # Libraries
        lib/gkit3/src/color.cpp
        lib/gkit3/src/files.cpp
        lib/gkit3/src/image_io.cpp
        lib/gkit3/src/mat.cpp
        lib/gkit3/src/mesh_io.cpp
        lib/gkit3/src/vec.cpp
)

set(INCLUDES
        include
        # Libraries
        lib/gkit3/src
)

# Executables
add_executable(Analyse src/analyse.cpp ${SOURCES})
target_include_directories(Analyse PUBLIC ${INCLUDES})
target_link_libraries(Analyse ${OpenCV_LIBS})

add_executable(Synthese src/synthese.cpp ${SOURCES})
target_include_directories(Synthese PUBLIC ${INCLUDES})

target_include_directories(Synthese PUBLIC ${OpenCV_INCLUDE_DIRS})
target_link_libraries(Synthese PUBLIC ${OpenCV_LIBS})

# OpenMP support
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Synthese PUBLIC OpenMP::OpenMP_CXX)
endif()

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)