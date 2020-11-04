project(FLS)

cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_C_COMPILER /usr/local/bin/gcc) 
set(CMAKE_CXX_COMPILER /usr/local/bin/g++) 
set(CMAKE_CXX_FLAGS "-g -fopenmp -Ofast")
set(CMAKE_CXX_STANDARD 14) 

set(VTK_DIR "/home/sci/ssane/packages/VTK-build")
#set(VTK_DIR "/home/users/ssane/VTK-6.3.0")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(FLS main.cxx)

  target_link_libraries(FLS PRIVATE ${VTK_LIBRARIES} ${MPI_CXX_LIBRARIES} -fopenmp)
