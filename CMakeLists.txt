cmake_minimum_required(VERSION 3.16)
project(e2egpu_kc LANGUAGES CXX CUDA)

set(CMAKE_CXX_FLAGS "-march=native -O3")
set(CMAKE_CUDA_FLAGS "-march=native -O3")
#set(CMAKE_CXX_FLAGS "-march=native -g -Wall")
#set(CMAKE_CUDA_FLAGS "-march=native -g -Wall")
set(CMAKE_CUDA_ARCHITECTURES "75")

set(SOURCES ${SOURCES}
  src/io.cpp
  src/r2kmer.cpp
)

add_executable(${PROJECT_NAME} src/main.cpp ${SOURCES})
