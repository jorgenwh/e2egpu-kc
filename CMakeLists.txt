cmake_minimum_required(VERSION 3.16)
project(e2egpu_kc LANGUAGES CXX CUDA)

find_package(CUDA REQUIRED)

set(SOURCES ${SOURCES}
  src/io.cpp
  src/r2kmer.cpp
  src/kernels.cu
  src/hashtable.cu
)

set(CMAKE_CXX_FLAGS "-march=native -O3")
set(CMAKE_CUDA_FLAGS "-march=native -O3")
set(CMAKE_CUDA_ARCHITECTURES "75")

add_compile_definitions(__CUDA_ERROR_CHECK__)
add_compile_definitions(__USE_MURMUR_HASH__)

add_executable(${PROJECT_NAME} src/main.cpp ${SOURCES})
