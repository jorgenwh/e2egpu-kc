cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_FLAGS "-march=native -O3")
set(CMAKE_CUDA_FLAGS "-march=native -O3")
set(CMAKE_CUDA_ARCHITECTURES "75")

set(SOURCES ${SOURCES}
  src/io.cpp
  src/r2kmer.cpp
  src/kernels.cu
  src/hashtable.cu
)

add_compile_definitions(__CUDA_ERROR_CHECK__)
add_compile_definitions(__USE_MURMUR_HASH__)

find_package(CUDA REQUIRED)

if(${py} MATCHES "True")
  project(e2epy_backend LANGUAGES CXX CUDA)
  find_package(pybind11 REQUIRED)
  set(SOURCE_DIR "e2epy/backend")
  include_directories(${SOURCE_DIR})
  pybind11_add_module(${PROJECT_NAME} ${SOURCES} "${SOURCE_DIR}/bindings.cpp")
else()
  project(e2egpu_kc LANGUAGES CXX CUDA)
  add_executable(${PROJECT_NAME} src/main.cpp ${SOURCES})
endif()

