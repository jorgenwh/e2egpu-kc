cmake_minimum_required(VERSION 3.16)
project(e2epy_backend LANGUAGES CXX CUDA)

find_package(CUDA REQUIRED)
find_package(pybind11 REQUIRED)

set(SOURCE_DIR "e2epy/backend")
include_directories(${SOURCE_DIR})

set(SOURCES ${SOURCES}
  "${SOURCE_DIR}/io.cpp"
  "${SOURCE_DIR}/kernels.cu"
  "${SOURCE_DIR}/hashtable.cu"
)

set(CMAKE_CXX_FLAGS "-march=native -O3")
set(CMAKE_CUDA_FLAGS "-march=native -O3")
set(CMAKE_CUDA_ARCHITECTURES "75")

if(${CUDA_ERROR_CHECK} MATCHES "True")
  add_compile_definitions(__CUDA_ERROR_CHECK__)
endif()
if(${USE_MURMUR_HASH} MATCHES "True")
  add_compile_definitions(__USE_MURMUR_HASH__)
endif()

pybind11_add_module(${PROJECT_NAME} ${SOURCES} "${SOURCE_DIR}/bindings.cpp")
