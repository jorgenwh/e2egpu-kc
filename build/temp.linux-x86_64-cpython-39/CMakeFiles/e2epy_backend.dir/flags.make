# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# compile CUDA with /usr/local/cuda-11.7/bin/nvcc
# compile CXX with /usr/bin/c++
CUDA_DEFINES = -D__CUDA_ERROR_CHECK__ -D__USE_MURMUR_HASH__ -De2epy_backend_EXPORTS

CUDA_INCLUDES = --options-file CMakeFiles/e2epy_backend.dir/includes_CUDA.rsp

CUDA_FLAGS = -march=native -O3 -O3 -DNDEBUG --generate-code=arch=compute_75,code=[compute_75,sm_75] -Xcompiler=-fPIC -Xcompiler=-fvisibility=hidden -std=c++14

CXX_DEFINES = -D__CUDA_ERROR_CHECK__ -D__USE_MURMUR_HASH__ -De2epy_backend_EXPORTS

CXX_INCLUDES = -I/home/jorgen/projects/e2egpu-kc/e2epy/backend -isystem /tmp/pip-build-env-a6uc74kr/overlay/include -isystem /home/jorgen/miniconda3/include/python3.9

CXX_FLAGS = -march=native -O3 -O3 -DNDEBUG -fPIC -fvisibility=hidden -flto -fno-fat-lto-objects
