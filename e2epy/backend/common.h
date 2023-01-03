#ifndef COMMON_H_
#define COMMON_H_

#ifndef SDIV
  #define SDIV(x,y)(((x)+(y)-1)/(y))
#endif

#include <stdio.h>
#include <cuda_runtime.h>

#define cuda_errchk(err) { cuda_errcheck(err, __FILE__, __LINE__); }

inline void cuda_errcheck(cudaError_t code, const char *file, int line, bool abort=true) 
{
#ifdef __CUDA_ERROR_CHECK__
  if (code != cudaSuccess) 
  {
    switch (code) 
    {
      case 2:
        fprintf(stderr, "CUDA out of memory error in %s at line %d\n", file, line);
        break;
      default:
        fprintf(stderr, "CUDA assert: '%s', in %s, at line %d\n", cudaGetErrorString(code), file, line);
    }
    exit(code);
  }
#endif // __CUDA_ERROR_CHECK__
}

static const uint64_t kEmpty = 0xFFFFFFFFFFFFFFFF;

#endif // COMMON_H_
