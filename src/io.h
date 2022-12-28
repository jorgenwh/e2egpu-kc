#ifndef IO_H_
#define IO_H_

#include <fstream>
#include <inttypes.h>

#include "common.h"

class FastaReader
{
public:
  FastaReader(const char *filename);
  ~FastaReader();

  int read_chunk(char **buffer, const int num_reads);
  bool done() const { return done_m; }
private:
  std::ifstream *file_obj_m;
  uint64_t file_length_m;
  uint64_t bytes_read_m = 0;
  bool done_m = false;
};

#endif // IO_H_
