#include <fstream>
#include <inttypes.h>

#include "common.h"
#include "io.h"

FastaReader::FastaReader(const char *filename)
{
  file_obj_m = new std::ifstream(filename);

  file_obj_m->seekg(0, file_obj_m->end);
  file_length_m = file_obj_m->tellg();
  file_obj_m->seekg(0, file_obj_m->beg);
}

FastaReader::~FastaReader()
{
  file_obj_m->close();
  delete file_obj_m;
}

int FastaReader::read_chunk(char **buffer, 
    const int num_reads, const int header_length, const int read_length)
{
  *buffer = new char[num_reads*read_length];
  int reads = 0;
  
  while (true)
  {
    if (bytes_read_m >= file_length_m)
    {
      done_m = true;
      break;
    }
    if (reads >= num_reads)
    {
      break;
    }

    // Skip header
    file_obj_m->seekg(header_length + 1, std::ios::cur);
    bytes_read_m += header_length + 1;

    // Read DNA read
    file_obj_m->read(&(*buffer)[reads*read_length], read_length);
    file_obj_m->seekg(1, std::ios::cur);
    bytes_read_m += read_length + 1;

    reads++;
  }

  return reads;
}
