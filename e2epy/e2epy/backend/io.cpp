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

int FastaReader::read_chunk(char **buffer, const int num_reads)
{
  *buffer = new char[num_reads*READ_LENGHT];
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
    file_obj_m->seekg(HEADER_LENGTH + 1, std::ios::cur);
    bytes_read_m += HEADER_LENGTH + 1;

    // Read DNA read
    file_obj_m->read(&(*buffer)[reads*READ_LENGHT], READ_LENGHT);
    file_obj_m->seekg(1, std::ios::cur);
    bytes_read_m += READ_LENGHT + 1;

    reads++;
  }

  return reads;
}
