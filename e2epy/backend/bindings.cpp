#include <inttypes.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "hashtable.h"

namespace py = pybind11;

PYBIND11_MODULE(e2epy_backend, m) 
{
  m.doc() = "Documentation for the e2epy backend module";

  py::class_<HashTable>(m, "HashTable")
    .def(py::init([](py::array_t<uint64_t> &np_keys, const int capacity)
    {
      const uint64_t *keys = np_keys.data(); 
      const int size = np_keys.size();
      const bool keys_on_device = false;
      return new HashTable(keys, keys_on_device, size, capacity);
    }))
    .def(py::init([](long keys_ptr, const int size, const int capacity)
    {
      const uint64_t *keys = reinterpret_cast<uint64_t *>(keys_ptr);
      const bool keys_on_device = true;
      return new HashTable(keys, keys_on_device, size, capacity);
    }))
    .def("count", [](HashTable &self, 
          const std::string &filename, const int header_length, const int read_length, 
          const int reads_per_chunk, const int kmer_size)
    {
      self.count(filename.c_str(), header_length, read_length, reads_per_chunk, kmer_size);
    })
    .def("count_fasta_chunks", [](HashTable &self, 
          const std::string &filename, const int chunk_size, 
          const int header_length, const int read_length, const int kmer_size)
    {
      self.count_fasta_chunks(
          filename.c_str(), chunk_size, header_length, read_length, kmer_size);
    })
    .def("lookup", [](HashTable &self, py::array_t<uint64_t> &np_keys)
    {
      py::buffer_info buf = np_keys.request();
      const uint64_t *keys = np_keys.data();
      const int size = np_keys.size();

      auto ret = py::array_t<uint32_t>(buf.size);
      uint32_t *values = ret.mutable_data();

      self.lookup(keys, values, size);
      return ret;
    })
    ;
}
