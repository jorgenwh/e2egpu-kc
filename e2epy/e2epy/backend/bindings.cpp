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
    ;
}
