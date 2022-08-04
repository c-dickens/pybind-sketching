#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Need this for automatic type conversion from std::vector back to python lists.
#include "count_min_sketch.h"

namespace py = pybind11;

PYBIND11_MODULE(count_min_sketch, m) {

    py::class_<CountMinSketch>(m, "CountMinSketch")
        .def(py::init<uint64_t, uint64_t, uint64_t>(), "Initializes the CountMin sketch data structure")
        .def("get_num_hashes", &CountMinSketch::get_num_hashes,
        "Returns the number of hash functions used in the sketch.")
        .def("get_num_buckets", &CountMinSketch::get_num_buckets,
        "Returns the number of buckets the hash functions map into.")
        .def("get_seed", &CountMinSketch::get_seed,
        "Returns the seed used for hashing via the PRNG.  Necessary for the merge step.")
        .def("get_table_shape", &CountMinSketch::get_table_shape,
        "Returns the shape of the sketch table as a python tuple.")
        .def("get_config", &CountMinSketch::get_config,
        "Returns the sketch configuration as a python list")
        .def("get_table", &CountMinSketch::get_table,
        "Returns a copy of the sketch table as a list of python lists.") // Returning np arrays is harder so keep list.
        .def("get_total_weight", &CountMinSketch::get_total_weight,
        "Returns the total mass added to the stream.")
        .def("update", (void (CountMinSketch::*)(uint64_t, int64_t)) &CountMinSketch::update, py::arg("item"), py::arg("weight")=1,
        "Updates the sketch with the item")
        .def("update", (void (CountMinSketch::*)(const std::string&, int64_t)) &CountMinSketch::update, py::arg("item"), py::arg("weight")=1,
        "Updates the sketch with the item")
        .def("get_estimate", (int64_t (CountMinSketch::*)(uint64_t)) &CountMinSketch::get_estimate, py::arg("item"),
        "Returns the frequency estimate for the given item.")
        .def("get_estimate", (int64_t (CountMinSketch::*)(const std::string&)) &CountMinSketch::get_estimate, py::arg("item"),
        "Returns the frequency estimate for the given item.")
        .def("get_upper_bound", &CountMinSketch::get_upper_bound, py::arg("item"),
        "Returns the frequency estimate upper bound for the given item.")
        .def("get_lower_bound", &CountMinSketch::get_upper_bound, py::arg("item"),
        "Returns the frequency estimate lower bound for the given item.")
        .def_static("suggest_num_buckets", &CountMinSketch::suggest_num_buckets, py::arg("relative_error"),
        "Returns the smallest number of buckets that will guarantee the given relative_error.")
        .def_static("suggest_num_hashes", &CountMinSketch::suggest_num_hashes, py::arg("confidence"),
        "Returns the smallest number of hash functions that will guarantee the given confidence (1-failure probability.")
        .def("get_relative_error", &CountMinSketch::get_relative_error,
        "Returns the error upper bound of the sketch estimates.")
        .def("get_confidence", &CountMinSketch::get_confidence,
        "Returns the confidence (between 0.0 and 1.0) with which the estimates are returned.")
        .def("merge", &CountMinSketch::merge, py::arg("other_count_min"),
        "Returns other_count_min sketch into the current CountMinSketch object.")
        .def("__add__", &CountMinSketch::merge, py::arg("other_count_min"),
        "Returns other_count_min sketch into the current CountMinSketch object."\
                "Note that this operation acts in place!")
        .def("to_string", &CountMinSketch::to_string,
        "Produces a string summary of the sketch")
        .def("__str__", &CountMinSketch::to_string,
        "Produces a string summary of the sketch");
}