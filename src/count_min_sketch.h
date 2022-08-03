//
// Created by Charlie Dickens on 03/08/2022.
//

#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Need this for automatic type conversion from std::vector back to python lists.
#include "pybind11/numpy.h"
#include <vector>
#include <cmath>
#include <random>
#include <iostream>

#ifndef SRC_COUNT_MIN_SKETCH_H
#define SRC_COUNT_MIN_SKETCH_H


class CountMinSketch{
public:
    uint64_t num_hashes, num_buckets, seed ;
    CountMinSketch(uint64_t num_hashes, uint64_t num_buckets, uint64_t seed) ;
    const uint64_t get_num_hashes() const { return num_hashes; }
    const uint64_t get_num_buckets() const { return num_buckets; }
    const uint64_t get_seed() const { return seed; } // nb will need this for merging.
    std::pair<uint64_t, uint64_t> get_table_shape() const {return {get_num_hashes(), get_num_buckets()} ; } ;
    std::vector<uint64_t> get_config() ;
    std::vector<std::vector<int64_t>> get_table() ;
    int64_t get_total_weight() {return total_weight ; }
    void update(int64_t item, int64_t weight=1) ;
    int64_t get_estimate(uint64_t item) ;
    static uint64_t suggest_num_buckets(float relative_error) ;
    static uint64_t suggest_num_hashes(float confidence) ;
    int64_t get_upper_bound(uint64_t item) ;
    int64_t get_lower_bound(uint64_t item) ;
    const float get_relative_error() ;
    const float get_confidence() ;
    void merge(CountMinSketch &sketch) ;
    std::string to_string() ;




private:
    float epsilon ; // Error parameter
    int64_t total_weight = 0 ; // This tracks how much weight has been added to the stream.
    std::vector<std::vector<int64_t>> table;
    uint64_t mersenne_exponent = 31 ; // NB move this to a private value
    uint64_t large_prime = (1 << mersenne_exponent) - 1 ; // nb change this to a mersenne prime
    std::vector<uint64_t> a_hash_params, b_hash_params ; // Keep these public for preliminary testing

    // Functions
    void set_hash_parameters() ;
    uint64_t get_bucket_hash(uint64_t item, uint64_t a, uint64_t b) ;
    int64_t get_size_estimate_in_bits() ;
};


#endif //SRC_COUNT_MIN_SKETCH_H
