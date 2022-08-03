//
// Created by Charlie Dickens on 03/08/2022.
//

#include "count_min_sketch.h"

CountMinSketch::CountMinSketch(uint64_t num_hashes, uint64_t num_buckets, uint64_t seed) :
        num_hashes(num_hashes), num_buckets(num_buckets), seed(seed) {
    table.resize(num_hashes, std::vector<int64_t>(num_buckets));
    set_hash_parameters() ;
} ;

void CountMinSketch::set_hash_parameters(){
    /* Sets the array containing a and b parameters for hashing.
     * a_hash_params contains all values of a for the hashing (see ::get_bucket_hash)
     * a_hash_params contains all values of b
     * We generate 2*num_hashes from CountingSketch::init_hash_parameters and then iterate through to get
     * the values of a and b for the respective arrays.
     * The line ``std::uniform_int_distribution<uint64_t> uniform_random_ints(1, large_prime - 1);``
     * returns uniform ints on [0, large_prime - 1] as per http://dimacs.rutgers.edu/~graham/pubs/papers/cmsoft.pdf
     * "Sketch Internals" page 3 but I think the [1, p-1] range is wrong as we need to operate over the integers
     * mod p.  Some details on this can be found in
     * Section 2.3 https://people.csail.mit.edu/ronitt/COURSE/S12/handouts/lec5.pdf
     *
     */
    std::default_random_engine rng(seed);
    std::uniform_int_distribution<uint64_t> uniform_random_ints(0, large_prime - 1);
    a_hash_params.reserve(num_hashes) ;
    b_hash_params.reserve(num_hashes) ;

    // By construction, every element in a_hash_params and b_hash_params is guaranteed to be in [0, large_prime - 1]
    for(uint64_t i=0 ; i < num_hashes ; ++i){
        //cout << uniform_random_ints(rng) << std::endl ;
        a_hash_params[i] = uniform_random_ints(rng);
        b_hash_params[i] = uniform_random_ints(rng);
    }
} ;  // End set_hash_parameters()

uint64_t CountMinSketch::get_bucket_hash(uint64_t item, uint64_t a, uint64_t b){
    /*
     * Performs bucket hashing, that is, for a given item and a given row in the sketch,
     * this function selects the bucket, meaning the column index of the sketch table.
     *
     * h = (a * x + b) % large_prime
     * h = h % nbuckets
     * This is executed once for each of the hash functions
     */
    //return uint64_t ((a * item + b) % large_prime) % num_buckets ;
    uint64_t hash = ((a * item + b) % large_prime) % num_buckets ;
    return (hash >= 0) ? hash : hash + large_prime ;
} ; // End get_bucket_hash()

std::vector<uint64_t>CountMinSketch::get_config() {
    std::vector<uint64_t> config(3) ;
    config = {get_num_hashes(), get_num_buckets(), get_seed() } ;
    return config ;
} ;

std::vector<std::vector<int64_t>> CountMinSketch::get_table() {
    /*
    * Returns a copy of the sketch.
    */
    std::vector<std::vector<int64_t>> sketch(num_hashes, std::vector<int64_t>(num_buckets));
    for(uint64_t i=0; i<num_hashes; i++){
        for(uint64_t j=0; j<num_buckets; j++){
            sketch[i][j] = table[i][j] ;
        }
    }
    return sketch;
} ;

void CountMinSketch::update(int64_t item, int64_t weight){
    /*
    * Updates the sketch with the item and a default weight of 1
    * iterates through the number of hash functions and gets the bucket index.
    * Then increment the sketch table at index (row, column) where row is one of the the hash
    * functions that are iterated over, and column is the corresponding bucket index.
    * Finally, increment the sketch table with the weight associated to the item.
    *
    * TODO: We can improve the update time by removing the modulus operations as described in
    * Section 3: http://dimacs.rutgers.edu/~graham/pubs/papers/cmsoft.pdf
    */
    if(item < 0) {
        throw std::invalid_argument( "Item identifier must be nonnegative." ) ;
    };

    for(uint64_t i=0; i < num_hashes; i++){
        uint64_t a = a_hash_params[i] ;
        uint64_t b = b_hash_params[i] ;
        uint64_t h = get_bucket_hash(item, a, b) ;
        table[i][h] += weight ;
    }
    total_weight += weight ;
} ; // End update function

int64_t CountMinSketch::get_estimate(uint64_t item) {
    /*
     * Returns the estimate from the sketch for the given item.
     * TODO:  Can we explore the estimator from this paper?
     * https://dl.acm.org/doi/10.1145/3219819.3219975
     */
    int64_t estimate = std::numeric_limits<int64_t>::max() ; // start arbitrarily large
    for(uint64_t i=0; i < num_hashes; i++){
        uint64_t a = a_hash_params[i] ;
        uint64_t b = b_hash_params[i] ;
        uint64_t h = get_bucket_hash(item, a, b) ;
        estimate = std::min(estimate, table[i][h]) ;
    }
    return estimate ;
} // end get_estimate()

int64_t CountMinSketch::get_upper_bound(uint64_t item) {
    /*
     * Returns the upper bound of the estimate as:
     * f_i - true frequency
     * est(f_i) - estimate frequency
     * f_i <= est(f_i)
     */
    return get_estimate(item) ;
} // End get_upper_bound()


int64_t CountMinSketch::get_lower_bound(uint64_t item) {
    /*
     * Returns the lower bound of the estimate as:
     * f_i - true frequency
     * est(f_i) - estimate frequency
     * f_i >= est(f_i) - epsilon*||f||_1 with ||f||_1 being the total weight in the sketch.
     */
    return get_estimate(item) - epsilon*total_weight ;
} // End get_lower_bound()

uint64_t CountMinSketch::suggest_num_hashes(float confidence){
    /*
     * Function to help users select a number of hashes for a given confidence
     * eg confidence is 1 - failure probability
     * failure probability == delta in the literature.
     * * TODO: Change this when update is improved
     */
    if(confidence < 0. || confidence > 1.0){
        throw std::invalid_argument( "Confidence must be between 0 and 1.0 (inclusive)." );
    }
    return ceil(log(1.0/(1.0 - confidence))) ;
} // End suggest_num_hashes

uint64_t CountMinSketch::suggest_num_buckets(float relative_error){
    /*
     * Function to help users select a number of buckets for a given error.
     * TODO: Change this when update is improved
     */
    if(relative_error < 0.){
        throw std::invalid_argument( "Relative error must be at least 0." );
    }
    return ceil(exp(1.0) / relative_error) ;
} // End suggest_num_buckets()

const float CountMinSketch::get_relative_error() {
    /*
     * Returns the relative error (epsilon) parameter when the sketch has been initialised
     */
    return exp(1.0) / num_buckets ;
} // End get_relative_error()

const float CountMinSketch::get_confidence() {
    /*
     * Returns the confiddenc= 1. - failure probability (confidence = 1 - delta) parameter when the sketch
     * has been initialised
     */
    return float(1.0 - exp(-(double)num_hashes));
} // End get_confidence()

void CountMinSketch::merge(CountMinSketch &sketch){
    /*
     * Merges this sketch into that sketch by elementwise summing of buckets
     */
    if(this == &sketch){
        throw std::invalid_argument( "Cannot merge a sketch with itself." );
    }

    bool same_sketch_config = (get_config() == sketch.get_config()) ;
    if(!same_sketch_config){
        throw std::invalid_argument( "Incompatible sketch config." );
    }

    // Iterate through the table and increment
    for(int i=0 ; i < num_hashes; i++){
        for(int j = -0; j < num_buckets; j++){
            table[i][j] += sketch.table[i][j] ;
        }
    }
    total_weight += sketch.total_weight ;
} // End merge()

std::string CountMinSketch::to_string(){
    /*
     * Converts the sketch information to a string format
     */
    std::ostringstream os;
    os          << "### CountMinSketch Summary:"  << std::endl ;
    os          << " Depth (number of hashes) : " << get_num_hashes()  << std::endl ;
    os          << " Width (number of buckets): " << get_num_buckets() << std::endl ;
    os          << " Relative Error (epsilon) : " << get_relative_error() << std::endl ;
    os          << " Confidence (1 - δ)       : " << get_confidence() << std::endl ;
    os          << " Total weight in sketch   : " << get_total_weight() << std::endl ;
    os          << " Seed                     : " << get_seed() << std::endl ;
    os          << " Size estimate (in bits)  : " <<  get_size_estimate_in_bits() << std::endl ;
    os          << "### End CountMinSketch Summary"  << std::endl ;
    return os.str() ;
} // End to_string()

int64_t CountMinSketch::get_size_estimate_in_bits(){
    /* Estimates the size in bits of the sketch
     * Calculation
     * 64 bits for seed
     * 64 bits for the mersenne prime (technically we only use 32 bits here as currently written)
     * 64 bits for every entry in the a and b array which are num_hashes long
     * 64 bits for every entry in the table array, and there are num_hashes*num_buckets such entries
     * 64 * ( 1 + 1 + 2*num_hashes + num_hashes*num_buckets)
     *
     */
    return 64 * (1 + 1 + 2*num_hashes + num_hashes*num_buckets) ;
} // End get_size_estimate_in_bits
