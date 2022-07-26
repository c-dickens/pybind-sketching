import pytest
from count_min_sketch import CountMinSketch

@pytest.fixture
def fixed_sketch_setup():
    """
    A simple fixed sketch setup.
    :return:
    :rtype:
    """
    return [3, 10, 1]

@pytest.fixture()
def make_data():
    """
    Generates a simple synthetic dataset.
    """
    number_of_items = 10
    data = [i for i in range(number_of_items)]
    frequencies = [1 << (number_of_items - i) for i in range(number_of_items)]
    return data, frequencies

@pytest.fixture()
def make_count_min_with_seed(request):
    """
    The request is a list of [relative_error:float, confidence:float, seed:int]
    :param request:
    :type request:
    :return:
    :rtype:
    """
    [relative_error, confidence, seed] = request.param
    n_buckets = CountMinSketch.suggest_num_buckets(relative_error)
    n_hashes = CountMinSketch.suggest_num_hashes(confidence)

    # Return input parameters as this seems a simple way to get them into the scope of the function calling the fixture.
    all_parameters = {
        "sketch" : CountMinSketch(n_hashes, n_buckets, seed),
        "n_hashes" : n_hashes,
        "n_buckets" : n_buckets,
        "seed" : seed,
        "relative_error" : relative_error,
        "confidence" : confidence
    }
    return all_parameters