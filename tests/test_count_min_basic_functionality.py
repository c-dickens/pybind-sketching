from count_min_sketch import CountMinSketch
import pytest
import numpy as np
from fixtures import fixed_sketch_setup, make_data, make_count_min_with_seed

def test_string_methods(fixed_sketch_setup):
    sketch_params = fixed_sketch_setup
    c = CountMinSketch(*sketch_params)
    summary = c.to_string()
    print(summary)
    print(c)

def test_constructors(fixed_sketch_setup):
    """
    Tests the construction methods and basic funtionality such as getting bucket numbers or seed
    :return:
    :rtype:
    """
    sketch_params = [num_hashes, num_buckets, seed] = fixed_sketch_setup
    c = CountMinSketch(*sketch_params)  # Unpack the params list into the sketch constructor

    assert c.get_num_hashes() == num_hashes
    assert c.get_num_buckets() == num_buckets
    assert c.get_seed() == seed
    assert c.get_table_shape() == (num_hashes, num_buckets)
    assert c.get_config() == [num_hashes, num_buckets, seed]
    assert c.get_table() == [[0]*num_buckets]*num_hashes

    # This line tells pytest that the next block should raise a TypeError exception.
    with pytest.raises(ValueError):
        c.update(-1)
        raise ValueError


def test_single_updates(fixed_sketch_setup):
    """
    Tests updates are performed as expected under the hash functions.
    :return:
    :rtype:

    Nb. for this seed we know from inspection in the main.cpp that 1 is hashed to the buckets [9, 1, 3].
    No extra weight is added so each of these sketch entries should be equal to 1
    """
    c = CountMinSketch(*fixed_sketch_setup)
    assert c.get_estimate(1) == 0
    c.update(1)
    sk = c.get_table()
    assert sk[0][9] == 1
    assert sk[1][1] == 1
    assert sk[2][3] == 1
    c.update(1,9)
    sk = c.get_table()
    assert sk[0][9] == 10
    assert sk[1][1] == 10
    assert sk[2][3] == 10
    # no other items inserted so estimation should be exact.
    assert c.get_estimate(1) == 10
    c.update(1, -10) # Remove all weight so the buckets should now have a zero.
    sk = c.get_table()
    assert sk[0][9] == 0
    assert sk[1][1] == 0
    assert sk[2][3] == 0

def test_parameter_suggestions():
    """
    Tests the suggesttion functions return the correct sketch size parameters.
    :return:
    :rtype:
    """

    # Test number of bucket suggestions: w = ceil(e / epsilon)
    errors = [0.2, 0.1, 0.05, 0.01, 1E-5]
    buckets = [14, 28, 55, 272, 271829]
    for e, b in zip(errors, buckets):
        assert CountMinSketch.suggest_num_buckets(e) == b
    with pytest.raises(ValueError) as excinfo:
        CountMinSketch.suggest_num_buckets(-1.0)
        assert "Relative error must be at least 0." in str(excinfo.value)
        raise ValueError("Relative error must be at least 0.")


    # NTest number of hashes suggestions: d = ceil ln(1 / delta) )
    assert CountMinSketch.suggest_num_hashes(0.682689492) == 2
    assert CountMinSketch.suggest_num_hashes(0.954499736) == 4
    assert CountMinSketch.suggest_num_hashes(0.997300204) == 6
    assert CountMinSketch.suggest_num_hashes(1 - 1E-5) == 12
    for _ in [10.0, -1.0]:
        with pytest.raises(ValueError) as excinfo:
            CountMinSketch.suggest_num_hashes(_)
            assert "Confidence must be between 0 and 1.0 (inclusive)." in str(excinfo.value)
            raise ValueError("Confidence must be between 0 and 1.0 (inclusive).")


@pytest.mark.parametrize("make_count_min_with_seed", [[0.1, 0.99, 1]], indirect=True)
def test_single_stream(make_data, make_count_min_with_seed):
    """
    Tests that the frequency estimates are within the bounds on a single stream.
    :return:
    :rtype:
    """
    data, frequencies = make_data
    c = make_count_min_with_seed["sketch"]

    # populate the sketch
    for item, weight in zip(data, frequencies):
        c.update(item, weight)

    # Output the results:
    print("{:<4} {:>6} {:>9} {:>6} {:>6}".format('Item', 'Weight', 'Estimate', 'Upper', 'Lower'))
    for item, weight in zip(data, frequencies):
        est = c.get_estimate(item)
        upp = c.get_upper_bound(item)
        low = c.get_lower_bound(item)
        print("{:<4} {:>6} {:>9} {:>6} {:>6}".format(item, weight, est, upp, low))
        assert est <= upp
        assert est >= low


@pytest.mark.parametrize("make_count_min_with_seed", [[0.1, 0.99, 1]], indirect=True)
def test_merge(make_count_min_with_seed):
    """
    Tests the merging capability of the sketh.
    :return:
    :rtype:
    """
    s = make_count_min_with_seed["sketch"]
    n_hashes = make_count_min_with_seed["n_hashes"]
    n_buckets = make_count_min_with_seed["n_buckets"]
    seed = make_count_min_with_seed["seed"]

    s_config = s.get_config()
    assert s_config[0] == n_hashes
    assert s_config[1] == n_buckets
    assert s_config[2] == seed

    with pytest.raises(ValueError) as excinfo:
        s.merge(s)
        assert "Cannot merge a sketch with itself." in str(excinfo.value)
        raise ValueError("Cannot merge a sketch with itself.")

    # Generate sketches that we cannot merge into ie they disagree on at least one of the config entries
    # TODO: implement more general merge procedure so different hashes or buckets are permitted.
    s1 = CountMinSketch(n_hashes+1, n_buckets, seed)
    s2 = CountMinSketch(n_hashes, n_buckets+1, seed)
    s3 = CountMinSketch(n_hashes, n_buckets, seed+1)
    for sk in [s1, s2, s3]:
        with pytest.raises(ValueError) as excinfo:
            s.merge(sk)
            assert "Incompatible sketch config." in str(excinfo.value)
            raise ValueError("Incompatible sketch config.")

    # Passing cases
    t = CountMinSketch(*s_config)
    data = list(range(5))
    for d in data:
        s.update(d)
        t.update(d)
    estimates = {
        "s": [s.get_estimate(d) for d in data],
        "t": [s.get_estimate(d) for d in data]
    }
    original_sketch = np.array(s.get_table())
    s.merge(t)
    assert s.get_total_weight() == 2*t.get_total_weight()
    for i, d in enumerate(data):
        assert s.get_estimate(d) == (estimates["s"][i] + estimates["t"][i])
    # Check the merged array is the old one plus the t array.
    assert np.all(np.array(s.get_table()) == (original_sketch + np.array(t.get_table())))

@pytest.mark.parametrize("make_count_min_with_seed", [[0.1, 0.99, 1]], indirect=True)
def test_merge_overload(make_count_min_with_seed):
    """
    Tests the merging capability of the sketh.
    :return:
    :rtype:
    """
    s = make_count_min_with_seed["sketch"]
    n_hashes = make_count_min_with_seed["n_hashes"]
    n_buckets = make_count_min_with_seed["n_buckets"]
    seed = make_count_min_with_seed["seed"]
    s_config = s.get_config()
    assert s_config[0] == n_hashes
    assert s_config[1] == n_buckets
    assert s_config[2] == seed

    with pytest.raises(ValueError) as excinfo:
        s + s
        assert "Cannot merge a sketch with itself." in str(excinfo.value)
        raise ValueError("Cannot merge a sketch with itself.")

    # Generate sketches that we cannot merge into ie they disagree on at least one of the config entries
    # TODO: implement more general merge procedure so different hashes or buckets are permitted.
    s1 = CountMinSketch(n_hashes+1, n_buckets, seed)
    s2 = CountMinSketch(n_hashes, n_buckets+1, seed)
    s3 = CountMinSketch(n_hashes, n_buckets, seed+1)
    for sk in [s1, s2, s3]:
        with pytest.raises(ValueError) as excinfo:
            s + sk
            assert "Incompatible sketch config." in str(excinfo.value)
            raise ValueError("Incompatible sketch config.")

    # Passing cases
    t = CountMinSketch(*s_config)
    data = list(range(5))
    for d in data:
        s.update(d)
        t.update(d)
    estimates = {
        "s" : [s.get_estimate(d) for d in data],
        "t" : [s.get_estimate(d) for d in data]
    }
    original_sketch = np.array(s.get_table())
    s + t
    assert s.get_total_weight() == 2*t.get_total_weight()
    for i, d in enumerate(data):
        assert s.get_estimate(d) == (estimates["s"][i] + estimates["t"][i])
    # Check the merged array is the old one plus the t array.
    assert np.all(np.array(s.get_table()) == (original_sketch + np.array(t.get_table())))
