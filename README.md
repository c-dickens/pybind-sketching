pybind_sketching
==============

# CHARLIE'S EDITS
The environment to activate is `conda activate pybind-example`
Note that the install command should be `pip install .`
1. Changed module name to pybind_sketching
2. After updating the c++ code, one should run the `./install_and_test.sh` executable script


## To Do:
*Why would we want to use a CountMin rather than a frequent items sketch?*
0. Check that the estimation functions should return uint or int
0a. Add a streaming heap to keep a list of the current heavy items?   
1. Separate the source code into header and source files
3. Add more unit tests for other queries (point, range, heavy hitter, quantile).  See 
   [here](https://barnasahadotcom.files.wordpress.com/2016/01/lec3-haritha-1.pdf) and 
   [here](https://crahen.github.io/algorithm/stream/count-min-sketch-point-query.html) for handy overviews.
4. Add string functionality
5. Add serialization?
7. Add a `get_numpy_table()`.
8. Use a better estimator as show [here](https://dl.acm.org/doi/10.1145/3219819.3219975)
9. Measure update time and merge time.
10. Optimise the performance using bitwise operations for the modular arithmetic.
11. Scatter plot of error vs item frequency (which should be somewhat inverse of count vs frequency)
12. Add catch into the src directory for direct c++ testing?





[![Gitter][gitter-badge]][gitter-link]

|      CI              | status |
|----------------------|--------|
| Linux/macOS Travis   | [![Travis-CI][travis-badge]][travis-link] |
| MSVC 2015            | [![AppVeyor][appveyor-badge]][appveyor-link] |
| conda.recipe         | [![Conda Actions Status][actions-conda-badge]][actions-conda-link] |
| pip builds           | [![Pip Actions Status][actions-pip-badge]][actions-pip-link] |
| [`cibuildwheel`][]   | [![Wheels Actions Status][actions-wheels-badge]][actions-wheels-link] |

[gitter-badge]:            https://badges.gitter.im/pybind/Lobby.svg
[gitter-link]:             https://gitter.im/pybind/Lobby
[actions-badge]:           https://github.com/pybind/python_example/workflows/Tests/badge.svg
[actions-conda-link]:      https://github.com/pybind/python_example/actions?query=workflow%3A%22Conda
[actions-conda-badge]:     https://github.com/pybind/python_example/workflows/Conda/badge.svg
[actions-pip-link]:        https://github.com/pybind/python_example/actions?query=workflow%3A%22Pip
[actions-pip-badge]:       https://github.com/pybind/python_example/workflows/Pip/badge.svg
[actions-wheels-link]:     https://github.com/pybind/python_example/actions?query=workflow%3AWheels
[actions-wheels-badge]:    https://github.com/pybind/python_example/workflows/Wheels/badge.svg
[travis-link]:             https://travis-ci.org/pybind/python_example
[travis-badge]:            https://travis-ci.org/pybind/python_example.svg?branch=master&status=passed
[appveyor-link]:           https://ci.appveyor.com/project/wjakob/python-example
<!-- TODO: get a real badge link for appveyor -->
[appveyor-badge]:          https://travis-ci.org/pybind/python_example.svg?branch=master&status=passed

An example project built with [pybind11](https://github.com/pybind/pybind11).
This requires Python 3.6+; for older versions of Python, check the commit
history.

Installation
------------

 - clone this repository
 - `pip install ./python_example`

CI Examples
-----------

There are examples for CI in `.github/workflows`. A simple way to produces
binary "wheels" for all platforms is illustrated in the "wheels.yml" file,
using [`cibuildwheel`][]. You can also see a basic recipe for building and
testing in `pip.yml`, and `conda.yml` has an example of a conda recipe build.


Building the documentation
--------------------------

Documentation for the example project is generated using Sphinx. Sphinx has the
ability to automatically inspect the signatures and documentation strings in
the extension module to generate beautiful documentation in a variety formats.
The following command generates HTML-based reference documentation; for other
formats please refer to the Sphinx manual:

 - `cd python_example/docs`
 - `make html`

License
-------

pybind11 is provided under a BSD-style license that can be found in the LICENSE
file. By using, distributing, or contributing to this project, you agree to the
terms and conditions of this license.

Test call
---------

```python
import python_example
python_example.add(1, 2)
```

[`cibuildwheel`]:          https://cibuildwheel.readthedocs.io
