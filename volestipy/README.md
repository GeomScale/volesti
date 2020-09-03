# volestipy

## Installation

### Dependencies

To compile the Python wrapper you first need to get the [liblpsolve55.so](https://sourceforge.net/projects/lpsolve/) library. Let's assume you have it under `/usr/lib/lpsolve/`.
You also need `cython`, `numpy` and `setuptools`. In `debian` systems you can get then by
```
sudo apt-get install cython python-numpy python-setuptools
```

You also need to get the [Gurobi solver](https://www.gurobi.com/).
Through the [Download center](https://www.gurobi.com/downloads/) you may download its Python interface along with a license.
Without a license *gurobipy* will not be able to perform neither will *volestipy*. 


### Install *volestipy*

After getting the dependencies, download the *volesti* repository from GitHub and run the following commands.

```
cd volestipy
LDFLAGS="-L/usr/lib/lp_solve/" python3 setup.py install --user
```

## Run an example
Before you use the python library make sure the `liblpsolve55.so` library is present in `LD_PATH_LIBRARY`, for example by adding the following line to your `.bashrc` file:

```
export LD_LIBRARY_PATH=/usr/lib/lpsolve:$LD_LIBRARY_PATH
```
You should now be able to run the two test files, e.g.:
```
python3 tests/test_read_bigg_files.py
```

A test is available for each ```volestipy``` function.