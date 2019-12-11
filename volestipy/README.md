# volestipy

## Installation
To compile the Python wrapper you have first need to get  [liblpsolve55.so](https://sourceforge.net/projects/lpsolve/) library. Let's assume you have it under `/usr/lib/lpsolve/`. 
You also need `cython`, `numpy` and `setuptools`. In `debian` systems you can get then by 
```
sudo apt-get install cython python-numpy python-setuptools
```

Then execute then the following commands.

```
cd volestipy
LDFLAGS="-L/usr/lib/lpsolve/" python setup.py build_ext -i
python setup.py install -f
```

## Run an example
Before you use the python library make sure the `liblpsolve55.so` library is present in `LD_PATH_LIBRARY`, for example by adding the following line to your `.bashrc` file:
```
export LD_LIBRARY_PATH=/usr/lib/lpsolve:$LD_LIBRARY_PATH

```
You shoudl now be able to run the two test files, e.g.:

```
python test1.py
```
