# volestipy

## Installation
To compile the Python wrapper you have first need to get  [liblpsolve55.so](https://sourceforge.net/projects/lpsolve/) library. Let's assume you have it under `/usr/lib/lpsolve/`. Execute then the following commands.

```
cd volestipy
LDFLAGS="-L/usr/lib/lpsolve/" python setup.py build_ext -i
python setup.py install -f
```

## Run an example

```
python test1.py
```
