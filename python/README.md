## Compilation

Currently, we assume that [https://github.com/facebookresearch/faiss](faiss) is installed in the system. Alternatively you can modify the `compile` script to point it to its library file.
To generate the python port, do the following:

```bash
$ virtualenv env (ideally have python3.6+)
$ source env/bin/activate
$ pip install pybind11 numpy
$ ./compile
```

## Running

Compilation generates a `volesti*.so` file which you should add to your `PYTHONPATH` or have it in your working directory. Then you can use the exposed `volesti` functions as shown in the  `test.py` file.
