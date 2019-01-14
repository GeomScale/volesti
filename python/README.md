## Compilation

Currently, we assume that [https://github.com/facebookresearch/faiss](faiss) is installed in the system. Alternatively you can modify the `compile` script to point it to its library file.

## Running

Compilation generates a `volesti.cpython-*-linux-gnu.so` file which you should add to your `PYTHONPATH` or have it in your working directory. Then you can use the exposed `volesti` functions as shown in the  `test.py` file.
