## Compilation
Build the example by running the following commands in this directory.

```bash
cmake . -DLP_SOLVE=_PATH_TO_LIB_FILE
make
```  
You have to specify the path to liblpsolve55.so/dll/dylib.  
For example: -DLP_SOLVE=/usr/lib/lpsolve/liblpsolve55.so

## Usage:
```bash
 ./order_polytope poset_data.txt
```  
where `poset_data.txt` is the file containing the poset data as follows:
- first line of file tells the number of elements in the poset data
- next `m` lines represent the order relations of the poset as a pair of indices, for example:  
a line containing `i j` means A<sub>i</sub> <= A<sub>j</sub>. ( 0 <= `i`, `j` < n)
    