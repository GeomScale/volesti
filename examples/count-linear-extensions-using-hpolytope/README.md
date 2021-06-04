## Compilation
In folder examples, first run cmake, to create the makefile:

```bash
cmake .
```

Then, in folder examples/count-linear-extensions-using-hpolytope compile and build using the makefile:

```bash
make
```

## Usage:
```bash
 ./volesti_lecount INSTANCE VOLUME_METHOD ROUNDING_METHOD 
```
_Note:_ (ROUNDING_METHOD is optional) 

**Example: for (volume method = sequence of balls, rounding method = SVD)**
```bash
 ./volesti_lecount instances/bipartite_0.5_008_0.txt sob SVD  
```
