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
 ./volesti_lecount INSTANCE VOLUME_METHOD ROUNDING_METHOD 
```
_Note:_ (ROUNDING_METHOD is optional) 

**Example: for (volume method = sequence of balls, rounding method = SVD)**
```bash
 ./volesti_lecount instances/bipartite_0.5_008_0.txt sob SVD  
```

## Sample instances:
Currently the example contains two poset instances available as adjacency matrix of their respective DAG's (**D**irected **A**cyclic **G**raph).
The instances are as follows:
- [Bipartite Graph](https://en.wikipedia.org/wiki/Bipartite_graph) with 8 elements  
![bipartite_image](images/bipartite_0.5_008_0.png)  
**Exact number of linear extensions: 1504**   
<hr>  

- [Bayesian Network for Andes model](https://www.bnlearn.com/bnrepository/) with 8 elements  
![bayesian_andes_image](images/bayesiannetwork_andes_008_0.png)  
**Exact number of linear extensions: 28**  
<hr>  

Both instances are generated and visualized using scripts given in this [repo](https://github.com/ttalvitie/le-counting-practice).
