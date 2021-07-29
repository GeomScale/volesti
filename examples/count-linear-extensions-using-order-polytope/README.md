## Compilation
In folder examples, first run cmake, to create the makefile:

```bash
cmake .
```

Then, in folder examples/count-linear-extensions-using-order-polytope compile and build using the makefile:

```bash
make
```

## Usage:
```bash
 ./volesti_lecount_order_polytope INSTANCE
```

**Example:
```bash
 ./volesti_lecount_order_polytope instances/bipartite_0.5_008_0.txt
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
