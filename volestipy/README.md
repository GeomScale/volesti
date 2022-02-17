# volestipy

## Installation

### Dependencies

To compile the Python wrapper you first need to get the [liblpsolve55.so](https://sourceforge.net/projects/lpsolve/) library. 
You may need to have a look [here](http://lpsolve.sourceforge.net/5.5/) to get this right. 
Let's assume you have it under `/usr/lib/lpsolve/`.

[//]: # (**Reminder**)
[//]: # (As I recall `lpsolve` is not that straightforward to get. )
[//]: # (I think it would be useful to describe a *how to get it* thoroughly. )


You also need `cython`, `numpy` and `setuptools`. We need to install pip to be able to download these modules. If you do not have pip already you can check it out [here](https://pip.pypa.io/en/stable/installation/)

To download the packages
```
pip install Cython
```
```
pip install numpy
```
```
pip install setuptools
```

Finally, the [Gurobi solver](https://www.gurobi.com/) needs to be installed as well.
Through the [Download center](https://www.gurobi.com/downloads/) you may download its Python interface along with a license; the latter is a text file called `gurobi.lic`.
Without a license *gurobipy* will not be able to perform neither will *volestipy*. 

You can follow the inscriptions described [here](https://support.gurobi.com/hc/en-us/articles/360044290292-Installing-Gurobi-for-Python) to get ```gurobipy```. 

Once you complete these steps, make sure that `gurobipy` is now available for your Python environment. 

```
haris@XPS-13-9343:~$ python3
Python 3.6.9 (default, Apr 18 2020, 01:56:04) 
[GCC 8.4.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import gurobipy
>>> gurobipy.
gurobipy.AttrConstClass(     gurobipy.GenExpr(            gurobipy.QConstr(            gurobipy.atexit              gurobipy.gurobi(             gurobipy.multidict(          gurobipy.resetParams(
gurobipy.Batch(              gurobipy.GurobiError(        gurobipy.QuadExpr(           gurobipy.bi                  gurobipy.help(               gurobipy.numbers             gurobipy.setParam(
gurobipy.CallbackClass(      gurobipy.Iterable(           gurobipy.SOS(                gurobipy.dis                 gurobipy.inspect             gurobipy.operator            gurobipy.sys
gurobipy.Column(             gurobipy.LinExpr(            gurobipy.StatusConstClass(   gurobipy.disposeDefaultEnv(  gurobipy.itertools           gurobipy.or_(                gurobipy.system(
gurobipy.Constr(             gurobipy.MLinExpr(           gurobipy.TempConstr(         gurobipy.exprfactory(        gurobipy.izip(               gurobipy.os                  gurobipy.tupledict(
gurobipy.Env(                gurobipy.MQuadExpr(          gurobipy.Var(                gurobipy.exprfactory_iter(   gurobipy.logging             gurobipy.paramHelp(          gurobipy.tuplelist(
gurobipy.ErrorConstClass(    gurobipy.MVar(               gurobipy.abs_(               gurobipy.fnmatch             gurobipy.math                gurobipy.quicksum(           gurobipy.types
gurobipy.GRB(                gurobipy.Model(              gurobipy.all_(               gurobipy.gc                  gurobipy.max_(               gurobipy.re                  gurobipy.writeParams(
gurobipy.GRBStringIO(        gurobipy.ParamClass(         gurobipy.and_(               gurobipy.getParamInfo(       gurobipy.min_(               gurobipy.read(               
gurobipy.GenConstr(          gurobipy.ParamConstClass(    gurobipy.any_(               gurobipy.glob                gurobipy.models(             gurobipy.readParams(         
>>> gurobipy.
```


### Install *volestipy*

After getting the dependencies, download the *volesti* repository from GitHub and run the following commands.

```
cd volestipy
LDFLAGS="-L/usr/lib/lpsolve/" python3 setup.py install --user
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
