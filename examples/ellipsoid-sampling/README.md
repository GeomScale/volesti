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
 ./ellipsoid2d-sampling > ellipsoid.txt
```
After this, you can run `python3 plot_pts.py` to plot the sampled points.

## Sample instance:
The ellipsoid is in the form: `(x-c)' A (x-c) <= 1`, where `A = LL'`, currently there is only sample instance (number of sample points = 1000)
```
L = [0.5, 0,
     1.5, 1.0];
```
This makes
```
A = [0.25, 0.75,
     0.75, 3.25];
```

**Sampled points:**
![sampled_points](sampled_points.png)
