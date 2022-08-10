# Running the Examples

You have to install the latest mkl
[MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html)
For example, it is installed in the directory "/opt/intel/oneapi/mkl/2022.1.0". 

Change the environment variable:
```bash
export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64
```
Change the variable in CMakeLists.txt
```
set(MKLROOT /opt/intel/oneapi/mkl/2022.1.0)
```

Build the example by running the following commands in this directory.

```bash
cmake . -DLP_SOLVE=_PATH_TO_LIB_FILE
make
```  
You have to specify the path to liblpsolve55.so/dll/dylib.  
For example: -DLP_SOLVE=/usr/lib/lpsolve/liblpsolve55.so

Now, you can the HMC and ODE examples and generate their outputs
```
./simple_hmc >hmc_samples.txt
./simple_ode >ode_points.txt
```

For plotting the results there are two scripts
 * `plot_hmc.py` for plotting the outputs of the HMC sampler
 * `plot_ode.py` for plotting the outputs of the ODE solver

These scripts have requirements listed in `requirements.txt`. Install via
```
pip install -r requirements.txt
```

Plot the results as
```
  python plot_hmc.py <hmc_samples.txt
  python plot_ode.py <ode_points.txt
```

You can also use the `--limits` argument to set {x, y, z}-limits for the joint plots.

## Example Results

We sample with a negative log-density equal to  `f(x) = 2 * x^T x + sum(x)`, using 3 integrator steps
with step size 0.23 and draw 50000 samples, buring 25000.

**Marginals**

![Distributions](hmc_plot_marginals.png)

**Samples**

![Samples](hmc_plot_samples.png)

**Scatter Plot**

![Scatter](hmc_plot_scatter.png)
