
## Installation of the required libraries
C++ 17 is required

### Background information
Autopoint is implemented using [Autodiff Library](
https://github.com/autodiff/autodiff/), after experiments of comparing various automatic differentiating libraries.Read [here](https://gist.github.com/zhanggiene/8471601fa25ba9db90303661b0e2237b) to find out more!

### Steps to run the three examples

For plotting the results:
 * `examples/python_utilities/plot_samples.py` for plotting the outputs of the HMC sampler

The script has requirements listed in `requirements.txt`. Install via
```
pip install -r requirements.txt
```

You can also use the `--limits` argument to set {x, y, z}-limits for the joint plots.


### 1. User defined function where the pdf function is explicitly defined by the user. 

``` 
./userDefinedFunction_autopoint > userDefinedFunction_autopoint.txt

python ../python_utilities/plot_samples.py < userDefinedFunction_autopoint.txt
```
Here is the histogram of samples from user defined distribution.It is the same as the real pdf of the distribution when -1<x<1 since sampling is done is in truncated space.
![pdf](./1_2.png)

here is the density plot of the real distribution.
![1.2](./true_pdf_userdefined_func.png)

### 2. The pdf is defined using a mixture of gaussian distribution with data
1. generate data from the mixture of gaussian distribution with mean of dimension 1. 
```
python3 generate_multidimensional_gaussian_mixture_data.py -d=1 > data.txt
```
2. Data from the true distribution is used to define the pdf. 
```
./Gaussian_mixture_autopoint > Gaussian_mixture.txt
```
3. Plot the graph
```
python ../python_utilities/plot_samples.py < Gaussian_mixture.txt
```
![gaussian_mixture](./2.2.png)

### 3.The pdf is defined using a mixture of multi-dimensional gaussian distribution with multi-dimensional data. 

1. Generate data points from mixture of multi-dimensional gaussian distribution
d can be change to any integer number, but it should be the same as unsigned int dim variable in MultidimensionalGaussian_mixture_autopoint.cpp
```
python generate_multidimensional_gaussian_mixture_data.py -d=4 > data.txt 
```
2. Conduct hmc sampling

```
./MultidimensionalGaussian_mixture_autopoint
```
### Points to take note
it is advised to use NutsHamiltonianMonteCarloWalk instead of HamiltonianMonteCarloWalk.
