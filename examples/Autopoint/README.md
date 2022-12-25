
## Installation of the required libraries
C++ 17 is required

The offical installation guide for Autodiff can be found [here](
https://autodiff.github.io/installation/ ). Here are some useful links for installing some prerequisites for Autodiff library([Catch2 issue1](https://github.com/etaler/Etaler/issues/33),[Catch2 issue2](https://stackoverflow.com/questions/65098604/catch2-installation-on-ubuntu-20-04-include-catch2-catch-hpp
),[libeigen3-dev](https://zoomadmin.com/HowToInstall/UbuntuPackage/libeigen3-dev),[pybind11](https://stackoverflow.com/questions/46961942/pybind11-linux-building-tests-failure-could-not-find-package-configuration-fi ))


### Background information
Autopoint is implemented using AUtodiff Library , after experiments of comparing various automatic differentiating libraries.Read [here](https://gist.github.com/zhanggiene/8471601fa25ba9db90303661b0e2237b) to find out more!

### Steps to run the three examples




### 1. User defined function where the pdf function is explicitly defined by the user. 

``` 
./userDefinedFunction_autopoint > userDefinedFunction_autopoint.txt

python plot_hmc.py < userDefinedFunction_autopoint.txt
```
Here is the histogram of samples from user defined distribution.It is the same as the real pdf of the distribution when -1<x<1 since sampling is done is in truncated space.
![pdf](./1_2.png)

here is the density plot of the real distribution.
![1.2](./true_pdf_userdefined_func.png)

### 2. The pdf is defined using a mixture of gaussian distribution with data
1. generate data from the mixture of gaussian distribution
```
python generate_gaussian_mixture_data.py > data.txt
```
2. Data from the true distribution is used to define the pdf. 
```
./Gaussian_mixture_autopoint > Gaussian_mixture.txt
```
3. Plot the graph
```
python plot_hmc.py < Gaussian_mixture.txt
```
![gaussian_mixture](./2.2.png)

### 3.The pdf is defined using a mixture of multi-dimensional gaussian distribution with multi-dimensional data. 

1. Generate data points from mixture of multi-dimensional gaussian distribution
```
python generate_multidimensional_gaussian_mixture_data.py > data.txt 
```
2. Conduct hmc sampling

```
./MultidimensionalGaussian_mixture_autopoint
```
### Points to take note
it is advised to use NutsHamiltonianMonteCarloWalk instead of HamiltonianMonteCarloWalk.
