## Compilation

Create a build directory.
Build the example by running the following commands in the build directory you created.

```bash
cmake
make
```
You might have to specify the path to liblpsolve55.so/dll/dylib. Try these:

```bash
cmake . -DLP_SOLVE=_PATH_TO_LIB_FILE
make
``` 
TO find where this path is Try this command

   find /usr -name "liblpsolve55.so" 2>/dev/null  
   For example: -DLP_SOLVE=/usr/lib/lpsolve/liblpsolve55.so

If you dont find it try downloading it. Most common way is:

  sudo apt update
  sudo apt install lpsolve55 lpsolve55-dev

## Usage:
```bash
 ./general_sampling_time_limit
```

***How to use general_sampling_time_limit.cpp***

# Inside main:
  ->  comment or uncomment what method you want to run.
  ->  There is a dim vector. There define the dimensions you want to run the sampler for.
  ->  In addition, select what polytope you want to sample from. Comment or uncomment one already there or 
     create your own.
  ->  Choose an angle to rotate your polytope. Do this to avoid some methods taking advantage of a well 
     positioned polytope.

# Before main:
  ->  In function compute_batch_size, choose the batch size for each method. Some methods takes more points to
     deliver required ESS so increase or decrease based on how long you want to wait.
  ->  In function compute_ess and line 188 comment or uncomment to print ESS for every batch.
     This will help a lot with debugging and choosing the correct batch size. If you see ESS increasing by 10
     each batch, then you need to increase batch size or change walk_len.
  ->  In function set_walk_len, choose the walk_len parameter for each method. Based on experimental and theoretical 
     results, the values defined look ok. Walk len has to do with how many point the generator skips in order to 
     increase mix of the samples and their independence.

# Notes

A dynamic batch size was added that decides the batch size alone, without user input. It starts with 10 times the target ESS samples and adjusts 
the next batch size according to the ess per sample ratio and the remaining ess left. Right now the code uses this dynamic batch size.
I wll add a fail safe if ess is the same after a batch size.

To ask for a specific number of samples, set target ESS = 1 and choose batch size as big as the samples you want to generate.
 
