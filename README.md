## Practical volume computation, supplementary code for paper "Practical volume estimation by a new annealing schedule for cooling convex bodies" submitted to SoCG 2019

#### 0. Compile code
- Clone repository and switch to branch v_poly:  
```
git checkout v_poly
```
- Save `liblpsolve55.so` in folder `/usr/lib`. You will find it in `/external` folder.
- In folder `/test` compile C++ code by running:  
```
cmake .  
make  
```
#### 1. H-polytopes (Table 1):  
- To generate a unit cube in `dim` dimension and estimate the volume:  
```
./generate -cube -h -d dim
./vol -f1 cube_dim.ine -ban
```
For example:  
```
./generate -cube -h -d 20
./vol -f1 cube_20.ine -ban
```
- To generate a unit simplex in `dim` dimension and estimate the volume:  
```
./generate -simplex -h -d dim
./vol -f1 simplex_dim.ine -ban
```
- To generate a random H-polytope in `dim` with `k` facets dimension and estimate the volume:  
```
./generate -rh -d dim -m k
./vol -f1 random_h_poly_dim_k.ine -ban
```

#### 2. V-polytopes (Table 2):  
- To generate a cross polytope in `dim` dimension and estimate the volume:  
```
./generate -cross -v -d dim
./vol -f2 cross_dim.ext -ban
```
- To generate a unit simplex in `dim` dimension:  
```
./generate -simplex -v -d dim
./vol -f2 simplex_dim.ext -ban
```
- To generate a unit cube in `dim` dimension and estimate the volume:  
```
./generate -cube -v -d dim
./vol -f2 cube_dim.ext -ban
```
- To generate a random V-polytope in `dim` dimension with k vertices and estimate the volume:  
```
./generate -rv -d dim -m k
./vol -f2 random_v_poly_dim_k.ext -ban -r
```
Note: For random V-polytopes use the flag `-r` to round the polytope.

#### 3. Zonotopes (Table 3):  
- You can generate a random zonotope in dimension `dim` with `k` generators by running:  
```
./generate -zonotope -d dim -m k
```
- Estimate the volume using balls in MMC:  
```
./vol -f3 zonotope_dim_k.ext -ban
```
- Estimate the volume using h-polytopes in MMC:  
```
./vol -f3 zonotope_dim_k.ext -hpoly
```
- For example the following commands:  
```
./generate -zonotope -d 10 -m 15
./vol -f3 zonotope_10_15.ext -hpoly
```
Will generate a random 10-dimensional zonotope with 15 generators and estimate the volume by using h-polytopes in MMC.  
- To compute the exact volume run:  
```
./vol -f3 zonotope_10_15.ext -exact_zono
```

Note: If you wish to give a specific polytope as input use `.ine` file for a H-polytope and `.ext` file for a V-polytope or a zonotope. Keep the same format as in the generated files.

#### 4. Flags

- For H-polytopes the default random walk is Coordinate Directions HnR. Use flag `-rdhr` to use Random Directions HnR:  
```
./vol -f1 cube_dim.ine -ban -rdhr
```
- For V-polytopes and zonotopes the default random walk is Random Directions HnR. To use Coordinate Directions HnR use the flag `-cdhr`. For example:  
```
./vol -f2 cube_dim.ext -ban -cdhr
```
- Use flag `-WalkL` to set the step of the HnR (the default value is 1). For example:  
```
./vol -f1 cross_dim.ext -ban -WalkL 5
```
Will set the step equals to 5.
- Use flag `-e` to set the error (the default value is `0.1`). For example:  
```
./vol -f1 zonotope_dim_k.ext -ban -e 0.2
```
- Use flag `-WinLen` to set the length of the sliding window (the default value is 4d^2+250). For example:  
```
./vol -f1 cross_dim.ext -ban -WinLen 500
```
Will set the window's length `n=500`.
- Use flags `-l` and `-u` to set the test values for testR (r) and testL (r+\delta) respectively. For example:  
```
./vol -f1 cube_dim.ine -ban -l 0.01 -u 0.015
```
Will define ratios between `0.01` and `0.015` with high probability.
- Use flag `-nuN` to set the number of points that are generated in each step of the annealing schedule, from the convex body P_i of the previous step. For example:  
```
./vol -f3 zonotope_dim_k.ext -ban -nuN 1600 10
```
Wil sample 1600 points in total and split them to 10 sub-lists. So the degrees of freedom in each t-test will be 9 = 10-1.  

#### 5. Test PCA over-aproximations of a zonotope

- To compute the ratio for the PCA approximation of a zonotope that is described in a `.ext` file, use flag `-pca` and run:  
```
./vol -f3 zonotope_dim_k.ext -hpoly -pca
```

