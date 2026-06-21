# Volesti Benchmark Suite

This project is a benchmarking suite for testing various random walk algorithms on high-dimensional convex polytopes found inside the volesti library

## Usage

Create a build folder, and run Cmake and Make:

```bash
mkdir build && cd build
cmake ..
make
```

## Command-Line Arguments

While not necessary, the executable accepts the following optional parameters

### `-h` or `--help`

Prints the help message and available options, then exits safely.

---

### `-c <path>` or `--config <path>`

Specify a custom path to your JSON configuration file.

- **Default:** `../config/walk_config.json`

---

### `-d <number>` or `--dim <number>`

Set the dimension for generating standard polytopes (e.g., Cube, Simplex).

- **Default:** The `"dimension"` value in your JSON config.
- **Note:** *Custom Polytopes* will use its own dimesnion based on the csv files.

---

### `-p <name>` or `--polytope <name>`

Choose the mathematical shape of the polytope.

- **Valid Options:** Cube, Simplex, Birkhoff, Cross, OrderPolytope, Custom

---

### `-w <name>` or `--walk <name>`

Specify exactly which random walk to run. Walk names are case-sensitive.

- **Valid Options:**  
  All, BallWalk, BilliardWalk, AcceleratedBilliardWalk, SparseBilliardWalk, CDHRWalk, RDHRWalk, DikinWalk, JohnWalk, VaidyaWalk, GaussianBallWalk, GaussianCDHRWalk, BilliardShakeAndBakeWalk, ShakeAndBakeWalk, BCDHRWalk, BRDHRWalk
                                                                  |

## Examples

### Run the default benchmark

(Uses the configuration file to determine the experiment parameters)

```bash
./benchmark_run
```

### Test a specific algorithm at a specific dimension and a specific Polytope

(Runs only the Accelerated Billiard Walk in 25 dimensions for the Simplex polytope)

```bash
./benchmark_run -p Simplex -w AcceleratedBilliardWalk -d 25
```

### Use your custom csv polytope

```bash
./benchmark_run -p Custom -w BilliardWalk
```

### 3. Use shorthand flags for a custom configuration

(Runs the Billiard Walk in 100 dimensions using a custom JSON file)

```bash
./benchmark_run -w BilliardWalk -d 100 -c ../config/walk_config.json
```

You actually don't have to use any argument at all. Just run the "benchmark_run" file and it will retrieve all the necessary parameters from the configuration file.

## Configuration (JSON)

Algorithm-specific parameters are managed entirely via the JSON configuration file. To change them, simply edit `config/walk_config.json`.
