# VolEsti ( volume computation and sampling library)

"""Benchmark Suite for VolEsti and Dingo Comparison
This script benchmarks the performance of the VolEsti C++ library against the Dingo Python wrapper
for random walk sampling in high-dimensional polytopes.
It measures the time taken by both libraries to perform a specified number of walk steps
in polytopes of varying dimensions and computes the overhead percentage of Dingo over Volesti.
"""

# To make the path to the volesti binary a command-line argument so it doesn't just work on maintainers machine, instead of hardcoding it.

import argparse
import os
import subprocess
import sys

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Volesti-Dingo Benchmarking Suite")
    
    # Add the argument for the binary path
    parser.add_argument(
        '--bin', 
        type=str, 
        required=True,
        help="Path to the compiled volesti sampler binary (e.g., ./build/bin/sampling)"
    )
    
    # Add optional arguments for flexibility
    parser.add_argument('--dim', type=int, default=10, help="Dimension of the polytope")
    parser.add_argument('--steps', type=int, default=1000, help="Number of walk steps")

    args = parser.parse_args()

    # Verify the path exists before running
    if not os.path.exists(args.bin):
        print(f"Error: Binary not found at {args.bin}")
        sys.exit(1)

    # Use args.bin in your subprocess call
    # Example:
    # subprocess.run([args.bin, "-i", "cube.ext", "-n", str(args.steps)])
    
    print(f"Starting benchmark using binary: {args.bin}")
    # ... rest of your logic ...

if __name__ == "__main__":
    main()



# Third-Party Imports

import dingo
import time
import numpy as np
import pandas as pd
import subprocess

def run_audit(dimension, walk_steps, walk_type="ball_walk"):
    # 1. Setup a standard H-polytope (Unit Cube)
    # Ax <= b
    A = np.vstack([np.eye(dimension), -np.eye(dimension)])
    b = np.ones(2 * dimension)

    # 2. Benchmark Dingo (Python Wrapper)
    start_py = time.perf_counter()
    # Assuming dingo has a sampling function
    py_samples = dingo.sample(A, b, walk_type=walk_type, steps=walk_steps)
    end_py = time.perf_counter()
    
    # 3. Benchmark Volesti (Pure C++)
    # We save the polytope to a temporary file to read in C++
    np.savetxt("temp_poly.txt", np.hstack([A, b.reshape(-1, 1)]))
    
    start_cpp = time.perf_counter()
    # Call the volesti binary via subprocess
    subprocess.run(["./volesti_sampler", "-i", "temp_poly.txt", "-n", str(walk_steps)], 
                   stdout=subprocess.DEVNULL)
    end_cpp = time.perf_counter()

    return {
        "Dimension": dimension,
        "Dingo_Time": end_py - start_py,
        "Volesti_Time": end_cpp - start_cpp,
        "Overhead_Percentage": ((end_py - start_py) / (end_cpp - start_cpp) - 1) * 100
    }

# Execute for multiple dimensions
results = [run_audit(d, 1000) for d in [10, 50, 100, 500]]
print(pd.DataFrame(results))