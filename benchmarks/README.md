# Volesti-Dingo Benchmarking Suite

This suite provides tools to measure the performance gap and data-transfer overhead between the C++ core (`volesti`) and the Python bindings (`dingo`).

## How to Run
1. Ensure `volesti` is built and the binaries are in your path.
2. Install dependencies: `pip install dingo-python pandas matplotlib`
3. Execute: `python3 benchmark_suite.py`

## Current Metrics
- **Throughput:** Samples per second.
- **Latency:** Time to first sample.
- **Overhead:** Percentage of time spent in the Python wrapper vs. C++ core.