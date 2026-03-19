# LAPJV 


The objective of this project is to optimize the famous LAPJV algorithm in C++. 
The current C/C++ implementation is by https://github.com/yongyanghz/LAPJV-algorithm-c.

## Rules

- Implementations must abide by the matrix.h API as inputs. You may change the implementation of this data struct
- You may refer to the [Munkres Algorithm](https://github.com/saebyn/munkres-cpp) for reference
- Bazel must be used

## Tooling
- Benchmarking: `bazel build --config=max --config=dbg //test:lapbench`
- Benchmark Recent Comparison: `bazel run //test:compare_benchmarks -- $(pwd)/results.json`
- Benchmark Graphs: `bazel run //test:vis_benchmark -- $(pwd)/results.json`
- CPU Profiling: `sudo perf record -g ./bazel-bin/test/lapbench --benchmark_filter="BM_LAPJVSolve<float>/512"`
- Profiling Flamegraph`perf script | stackcollapse-perf.pl | flamegraph.pl > flamegraph.svg`
