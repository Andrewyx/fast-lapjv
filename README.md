# LAPJV 


The objective of this project is to optimize the famous LAPJV algorithm in C++. 
The current C/C++ implementation is based off of https://github.com/yongyanghz/LAPJV-algorithm-c, and inspired by:
- https://github.com/google/highway
- https://github.com/src-d/lapjv
- https://agner.org/optimize/optimizing_cpp.pdf
- https://github.com/saebyn/munkres-cpp
Gemini was used in the creation of this project.

## Requirements
The project uses Bazel as its main build system. Add the following to your `MODULE.bazel` file to add this to your project.
```
http_archive(

)
```

## Usage
Compiling your project with the `-c opt` Bazel flag in your project will invoke the optimized versions of this library.


## Benchmarking
- Benchmark Recent Comparison: `bazel run //test:compare_benchmarks -- $(pwd)/results.json`
- Benchmark Graphs: `bazel run //test:vis_benchmark -- $(pwd)/results.json`

## Tracing

![Flamegraph](flamegraph/flamegraph.svg)

Trace details and symbols are only generated in debug mode, hence to generate a flamegraph use the following steps:
1) `bazel build --config=dbg //test:lapbench`
2) `perf record -g -o flamegraph/perf.data ./bazel-bin/test/lapbench --benchmark_filter="BM_LAPJVSolve<double>/512"`
3) `perf script -i flamegraph/perf.data | stackcollapse-perf.pl | flamegraph.pl > flamegraph/flamegraph.svg`
Of course, ensure that you have installed the [flamegraphs project](https://github.com/brendangregg/flamegraphv) from the 
famous Brendan Gregg. Furthermore, make sure to add `stackcollapse-perf.pl` and `flamegraph.pl` to your `PATH`. 
