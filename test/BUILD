load("@pip//:requirements.bzl", "requirement")

cc_test(
    name = "laptest",
    srcs = ["laptest.cpp"],
    deps = [
        "//src:lapjv",
        "//src:matrix",
        "@googletest//:gtest_main",
    ],
)

cc_binary(
    name = "lapbench",
    srcs = ["lapbench.cpp"],
    linkopts = ["-static"],
    linkstatic = True,
    deps = [
        "//src:lapjv",
        "@google_benchmark//:benchmark",
    ],
)

py_binary(
    name = "compare_benchmarks",
    srcs = ["compare_benchmarks.py"],
)

py_binary(
    name = "vis_benchmark",
    srcs = ["vis_benchmark.py"],
    deps = [
        requirement("matplotlib"),
    ],
)

cc_binary(
    name = "test_hwy",
    srcs = ["test_hwy.cpp"],
    deps = ["@highway//:hwy"],
)
