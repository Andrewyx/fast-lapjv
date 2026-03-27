cc_library(
    name = "matrix",
    hdrs = ["matrix.h", "matrix.cpp"],
    includes = ["."],
    visibility = ["//visibility:public"],
)

config_setting(
    name = "opt_mode",
    values = {"compilation_mode": "opt"},
)

config_setting(
    name = "dbg_mode",
    values = {"compilation_mode": "dbg"},
)

cc_library(
    name = "lapjv",
    srcs = ["lap.cpp"],
    hdrs = ["lap.h"],
    copts = select({
        ":dbg_mode": [
            "-g",
            "-fopenmp-simd",
        ],
        "//conditions:default": [
            "-O3",
            "-march=native",
            "-ffast-math",
            "-fopenmp-simd",
            "-flto",
        ],
    }),
    linkopts = select({
        ":dbg_mode": [],
        "//conditions:default": ["-flto"],
    }),
    includes = ["."],
    visibility = ["//visibility:public"],
    deps = [
        ":matrix",
        "@highway//:hwy",
    ],
)
