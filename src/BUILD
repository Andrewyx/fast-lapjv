cc_library(
    name = "matrix",
    hdrs = ["matrix.h", "matrix.cpp"],
    includes = ["."],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "lapjv",
    srcs = ["lap.cpp"],
    hdrs = ["lap.h"],
    includes = ["."],
    deps = [
        ":matrix",
        "@highway//:hwy",
    ],
    visibility = ["//visibility:public"],
)