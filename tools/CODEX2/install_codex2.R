#!/usr/bin/env Rscript

if (!require(CODEX2)) {
    user_lib = Sys.getenv("R_LIBS_USER")
    dir.create(path = user_lib, recursive = TRUE)
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", lib = user_lib)
    }
    if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools", lib = user_lib)
    }
    BiocManager::install("CODEX", version = "3.8", lib = user_lib)
    devtools::install_github("yuchaojiang/CODEX2/package", lib = user_lib, lib.loc = user_lib)
}

