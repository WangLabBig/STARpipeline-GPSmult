options(repos = c(
    TUNA = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
    ALIY = "https://mirrors.aliyun.com/CRAN/"
))

# Install BiocManager first (bigsnpr dependency)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(ask = FALSE, update = FALSE)

# Install bigsnpr >= 1.7.1 to match Ldpred2LDRefCache.R requirement
install.packages("bigsnpr", Ncpus = max(1, parallel::detectCores() - 1))

# Verify version meets minimum requirement
ver <- packageVersion("bigsnpr")
cat(sprintf("\nbigsnpr installed version: %s\n", ver))
if (ver < "1.7.1") {
    stop(sprintf("bigsnpr version %s is below required 1.7.1", ver))
} else {
    cat("bigsnpr version OK.\n")
}

# Pre-load to catch any runtime dependency issues early
library(bigsnpr)
cat("bigsnpr loaded successfully.\n")
