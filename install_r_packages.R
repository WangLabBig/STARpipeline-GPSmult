options(repos = c(
    TUNA = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
    ALIY = "https://mirrors.aliyun.com/CRAN/"
))

# Install in dependency order to isolate failures
pkg_groups <- list(
    # Group 1: base utilities with no tricky dependencies
    base = c(
        "pacman",
        "magrittr",
        "R.utils",
        "optparse",
        "foreach",
        "doParallel",
        "MASS"
    ),
    # Group 2: data manipulation
    data = c(
        "data.table",
        "Matrix",
        "dplyr",
        "tidyr",
        "readr",
        "purrr",
        "tibble",
        "stringr"
    ),
    # Group 3: stats
    stats = c(
        "pROC",
        "docopt"
    ),
    # # Group 4: rio dependencies first, then rio
    rio_deps = c(
        "haven",
        "foreign",
        "readxl",
        "openxlsx"
    )
)

all_pkgs <- unlist(pkg_groups, use.names = FALSE)

for (grp_name in names(pkg_groups)) {
    cat(sprintf("\n--- Installing group: %s ---\n", grp_name))
    install.packages(
        pkg_groups[[grp_name]],
        Ncpus = max(1, parallel::detectCores() - 1)
    )
}

if (!require("remotes")){
    install.packages("remotes")
}
remotes::install_github("gesistsa/rio")

install.packages("curl", repos = "https://jeroen.r-universe.dev")

remove.packages("rlang")
install.packages(
  "rlang",
  dependencies = TRUE
)

remove.packages("vctrs")
install.packages(
  "vctrs",
  dependencies = TRUE
)

install.packages(
  "rsq",
  dependencies = TRUE
)

# Verify all
cat("\n========== Installation check ==========\n")
failed <- c()
for (pkg in all_pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("  OK  : %s\n", pkg))
    } else {
        cat(sprintf("  FAIL: %s\n", pkg))
        failed <- c(failed, pkg)
    }
}

if (length(failed) > 0) {
    stop(sprintf("The following packages failed to install: %s",
                 paste(failed, collapse = ", ")))
} else {
    cat("All packages installed successfully.\n")
}
