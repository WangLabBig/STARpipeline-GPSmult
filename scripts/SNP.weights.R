#!/usr/bin/env Rscript
invisible(Sys.setenv(LANG = "C.UTF-8"))
invisible(Sys.setlocale("LC_ALL", "C.UTF-8"))
suppressWarnings(
suppressMessages(
library(data.table)
)
)

# ============================================================
# Usage:
#
# Rscript SNP.weights.R \
#     snpid.txt \
#     PRS.weights.txt \
#     output.txt
#
# PRS.weights.txt:
#
# BMI.AMR.HISLA     0.0326789
# BMI.AMR.MXBB      0.0733579
# T2D.EUR.GERA      0.00970096
#
# Column 1 = PRS file name (or prefix)
# Column 2 = mixing weight
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3){
    stop(
        paste(
            "Usage:",
            "Rscript SNP.weights.R",
            "snpid.txt",
            "PRS.weights.txt",
            "output.txt"
        )
    )
}

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

snp_file    <- args[1]
weight_file <- args[2]
out_file    <- args[3]

# ------------------------------------------------------------
# Read PRS-weight table
# ------------------------------------------------------------

weight_dt <- fread(
    weight_file,
    header = FALSE
)

if(ncol(weight_dt) < 2){
    stop(
        "PRS.weights.txt must contain at least 2 columns"
    )
}

prs_files <- weight_dt[[1]]
weights   <- suppressWarnings(
    as.numeric(weight_dt[[2]])
)

if(any(is.na(weights))){
    print(weight_dt)
    stop(
        "Some weights could not be converted to numeric."
    )
}

# automatically append .txt if missing

prs_files <- ifelse(
    grepl("\\.txt$", prs_files),
    prs_files,
    paste0(prs_files, ".txt")
)

# ------------------------------------------------------------
# Check PRS files
# ------------------------------------------------------------

missing_files <- prs_files[
    !file.exists(prs_files)
]

if(length(missing_files) > 0){
    stop(
        paste(
            "Missing files:",
            paste(missing_files, collapse = ", ")
        )
    )
}

# ------------------------------------------------------------
# Read reference SNP list
# ------------------------------------------------------------

ref <- fread(
    snp_file,
    header = TRUE
)

setnames(
    ref,
    names(ref)[1:2],
    c("rsID", "a1")
)

if(anyDuplicated(ref$rsID)){
    stop(
        "Duplicate rsID found in reference SNP file."
    )
}

# ------------------------------------------------------------
# Harmonize one PRS file
# ------------------------------------------------------------

process_prs <- function(prs_file, ref){

    prs <- fread(
        prs_file,
        header = TRUE
    )

    if(ncol(prs) < 3){
        stop(
            paste(
                prs_file,
                "must contain at least 3 columns"
            )
        )
    }

    setnames(
        prs,
        names(prs)[1:3],
        c("rsID", "a1", "beta")
    )

    if(anyDuplicated(prs$rsID)){
        stop(
            paste(
                "Duplicate rsID found in",
                prs_file
            )
        )
    }

    merged <- merge(
        ref,
        prs,
        by = "rsID",
        all.x = TRUE,
        sort = FALSE,
        suffixes = c("_ref", "_prs")
    )

    # Verify SNP order
    if(!identical(merged$rsID, ref$rsID)){
        stop(
            paste(
                "SNP order mismatch detected in",
                prs_file
            )
        )
    }

    # Missing SNPs -> beta = 0
    merged[
        is.na(beta),
        beta := 0
    ]

    # Flip beta if allele orientation differs
    merged[
        !is.na(a1_prs) &
        a1_ref != a1_prs &
        beta != 0,
        beta := -beta
    ]

    merged[
        ,
        .(
            rsID,
            beta
        )
    ]
}

# ------------------------------------------------------------
# Process all PRS files
# ------------------------------------------------------------

message(
    "Reading ",
    length(prs_files),
    " PRS files..."
)

prs_list <- lapply(
    prs_files,
    process_prs,
    ref = ref
)

# ------------------------------------------------------------
# Weighted sum
# ------------------------------------------------------------

final <- copy(ref)

final[
    ,
    final_beta := 0
]

for(i in seq_along(prs_list)){

    message(
        "Adding ",
        basename(prs_files[i]),
        " (weight = ",
        signif(weights[i], 6),
        ")"
    )

    final[
        ,
        final_beta :=
            final_beta +
            prs_list[[i]]$beta * weights[i]
    ]
}

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

fwrite(
    final[
        ,
        .(
            rsID,
            a1,
            final_beta
        )
    ],
    file = out_file,
    sep = "\t",
    quote = FALSE
)

message("")
message("Finished!")
message("Output: ", out_file)
message("Number of SNPs: ", nrow(final))
message("Number of PRS files: ", length(prs_files))
