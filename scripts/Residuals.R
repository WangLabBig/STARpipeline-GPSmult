invisible(Sys.setenv(LANG = "C.UTF-8"))
invisible(Sys.setlocale("LC_ALL", "C.UTF-8"))
'
Extract residulas from linear regression.

Notes:
    1. Read data from stdin and output results to stdout.
    2. Remove rows with NA.
    3. Add one column with header as residual.

Usage:
    Residuals.R -f fomular [-t title]

Options:
    -f string   Regression formular.for example: PRS ~ PC1 + PC2.
    -t string   Output residual name.
    -h --help   Show this screen.
    --version   Show version.
' -> doc
## Auto-detect and install needed packages.
## Written by someone
options (warn = -1)
#if (!require("pacman")) install.packages("pacman")
suppressMessages(library(pacman))
pacman::p_load(docopt,data.table)

opts <- docopt(doc,version="V1_2023.02.15")

reg.formula <- opts$f
title  <- if(is.null(opts$t)) 'residual' else opts$t
df = na.omit(read.table(file("stdin"),header = T,check.names=F))

reg.fit <- lm(as.formula(reg.formula),data=df)
df[[title]] <- residuals(reg.fit)

write.table(df ,"",row.names=F,col.names=T,quote=F)
