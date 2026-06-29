invisible(Sys.setenv(LANG = "C.UTF-8"))
invisible(Sys.setlocale("LC_ALL", "C.UTF-8"))
'
Calculate pearson correlation for two variables.

Notes:
    1. Read data from stdin and output results to stdout.
    2. Rows with NA (missing value) will be ignored.

Usage:
    PearsonR2.R -f variable1 -s variable2

Options:
    -f string       variable1
    -s string       variable2
    -h --help       Show this screen.
    --version       Show version.
' -> doc

## Auto-detect and install needed packages.
## By Qiuli Chen
options (warn = -1)
#if (!require("pacman")) install.packages("pacman")
suppressMessages(library(pacman))
pacman::p_load(docopt,data.table,pROC,dplyr,rsq)

opts <- docopt(doc,version="V1_2024.04.28")

#options(digits=5)
df <- na.omit(read.table(file("stdin"),header = T,check.names=F))
v1 <- opts$f
v2 <- opts$s

pr <- cor(df[v1],df[v2],method="pearson")
pr <- formatC(pr, digits = 6, format = "g")
cat(c('Pearson_Correlation:', pr, '\n'),sep='\t')
