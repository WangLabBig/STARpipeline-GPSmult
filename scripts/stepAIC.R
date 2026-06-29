invisible(Sys.setenv(LANG = "C.UTF-8"))
invisible(Sys.setlocale("LC_ALL", "C.UTF-8"))
'
Feature selection by stepAIC.

Notes:
    1. Read data from stdin and output results to stdout.
    2. A glm model was fit, supported "family": https://www.statmethods.net/advstats/glm.html

Usage:
    stepAIC.R -f regression_formula -m family_link_function

Options:
    -f string       Regression formula.eg:CAD~PRS+age+sex+PC.
    -m string       Family regression type that glm supports: binomial,gaussian,poisson ect.
    -h --help       Show this screen.
    --version       Show version.
' -> doc

## Auto-detect and install needed packages.
## Written by someone
options (warn = -1)
#if (!require("pacman")) install.packages("pacman")
suppressMessages(library(pacman))
pacman::p_load(docopt,data.table,MASS,dplyr)
opts <- docopt(doc,version="V1_2023.11.08")
df <- read.table(file("stdin"),header = T,check.names=F)
reg.formula <- opts$f
print("Regression formula:")
print(reg.formula)
#print("-------------------------")
#print("Summary of full model:")
full.model <- glm(as.formula(reg.formula),family=opts$m,data=df)
#summary(full.model)

# StepWise Regression
step.model <- stepAIC(full.model, direction = "both",
                      trace = FALSE)
#print("-------------------------")
print("Summary of stepwise model:")
summary(step.model)
