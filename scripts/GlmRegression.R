invisible(Sys.setenv(LANG = "C.UTF-8"))
invisible(Sys.setlocale("LC_ALL", "C.UTF-8"))
'
Regression analysis using glm.

Notes:
    1. Read data from stdin and output results to stdout.
    2. Rows with NA (missing value) will be ignored.
    3. A glm model was fit, supported "family": https://www.statmethods.net/advstats/glm.html

Usage:
    GlmReression.R -f regression_formula -m family_link_function [-n base_model_formula] [-a AUC] [-i beta_95%CI] [-r ModelFitRsquare] [-t PhenotypeName]

Options:
    -f string       Regression formular of Full model.eg:CAD~PRS+age+sex+PC.
    -m string       Family regression type that glm supports: binomial,gaussian,poisson ect.
    -n string       Regression formular of base model.eg:CAD~age+sex+PC
    -a string       Calculate model AUC, can use if predict binary variable outcome.
    -t string       Phenotype name, necessary if -a is provided.
    -i string       Calculate coefficients 95% CI.
    -r string       Output the model fit r square. Provide anything if "yes".
    -h --help       Show this screen.
    --version       Show version.
' -> doc

## Auto-detect and install needed packages.
## By Qiuli Chen
options (warn = -1)
#if (!require("pacman")) install.packages("pacman")
suppressMessages(library(pacman))
pacman::p_load(docopt,data.table,pROC,dplyr,rsq)

opts <- docopt(doc,version="V1_2023.11.02")

#options(digits=5)

outAUC <- opts$a
outCI <- opts$i
outR <- opts$r
N <- opts$n
t <- opts$t

df <- na.omit(read.table(file("stdin"),header = T,check.names=F))

reg.formula <- opts$f

if (!"adjNormPRS" %in% colnames(df)) {
    cat("(This parameter is not available)\n")
    quit(save = "no", status = 0)
}

target.values <- df[["adjNormPRS"]]
target.values <- target.values[!is.na(target.values)]

if (
    length(target.values) == 0 ||
    length(unique(target.values)) <= 1 ||
    all(target.values == 0)
) {
    cat("(It is not available, only zero value, skip)\n")
    quit(save = "no", status = 0)
}

cat("\n")
cat("=========================================================================\n")
cat("Regression Formula\n")
cat("=========================================================================\n")
cat(reg.formula,"\n")



glm.fit <- glm(as.formula(reg.formula),family=opts$m,data=df)
null.glm.fit <- glm(as.formula(N),family=opts$m,data=df)

# Output regression results  
#can output real p value.

cat("\n")
cat("=========================================================================\n")
cat("Regression Results\n")
cat("=========================================================================\n")

summary(glm.fit)$coefficients

# Calculate 95% CI
if(is.null(outCI) == F){
    x= confint(glm.fit)
    y = as.data.frame(x)
    y$CI='95CI'
    y
}

# Calculate AUC
if(is.null(outAUC) == F){
    pred_val <- predict(glm.fit, type='response') 
    roc_obj <- roc(df[[t]], pred_val)
    myauc <- auc(roc_obj)
    ci_auc <- ci.auc(roc_obj)
    cat("=========================================================================\n")
    cat("Model AUC\n")
    cat("=========================================================================\n")
    cat(c('AUC: ', myauc, '\n'),sep='\t')
    cat('AUC_95%_CI:', ci_auc[1],'-',ci_auc[3], '\n')
}

# Calculate model fit correlation
if(is.null(outR) == F){
    model_r2 <- rsq(glm.fit,adj=F)
    model_r2 <- formatC(model_r2, digits = 4, format = "g")

    cat("\n")
    cat("=========================================================================\n")
    cat("Model R Square\n")
    cat("=========================================================================\n")

    cat(c('Full_Model_Rsq: ', model_r2, '\n'),sep='\t')
}

# Calculate R2: R2(full_model) - R2(base_model)
if(is.null(N) == F){
    model_r2_null <- rsq(null.glm.fit,adj=F)
    model_r2_null <- formatC(model_r2_null, digits = 4, format = "g")
    cat(c('Base_Model_Rsq: ', model_r2_null, '\n'),sep='\t')
    diff <- as.numeric(model_r2)-as.numeric(model_r2_null)
    diff <- formatC(diff, digits = 6, format = "g")
    cat(c('Incremental_Model_Rsq: ', diff, '\n'),sep='\t')
}
