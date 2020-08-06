#!/usr/bin/env Rscript

# load arguments ---------------------------------------------------------
cat(R.version$version.string, "\n")

args <- commandArgs(TRUE)
inp.abundances.path <- args[[1]]
inp.metadata.path <- args[[2]]
condition <- args[[3]] # condition(s)
mc.samples <- as.integer(args[[4]]) # mc.samples
test <- args[[5]] # test
denom <- args[[6]] # denom
output <- args[[7]] # output

# load data ---------------------------------------------------------------
map <- read.delim(inp.metadata.path, check.names=FALSE, row.names=1) # conditions
otu <- read.delim(inp.abundances.path, check.names=FALSE, row.names=1) # reads

# load libraries ----------------------------------------------------------

source("/Users/Tyler/Desktop/ALDEx_bioc/R/AllClasses.R")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/AllGenerics.R")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/aldex.R")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/clr_corr.R")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/clr_effect.r")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/clr_function.r")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/clr_glm-tpq.R")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/clr_glm.r")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/clr_ttest.r")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/expectedDistance.aldex.r")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/glm_effect.r")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/iqlr_features.r")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/plot.aldex.r")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/plotFeature.aldex.r")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/progress.R")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/rdirichlet.r")
source("/Users/Tyler/Desktop/ALDEx_bioc/R/stats.fast.R")

#suppressWarnings(library(ALDEx2)) # <- NEED TO TEST THIS (Rather than use local directory sources above) *********************

# glm (if) or t-test (else) analysis --------------------------------------
if(test == 'glm'){
    #new glm addition (testing variations of conds)
    condsTotal <- ncol(t(otu))
    condsNumSplit <- length(unique(map[,condition]))
    condsSplitSize <- condsTotal/condsNumSplit
    covariates <- data.frame("A" = sample(0:1, condsTotal, replace = TRUE), "B" = c(rep(0, condsSplitSize), rep(1, condsSplitSize)))

    mm <- model.matrix(~ A + B, covariates)
    x <- aldex.clr(t(otu), mm, mc.samples=1, denom=denom)
    xglm <- aldex.glm(x)
    covariates_effect <- data.frame("A" = sample(0:1, condsTotal, replace = TRUE),"B" = c(rep(0, condsSplitSize), rep(1, condsSplitSize))) #,"Z" = sample(c(1,2,3), 8, replace=TRUE))
    mm_effect <- model.matrix(~ A + B, covariates_effect) # (~ A + Z + B, covariates_effect)                                               # DO I NEED THE Z FOR THE TEST????
    x_effect <- aldex.clr(t(otu), mm_effect, mc.samples=8, denom=denom)
    #x_effect <- aldex.clr(t(otu), mm, mc.samples=8, denom=denom)
    xeffect <- aldex.glm.effect(x_effect)
    fit <- data.frame(xeffect, xglm, check.names=F)
} else {
    fit <- aldex(t(otu), as.character(map[[condition]]), denom=denom, test=test, mc.samples=mc.samples)
}

# final processing --------------------------------------------------------
sfit <- as.data.frame(fit)
write.csv(sfit, file=output)
