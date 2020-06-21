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
suppressWarnings(library(ALDEx2))

# glm analysis ------------------------------------------------------------
condsTotal <- ncol(t(otu))
condsNumSplit <- length(unique(map[,condition]))
condsSplitSize <- condsTotal/condsNumSplit
if (test == 'glm'){
covariates <- data.frame("A" = sample(0:1, condsTotal, replace = TRUE),"B" = c(rep(0, condsSplitSize), rep(1, condsSplitSize)))
mm <- model.matrix(~ A + B, covariates)
x <- aldex.clr(t(otu), mm, mc.samples=1, denom="all")
x.tt <- aldex.glm(x)
covariates_effect <- data.frame("A" = sample(0:1, condsTotal, replace = TRUE),"B" = c(rep(0, condsSplitSize), rep(1, condsSplitSize)),"Z" = sample(c(1,2,3), condsTotal, replace=TRUE))
mm_effect <- model.matrix(~ A + Z + B, covariates_effect)
x_effect <- aldex.clr(t(otu), mm_effect, mc.samples=8, denom="all")
x.effect <- aldex.glm.effect(x_effect)
fit <- data.frame(x.effect, x.tt, check.names=F)
}

# t-test analysis ---------------------------------------------------------
else{
fit <- aldex(t(otu), as.character(map[[condition]]),
	     denom=denom, test=test, mc.samples=mc.samples)
}

# final processing --------------------------------------------------------
sfit <- as.data.frame(fit)
write.csv(sfit, file=output)
