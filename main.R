## ----echo=FALSE----------------------------------------------------------
library(Hmisc)
library(xtable)
library(memisc)
library(rms)
library(reshape2)
library(ggplot2)
library(plyr)
theme_set(theme_bw())

## ----01-load-------------------------------------------------------------
raw <- textConnection("
100 P 18 100 P 14 100 D 23 100 D 18 100 P 10 100 P 17 100 D 18 100 D 22
100 P 13 100 P 12 100 D 28 100 D 21 100 P 11 100 P  6 100 D 11 100 D 25
100 P  7 100 P 10 100 D 29 100 P 12 100 P 12 100 P 10 100 D 18 100 D 14
101 P 18 101 P 15 101 D 12 101 D 17 101 P 17 101 P 13 101 D 14 101 D  7
101 P 18 101 P 19 101 D 11 101 D  9 101 P 12 101 D 11 102 P 18 102 P 15
102 P 12 102 P 18 102 D 20 102 D 18 102 P 14 102 P 12 102 D 23 102 D 19
102 P 11 102 P 10 102 D 22 102 D 22 102 P 19 102 P 13 102 D 18 102 D 24
102 P 13 102 P  6 102 D 18 102 D 26 102 P 11 102 P 16 102 D 16 102 D 17
102 D  7 102 D 19 102 D 23 102 D 12 103 P 16 103 P 11 103 D 11 103 D 25
103 P  8 103 P 15 103 D 28 103 D 22 103 P 16 103 P 17 103 D 23 103 D 18
103 P 11 103 P -2 103 D 15 103 D 28 103 P 19 103 P 21 103 D 17 104 D 13
104 P 12 104 P  6 104 D 19 104 D 23 104 P 11 104 P 20 104 D 21 104 D 25
104 P  9 104 P  4 104 D 25 104 D 19
")
d <- scan(raw, what = "character")
rm(raw)
d <- as.data.frame(matrix(d, ncol = 3, byrow = TRUE))
names(d) <- c("center", "drug", "change")
d$change <- as.numeric(as.character(d$change))
d$drug <- relevel(d$drug, ref = "P")

## ----hamd-xyplot, fig.cap="Distribution of change scores in each centre", fig.width=6, fig.height=3----
p <- ggplot(data = d, aes(x = drug, y = change))
p <- p + geom_jitter(width = .2)
p <- p + geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey30")
p + facet_grid(~ center) + labs(x = "Drug type", y = "HAMD17 change") 

## ------------------------------------------------------------------------
fm <- change ~ drug + center
s <- summary(fm, data = subset(d, center %in% c("100","101","102")), 
             method = "cross", fun = smean.sd)

## ----echo=FALSE, results="asis"------------------------------------------
latex(s, file = "", title = "", caption = "Mean HAMD17 change by drug, center", 
      first.hline.double = FALSE, where = "!htbp", label = "tab:hamd-desc",
      insert.bottom = "Only 3 out of 5 centres are shown.", 
      table.env = TRUE, ctable = TRUE, size = "small", digits = 2)

## ----hamd-delta, fig.cap="Average difference between drug and placebo in each centre", fig.width=6, fig.height=3----
r <- ddply(d, "center", summarize, 
           delta = mean(change[drug == "D"]) - mean(change[drug == "P"]))
p <- ggplot(data = r, aes(x = center, y = delta))
p <- p + geom_point() + geom_hline(yintercept = 0, linetype = 2, colour = "grey30")
p + labs(x = "Center", y = "Difference D-P")

## ------------------------------------------------------------------------
fm <- change ~ drug * center

## ------------------------------------------------------------------------
replications(change ~ drug:center, data = d)

## ------------------------------------------------------------------------
options(contrasts = c("contr.sum", "contr.poly"))
m <- lm(fm, data = d)
anova(m)

## ------------------------------------------------------------------------
car::Anova(m, type = "II")

## ------------------------------------------------------------------------
car::Anova(m, type ="III")

## ------------------------------------------------------------------------
drop1(m, scope = ~ ., test = "F")

## ----echo=FALSE----------------------------------------------------------
print(xtable(anova(m)), file = "s1.tex", floating = FALSE, booktabs = TRUE)
print(xtable(car::Anova(m, type = "II")), file = "s2.tex", floating = FALSE, booktabs = TRUE)
print(xtable(car::Anova(m, type = "III")[-1,]), file = "s3.tex", floating = FALSE, booktabs = TRUE)

## ------------------------------------------------------------------------
D <- model.matrix(m)                            ## design matrix
bhat <- solve(t(D) %*% D) %*% t(D) %*% d$change ## beta parameters
get.ss <- function(C) {
  require(MASS)
  teta <- C%*%bhat
  M <- C %*% ginv(t(D)%*%D) %*% t(C)
  SSH <- t(teta) %*% ginv(M) %*% teta
  return(as.numeric(SSH))
}
## SS(drug|center,drug:center)
get.ss(matrix(c(0,1,0,0,0,0,0,0,0,0), nrow = 1, ncol = 10)) 

## ----02-load-------------------------------------------------------------
source("./urininc.R")
str(d)

## ------------------------------------------------------------------------
s <- summary(change ~ group + strata, data = d, method = "cross", overall = FALSE)

## ----echo=FALSE, results="asis"------------------------------------------
latex(s, file = "", title = "", caption = "Mean change in number of incontinence episods by drug, strata", 
      first.hline.double = FALSE, where = "!htbp", label = "tab:urininc-desc",
      table.env = TRUE, ctable = TRUE, size = "small", digits = 3)

## ----urininc-density, fig.cap="Density estimates for the percent change in frequency of incontinence episodes", fig.width=9, fig.height=3----
p <- ggplot(data = d, aes(x = change, colour = group))
p <- p + geom_line(stat = "density", adjust = 1.2) + facet_grid(~ strata)
p + scale_x_continuous(limits = c(-100, 150)) + labs(x = "Percent change", y = "Density")

## ------------------------------------------------------------------------
library(coin)
dc <- subset(d, complete.cases(d))
independence_test(change ~ group | strata, data = dc, 
                  ytrafo = function(data) trafo(data, numeric_trafo = rank, 
                                                block = dc$strata), 
                  teststat = "quad")

## ------------------------------------------------------------------------
m <- lm(change ~ group + strata, data = d)
car::Anova(m, type = "III")

## ----03-load-------------------------------------------------------------
varnames <- list(strata = 1:4,
                 status = c("Dead", "Alive", "Total"),
                 group = c("Experimental", "Placebo"))
                 
d <- array(c(33,49,48,80,185,169,156,130,218,218,204,210,
             26,57,58,118,189,165,104,123,215,222,162,241),
           dim = c(4,3,2), dimnames = varnames)

## ----eval=FALSE----------------------------------------------------------
## addmargins(d[,-3,], c(1,2))

## ------------------------------------------------------------------------
d <- d[,-3,]
dim(d)

## ------------------------------------------------------------------------
ftable(d)

## ----echo=FALSE, results="asis"------------------------------------------
toLatex(ftable(d, row.vars = 1, col.vars = c(3,2)))

## ----sepsis-dotplot, fig.cap="Proportion of patients who died by the end of the study", fig.width=6, fig.height=6----
dd <- as.data.frame(ftable(d))
r <- ddply(dd, c("strata", "group"), mutate, prop = Freq/sum(Freq))
p <- ggplot(subset(r, status == "Dead"), aes(x = prop, y = group))
p <- p + geom_point() + facet_wrap(~ strata, nrow = 2)
p + scale_x_continuous(limits = c(0,0.5)) + labs(x = "Proportion deads", y = "")

## ----sepsis-cotab, fig.cap="Conditional association plot", fig.width=8, fig.height=8, out.width=".75\\linewidth"----
library(vcd)
cotabplot(d, 1)

## ------------------------------------------------------------------------
n <- rbind(d[,1:2,1], d[,1:2,2])
rownames(n) <- NULL
n <- as.data.frame(n)
n$strata <- gl(4, 1)
n$group <- gl(2, 4, labels = c("Experimental", "Placebo"))
n$group <- relevel(n$group, ref = "Placebo")
m <- glm(cbind(Dead,Alive) ~ group + strata, data = n, family = binomial,
         contrasts = list(strata = "contr.SAS"))
car::Anova(m, type = "III")

## ------------------------------------------------------------------------
exp(confint(m))

