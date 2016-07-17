## ----echo=FALSE, cache=FALSE----------------------------------------------------------------------
library(Hmisc)  # also lattice, ggplot2
library(xtable)
library(memisc)
library(rms)
library(multcomp)
library(reshape2)
library(plyr)
library(hrbrthemes)
theme_set(theme_ipsum(base_family = "Bitstream Vera Sans", base_size = 11))
options(digits = 6, show.signif.stars = FALSE, width = 100)


## ----01-load--------------------------------------------------------------------------------------
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


## ----hamd-xyplot, fig.cap="Distribution of change scores in each centre", fig.width=8, fig.height=4----
p <- ggplot(data = d, aes(x = drug, y = change)) +
     geom_jitter(width = .2, color = grey(.3)) +
     geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "lightcoral") +
     facet_grid(~ center) +
     labs(x = "Drug type", y = "HAMD17 change")
p


## -------------------------------------------------------------------------------------------------
fm <- change ~ drug + center
s <- summary(fm, data = subset(d, center %in% c("100","101","102")),
             method = "cross", fun = smean.sd)


## ----echo=FALSE, results="asis"-------------------------------------------------------------------
latex(s, file = "", title = "", caption = "Mean HAMD17 change by drug, center",
      first.hline.double = FALSE, where = "!htbp", label = "tab:hamd-desc",
      insert.bottom = "Only 3 out of 5 centres are shown.",
      table.env = TRUE, ctable = TRUE, size = "small", digits = 2)


## ----hamd-delta, fig.cap="Average difference between drug and placebo in each centre", fig.width=8, fig.height=4----
m <- with(d, Hmisc::summarize(change, llist(drug, center), mean))
r <- aggregate(change ~ center, m, diff)
p <- ggplot(data = r, aes(x = center, y = change)) +
     geom_point() +
     geom_hline(yintercept = 0, linetype = 2, colour = grey(.3)) +
     labs(x = "Center", y = "Difference D-P")
p


## -------------------------------------------------------------------------------------------------
fm <- change ~ drug * center
replications(change ~ drug:center, data = d)


## -------------------------------------------------------------------------------------------------
options(contrasts = c("contr.sum", "contr.poly"))
m <- lm(fm, data = d)
anova(m)


## ----eval=FALSE-----------------------------------------------------------------------------------
## anova(ols(fm, d))


## -------------------------------------------------------------------------------------------------
car::Anova(m, type = "II")


## -------------------------------------------------------------------------------------------------
car::Anova(m, type ="III")


## -------------------------------------------------------------------------------------------------
drop1(m, scope = ~ ., test = "F")


## ----echo=FALSE-----------------------------------------------------------------------------------
print(xtable(anova(m)[,c(2,1,4,5)]), file = "s1.tex", floating = FALSE, booktabs = TRUE)
print(xtable(car::Anova(m, type = "II")), file = "s2.tex", floating = FALSE, booktabs = TRUE)
print(xtable(car::Anova(m, type = "III")[-1,]), file = "s3.tex", floating = FALSE, booktabs = TRUE)


## -------------------------------------------------------------------------------------------------
D <- model.matrix(m)                            ## design matrix
bhat <- solve(t(D) %*% D) %*% t(D) %*% d$change ## beta parameters
get.ss <- function(C) {
  require(MASS)
  teta <- C %*% bhat
  M <- C %*% ginv(t(D) %*% D) %*% t(C)
  SSH <- t(teta) %*% ginv(M) %*% teta
  return(as.numeric(SSH))
}
## SS(drug|center,drug:center)
get.ss(matrix(c(0,1,0,0,0,0,0,0,0,0), nrow = 1, ncol = 10))


## -------------------------------------------------------------------------------------------------
summary(glht(m, mcp(drug = "Tukey", interaction_average = TRUE)))


## -------------------------------------------------------------------------------------------------
library(nlme)
m <- lme(change ~ drug, data = d, random = ~ 1 | center/drug)
summary(m)


## -------------------------------------------------------------------------------------------------
library(lme4)
library(lmerTest)
m <- lmer(change ~ drug + (1 | center/drug), data = d)
summary(m, ddf = "Satterthwaite")  # this is the default


## -------------------------------------------------------------------------------------------------
library(QualInt)
with(d, qualint(change, drug, center, test = "LRT"))


## ----hamd-ibga, echo=FALSE, fig.cap="Average differences between drug and placebo stratified by centres", fig.width=6, fig.height=3----
r <- with(d, qualint(change, drug, center, test = "IBGA"))
plot(r) # + scale_colour_ipsum()


## ----02-load--------------------------------------------------------------------------------------
source("R/urininc.R")
str(d)


## -------------------------------------------------------------------------------------------------
s <- summary(change ~ group + strata, data = d, method = "cross", overall = FALSE)


## ----echo=FALSE, results="asis"-------------------------------------------------------------------
latex(s, file = "", title = "", caption = "Mean change in number of incontinence episods by drug, strata",
      first.hline.double = FALSE, where = "!htbp", label = "tab:urininc-desc",
      table.env = TRUE, ctable = TRUE, size = "small", digits = 3)


## ----urininc-density, fig.cap="Density estimates for the percent change in frequency of incontinence episodes", fig.width=9, fig.height=3----
p <- ggplot(data = d, aes(x = change, color = group)) +
     geom_line(stat = "density", adjust = 1.2) +
     facet_wrap(~ strata, ncol = 3) +
     scale_color_manual("Group", values = c("steelblue", "orange")) +
     scale_x_continuous(limits = c(-100, 150)) +
     labs(x = "Percent change", y = "Density")
p


## -------------------------------------------------------------------------------------------------
library(coin)
dc <- subset(d, complete.cases(d))
independence_test(change ~ group | strata, data = dc,
                  ytrafo = function(data) trafo(data, numeric_trafo = rank,
                                                block = dc$strata),
                  teststat = "quad")


## -------------------------------------------------------------------------------------------------
m <- lm(change ~ group + strata, data = d)
car::Anova(m, type = "III")


## ----03-load--------------------------------------------------------------------------------------
varnames <- list(strata = 1:4,
                 status = c("Dead", "Alive", "Total"),
                 group = c("Experimental", "Placebo"))
d <- array(c(33,49,48,80,185,169,156,130,218,218,204,210,
             26,57,58,118,189,165,104,123,215,222,162,241),
           dim = c(4,3,2), dimnames = varnames)


## ----eval=FALSE-----------------------------------------------------------------------------------
## addmargins(d[,-3,], c(1,2))


## -------------------------------------------------------------------------------------------------
d <- d[,-3,]
dim(d)


## -------------------------------------------------------------------------------------------------
ftable(d)


## ----echo=FALSE, results="asis"-------------------------------------------------------------------
toLatex(ftable(d, row.vars = 1, col.vars = c(3,2)), digits = 0)


## ----sepsis-dotplot, fig.cap="Proportion of patients who died by the end of the study", fig.width=8, fig.height=4----
dd <- as.data.frame(ftable(d))
r <- ddply(dd, c("strata", "group"), mutate, prop = Freq/sum(Freq))
p <- ggplot(subset(r, status == "Dead"), aes(x = prop, y = group)) +
     geom_point() + facet_wrap(~ strata, nrow = 2) +
     scale_x_continuous(limits = c(0,0.5)) +
     labs(x = "Proportion deads", y = "")
p


## -------------------------------------------------------------------------------------------------
n <- rbind(d[,1:2,1], d[,1:2,2])
rownames(n) <- NULL
n <- as.data.frame(n)
n$strata <- gl(4, 1)
n$group <- gl(2, 4, labels = c("Experimental", "Placebo"))
n$group <- relevel(n$group, ref = "Placebo")


## -------------------------------------------------------------------------------------------------
m <- glm(cbind(Dead,Alive) ~ group + strata, data = n, family = binomial,
         contrasts = list(strata = "contr.SAS"))
car::Anova(m, type = "III")


## -------------------------------------------------------------------------------------------------
exp(confint(m))


## -------------------------------------------------------------------------------------------------
pvals <- c(0.047, 0.0167, 0.015)  ## scenario 1
p.adjust(pvals, method = "bonferroni")


## -------------------------------------------------------------------------------------------------
f <- function(x) (1-(1-x)^length(x))
f(pvals)


## -------------------------------------------------------------------------------------------------
p.adjust(pvals, method = "holm")


## -------------------------------------------------------------------------------------------------
pvals <- c(0.047, 0.027, 0.015)  ## scenario 2
p.adjust(pvals, method = "holm")
p.adjust(pvals, method = "hommel")


## -------------------------------------------------------------------------------------------------
pvals <- c(0.053, 0.026, 0.017)  ## scenario 3
p.adjust(pvals, method = "hochberg")
p.adjust(pvals, method = "hommel")


## -------------------------------------------------------------------------------------------------
tmp <- matrix(c(0.25,10,0.58,0.29,10,0.71,0.35,
                0.5,10,0.62,0.31,10,0.88,0.33,
                0.75,10,0.51,0.33,10,0.73,0.36,
                1,10,0.34,0.27,10,0.68,0.29,
                2,10,-0.06,0.22,10,0.37,0.25,
                3,10,0.05,0.23,10,0.43,0.28),
                nrow = 6, byrow = TRUE)
colnames(tmp) <- c("time", "N0", "Mean0", "SD0", "N1", "Mean1", "SD1")
d <- melt(as.data.frame(tmp), id.vars = 1, measure.vars = c(3,4,6,7))


## -------------------------------------------------------------------------------------------------
r <- ddply(dcast(d, time ~ variable), "time", mutate,
           diff = Mean1 - Mean0, se = sqrt((1/10+1/10)*(SD0^2+SD1^2)/2))


## ----fev-xyplot, fig.cap="Treatment comparisons in the asthma study", fig.width=8, fig.height=4----
p <- ggplot(r, aes(x = time, y = diff))
p <- p + geom_line() + geom_point()
p <- p + geom_line(aes(x = time, y = diff - qt(0.95, 20-2) * se), linetype = 2)
p <- p + scale_x_continuous(breaks = seq(0, 3, by = 1))
p <- p + scale_y_continuous(breaks = seq(-0.2, 0.5, by = 0.1))
p <- p + geom_hline(aes(yintercept = 0))
p + labs(x = "Time (hours)", y = "Treatment difference (95% Lower CI)")


## ----eval = FALSE---------------------------------------------------------------------------------
## p <- p + geom_errorbar(aes(ymin = diff - qt(0.975, 20-2) * se,
##                            ymax = diff + qt(0.975, 20-2) * se),
##                        width = .1)

