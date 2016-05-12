## ----setup, echo=FALSE, include=FALSE------------------------------------
require(knitr)
opts_chunk$set(cache=TRUE, fig.align="center", message=FALSE, warning=FALSE, size="small",
               dev = "cairo_pdf", dev.args = list(family = "Bitstream Vera Sans"))
library(Hmisc)
library(rms)
library(reshape2)
library(ggplot2)
library(plyr)
theme_set(theme_bw())

## ------------------------------------------------------------------------
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

## ----fig.width=6, fig.height=3-------------------------------------------
p <- ggplot(data = d, aes(x = drug, y = change))
p <- p + geom_jitter(width = .2)
p <- p + geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey30")
p + facet_grid(~ center) + labs(x = "Drug type", y = "HAMD17 change") 

## ------------------------------------------------------------------------
fm <- change ~ drug + center
s <- summary(fm, data = subset(d, center %in% c("100","101","102")), 
             method = "cross", fun = smean.sd)

## ----echo=FALSE, results="asis"------------------------------------------
latex(s, file="", title="", caption="Mean HAMD17 change by drug, center", 
      first.hline.double=FALSE, where="!htbp",
      insert.bottom="Only 3 out of 5 centres are shown.", 
      table.env=TRUE, ctable = TRUE, size = "small", digits = 2)

## ----fig.width=6, fig.height=3-------------------------------------------
r <- ddply(d, "center", summarize, 
           delta = mean(change[drug == "D"]) - mean(change[drug == "P"]))
p <- ggplot(data = r, aes(x = center, y = delta))
p <- p + geom_point() + geom_hline(yintercept = 0, linetype = 2, colour = "grey30")
p + labs(x = "Center", y = "Difference D-P")

## ------------------------------------------------------------------------
m <- lm(change ~ drug * center, data = d)
anova(m)

## ------------------------------------------------------------------------
car::Anova(m, type = "II")

## ------------------------------------------------------------------------
car::Anova(m, type ="III")

## ------------------------------------------------------------------------
## options(contrasts = c("contr.sum", "contr.poly"))
drop1(m, scope = ~ ., test = "F")

