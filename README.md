## Classical analysis of clinical trials and medical data using R

Various examples taken from SAS and reproduced in R.

Be aware this is **raw and unedited** alpha release.

### Overview

Statistical analyses are based on the following textbook:

> A Dmitrienko, G Molenberghs, C Chuang-Stein, and W Offen. [Analysis of Clinical Trials Using SAS: A Practical Guide](https://www.sas.com/store/books/categories/usage-and-reference/analysis-of-clinical-trials-using-sas-a-practical-guide/prodBK_59390_en.html). SAS Institute Inc., Cary, NC, USA, 2005.

This book covers various topics: analysis of efficacy and safety study with continuous, categorical or censored endpoints, interim analysis, incomplete data, and multiple testing issues (including gatekeeping strategies). I may (or may not) include extra material from additional textbooks like [Common Statistical Methods for Clinical Research with SAS](https://www.sas.com/store/books/categories/usage-and-reference/common-statistical-methods-for-clinical-research-with-sas-examples-third-edition/prodBK_62004_en.html), [Modern Approaches to Clinical Trials Using SAS](https://www.sas.com/store/books/categories/usage-and-reference/modern-approaches-to-clinical-trials-using-sas-classical-adaptive-and-bayesian-methods/prodBK_67984_en.html), or [Pharmaceutical Statistics Using SAS](https://www.sas.com/store/books/categories/usage-and-reference/pharmaceutical-statistics-using-sas-a-practical-guide/prodBK_60622_en.html).

All analyses were carried out using the latest version of R (3.2.* series). The following R packages are used throughout the document: [ggplot2](http://ggplot2.org), [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html), [plyr](http://plyr.had.co.nz), [Hmisc](http://biostat.mc.vanderbilt.edu/wiki/Main/Hmisc), [rms](https://cran.r-project.org/web/packages/rms/index.html).

In addition to `Hmisc::latex`, the [xtable](https://cran.r-project.org/web/packages/xtable/index.html) package is occasionally used to produce nice looking tables.

The [R code](main.R) is available as a standalone file in the main repository.


### About the PDF

The final output is currently generated using the simplest Makefile I could imagine and [knitr](http://yihui.name/knitr/) which does all the work in the background. For the moment, a classic $\LaTeX$ template without much sophistication is used. I may update the default layout at some point (because nobody likes default templates).
