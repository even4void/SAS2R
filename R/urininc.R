## Time-stamp: <2016-05-18 17:25:04 chl>

raw <- textConnection("
Placebo  1  -86  -38   43 -100  289    0  -78   38  -80  -25
Placebo  1 -100 -100  -50   25 -100 -100  -67    0  400 -100
Placebo  1  -63  -70  -83  -67  -33    0  -13 -100    0   -3
Placebo  1  -62  -29  -50 -100    0 -100  -60  -40  -44  -14
Placebo  2  -36  -77   -6  -85   29  -17  -53   18  -62  -93
Placebo  2   64  -29  100   31   -6 -100  -30   11  -52  -55
Placebo  2 -100  -82  -85  -36  -75   -8  -75  -42  122  -30
Placebo  2   22  -82    .    .    .    .    .    .    .    .
Placebo  3   12  -68 -100   95  -43  -17  -87  -66   -8   64
Placebo  3   61  -41  -73  -42  -32   12  -69   81    0   87
Drug     1   50 -100  -80  -57  -44  340 -100 -100  -25  -74
Drug     1    0   43 -100 -100 -100 -100  -63 -100 -100 -100
Drug     1 -100 -100    0 -100  -50    0    0  -83  369  -50
Drug     1  -33  -50  -33  -67   25  390  -50    0 -100    .
Drug     2  -93  -55  -73  -25   31    8  -92  -91  -89  -67
Drug     2  -25  -61  -47  -75  -94 -100  -69  -92 -100  -35
Drug     2 -100  -82  -31  -29 -100  -14  -55   31  -40 -100
Drug     2  -82  131  -60     .   .    .    .    .    .    .
Drug     3  -17  -13  -55  -85  -68  -87  -42   36  -44  -98
Drug     3  -75  -35    7  -57  -92  -78  -69   -21 -14    .
")
d <- scan(raw, what = "character")
d <- matrix(d, ncol = 12, byrow = TRUE)
tmp <- as.numeric(as.character(as.vector(t(d[,-c(1,2)]))))
d <- data.frame(group = rep(c("Placebo", "Drug"), c(100, 100)),
                strata = rep(c(1:3,1:3), c(40,40,20,40,40,20)))
d$change <- tmp
rm(raw,tmp)
d$group <- relevel(d$group, ref = "Placebo")
d$strata <- factor(d$strata)
