/* HAMD17 study
 * This is a multicenter clinical trial comparing experimental drug vs. placebo
 * in patients with major depression disorder. The outcome is the change from
 * baseline after 9 weeks of acute treatment, and efficacy is measured using the
 * total score of the Hamilton depression rating scale (17 items).
 */

clear all

import delimited "HAMD17.csv"

encode drug, gen(drug_)
drop drug
rename drug_ drug

label variable drug "Drug type"
label variable change "HAMD17 change"

/* stripplot (N. Cox) available on SSC */
stripplot change, over(drug) by(center) vertical jitter(1 0)

table drug center, contents(count change mean change sd change) column format(%5.1f)


anova change drug##center, sequential

anova change drug##center
