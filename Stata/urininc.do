/* Urinary incontinence study
 * This is a subset of data collected in an RCT on urinary
 * incontinence where the primary endpoint was the percent change from baseline
 * of number of incontinence episodes per week over an 8-week period. Patients were
 * initially randomized into one of three strata depending on the baseline
 * frequency of incontinence episodes.
 */

clear all

/* maybe set your working directory near hear */
import delimited "urininc.csv"

encode group, gen(drug)

table drug strata, contents(count change mean change sd change) col format(%5.1f)

twoway (kdensity change if drug == 1, bw(50) range(-100 150)) || (kdensity change if drug == 2, bw(50) range(-100 150)), by(strata, cols(3)) xlabel(-100(50)150) legend(label(1 "Drug") label(2 "Placebo"))

/* search alignedranks and then net install */
alignedranks change, by(drug) strata(strata)

/* Type III ANOVA */
anova change drug strata

