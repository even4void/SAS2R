/* Severe Sepsis study
 * This is a placebo-controlled RCT examining the effect of an experimental
 * drug on 28-day all-cause mortality in patients with severe sepsis. Patients
 * were allocated to one of four strata depending on their APACHE II score.
 */

clear all

/* maybe set your working directory near hear */
import delimited "sepsis.csv"

encode group, gen(group_)
drop group
rename group_ group

gen total = dead+alive
order total, after(alive)

/* As of Stata 14, we should better use glm instead of blogit. */
blogit dead total group ib4.strata
