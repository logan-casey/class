clear all
set more off

use "../data/main_data.dta", clear

capture confirm variable regflag_0203
if _rc {
    di as err "Missing regflag_0203. Re-run clean_data.do."
    exit 111
}
capture confirm variable regflag_2010
if _rc {
    di as err "Missing regflag_2010. Re-run clean_data.do."
    exit 111
}

xtset gvkey fyear
cap mkdir "../output"

* Assignment variable used in regressions
gen v_0203 = L.assignment_v if regflag_0203 == 1
bysort gvkey: egen v_2009_val = max(cond(fyear == 2009, assignment_v, .))
gen v_2010 = v_2009_val if regflag_2010 == 1

* Display counts outside plotting range [-15, 15] $M
count if regflag_0203 == 1 & !missing(v_0203)
display as result "2002 sample nonmissing V: " r(N)
count if regflag_0203 == 1 & !missing(v_0203) & (v_0203 < -15 | v_0203 > 15)
display as result "2002 sample outside [-15,15] $M: " r(N)

count if regflag_2010 == 1 & !missing(v_2010)
display as result "2009 sample nonmissing V: " r(N)
count if regflag_2010 == 1 & !missing(v_2010) & (v_2010 < -15 | v_2010 > 15)
display as result "2009 sample outside [-15,15] $M: " r(N)

* Left panel: 2002 policy period sample (outcome years 2002/2003)
histogram v_0203 if regflag_0203 == 1 & !missing(v_0203) & inrange(v_0203, -15, 15), ///
    start(-15) ///
    width(0.25) frequency ///
    xline(0, lcolor(red) lpattern(dash)) ///
    xlabel(-15(5)15) ///
    xtitle("Assignment variable V ($M)") ///
    ytitle("Number of firms") ///
    title("2002 Policy Period") ///
    name(fig3_left, replace)

* Right panel: 2009 policy period sample (outcome years 2010/2011)
histogram v_2010 if regflag_2010 == 1 & !missing(v_2010) & inrange(v_2010, -15, 15), ///
    start(-15) ///
    width(0.25) frequency ///
    xline(0, lcolor(red) lpattern(dash)) ///
    xlabel(-15(5)15) ///
    xtitle("Assignment variable V ($M)") ///
    ytitle("Number of firms") ///
    title("2009 Policy Period") ///
    name(fig3_right, replace)

graph combine fig3_left fig3_right, col(2) ///
    title("Figure 3: Distribution of Sample Firms around Kink (V=0)")

graph export "../output/figure3_distribution.png", width(2200) replace
