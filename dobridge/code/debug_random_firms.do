clear all
set more off

* Usage:
*   do debug_random_firms.do [nfirms] [seed]
* Example:
*   do debug_random_firms.do 20 20260303

args nfirms seed
if "`nfirms'" == "" local nfirms = 10
if "`seed'"   == "" local seed   = 20260303

set seed `seed'
use "../data/main_data.dta", clear

display as text "==============================================="
display as text "Random Firm Debug Report"
display as text "nfirms = `nfirms' | seed = `seed'"
display as text "==============================================="

* Basic integrity checks
capture noisily isid gvkey fyear
if _rc {
    display as error "WARNING: (gvkey, fyear) is not unique."
}

* Random firm draw (one random number per firm)
bysort gvkey: gen byte firm_tag = (_n == 1)
gen double u = runiform() if firm_tag
bysort gvkey: egen double ufirm = max(u)

preserve
    keep if firm_tag
    gsort ufirm
    keep in 1/`nfirms'
    keep gvkey
    tempfile draw_firms
    save `draw_firms'
restore

merge m:1 gvkey using `draw_firms', keep(match) nogen

display as text " "
display as text "---- Selected firms ----"
levelsof gvkey, local(sample_gvkeys)
display as result "`sample_gvkeys'"

* Firm-level panel diagnostics
display as text " "
display as text "---- Firm-level panel diagnostics ----"
gen byte policy_obs_flag    = inlist(fyear,2001,2002,2008,2009)
gen byte refund_nonmiss_flag = !missing(potential_refund)
gen byte assign_nonmiss_flag = !missing(assignment_v)

preserve
    collapse ///
        (count) n_obs=fyear ///
        (min) min_year=fyear ///
        (max) max_year=fyear ///
        (sum) n_policy_obs=policy_obs_flag ///
        (sum) n_refund_nonmiss=refund_nonmiss_flag ///
        (sum) n_assign_nonmiss=assign_nonmiss_flag, by(gvkey)
    list, sep(0) noobs abbrev(20)
restore

* Window coverage checks by firm (useful for carryback construction)
display as text " "
display as text "---- Required-year coverage checks ----"
gen byte w_9602 = inrange(fyear, 1996, 2002)
gen byte w_0309 = inrange(fyear, 2003, 2009)
bysort gvkey: egen n_9602 = total(w_9602)
bysort gvkey: egen n_0309 = total(w_0309)

preserve
    keep gvkey n_9602 n_0309 firm_complete_0203 firm_complete_2010
    bysort gvkey: keep if _n == 1
    list gvkey n_9602 n_0309 firm_complete_0203 firm_complete_2010, sep(0) noobs
restore

* Detailed line-by-line panel view for carryback logic
display as text " "
display as text "---- Detailed yearly values (1996-2009) ----"
sort gvkey fyear
list gvkey fyear at taxinc profit loss taxrate ///
     loss_2002 loss_2003 loss_applied_2010 ///
     assign_v_2002 assign_v_2003 assign_v_2010 assignment_v ///
     refund_2002 refund_2003 refund_2010 potential_refund ///
     used_2008_loss regflag_0203 regflag_2010 ///
     if inrange(fyear, 1996, 2009), sepby(gvkey) noobs abbrev(20)

* Outcome-year mapping check for first-stage use
display as text " "
display as text "---- Outcome-year mapping check (t-1 policy values) ----"
xtset gvkey fyear
gen refund_lag = L.potential_refund
gen v_lag      = L.assignment_v
list gvkey fyear refund_lag v_lag regflag_0203 regflag_2010 ///
     if inlist(fyear, 2002, 2003, 2010, 2011), sepby(gvkey) noobs abbrev(20)

display as text "==============================================="
display as text "End Random Firm Debug Report"
display as text "==============================================="
