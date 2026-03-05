clear all
set more off

use "../data/main_data.dta", clear

capture confirm variable regflag_0203
if _rc {
    di as err "Missing regflag_0203. Re-run clean_data.do."
    exit 111
}
capture confirm variable sample_firm_0203
if _rc {
    di as err "Missing sample_firm_0203. Re-run clean_data.do."
    exit 111
}
capture confirm variable firm_complete_0203
if _rc {
    di as err "Missing firm_complete_0203. Re-run clean_data.do."
    exit 111
}
capture confirm variable firm_has_policy_loss_0203
if _rc {
    di as err "Missing firm_has_policy_loss_0203. Re-run clean_data.do."
    exit 111
}

xtset gvkey fyear

* Outcome-year rows used for 2002/2003 first-stage setup
gen byte row_0203 = inlist(fyear, 2002, 2003)

* Variables aligned to first-stage timing
gen refund_0203 = L.potential_refund if row_0203
gen v_0203      = L.assignment_v     if row_0203
gen byte has_timing_inputs = row_0203 & !missing(refund_0203, v_0203)

* Firm-level restrictions (as implemented)
gen byte pass_complete = (firm_complete_0203 == 1)
gen byte pass_losselig = (firm_has_policy_loss_0203 == 1)
gen byte pass_no_outlier = (sample_firm_0203 == 1) if pass_complete & pass_losselig

display as text "============================================================"
display as text "ATTRITION DECOMPOSITION: 2002/2003 FIRST-STAGE SAMPLE"
display as text "============================================================"

count if row_0203
local n0 = r(N)
display as result "Row base (fyear 2002/2003): " %9.0f `n0'

count if has_timing_inputs
local n1 = r(N)
display as result "After requiring nonmissing L.refund and L.assignment_v: " %9.0f `n1' ///
    "   drop = " %9.0f (`n0' - `n1')

count if has_timing_inputs & pass_complete
local n2 = r(N)
display as result "After firm_complete_0203==1: " %9.0f `n2' ///
    "   drop = " %9.0f (`n1' - `n2')

count if has_timing_inputs & pass_complete & pass_losselig
local n3 = r(N)
display as result "After firm_has_policy_loss_0203==1: " %9.0f `n3' ///
    "   drop = " %9.0f (`n2' - `n3')

count if has_timing_inputs & sample_firm_0203 == 1
local n4 = r(N)
display as result "After sample_firm_0203==1 (includes outlier rule): " %9.0f `n4' ///
    "   drop = " %9.0f (`n3' - `n4')

count if regflag_0203 == 1
local n5 = r(N)
display as result "Final regflag_0203 rows: " %9.0f `n5' ///
    "   delta vs prior = " %9.0f (`n4' - `n5')

display as text "------------------------------------------------------------"
display as text "Firm-level diagnostics"

egen tag_firm = tag(gvkey)

count if tag_firm
local f0 = r(N)
display as result "Unique firms in dataset: " %9.0f `f0'

count if tag_firm & pass_complete
local f1 = r(N)
display as result "Firms passing complete-history rule: " %9.0f `f1'

count if tag_firm & pass_complete & pass_losselig
local f2 = r(N)
display as result "Firms passing complete+loss-eligibility: " %9.0f `f2'

count if tag_firm & sample_firm_0203 == 1
local f3 = r(N)
display as result "Firms in sample_firm_0203: " %9.0f `f3'
display as result "Implied firm drop at outlier stage: " %9.0f (`f2' - `f3')

display as text "------------------------------------------------------------"
display as text "Shares relative to row base"
display as result "Timing-input share: " %6.3f (`n1'/`n0')
display as result "Complete-history share: " %6.3f (`n2'/`n0')
display as result "Loss-eligible share: " %6.3f (`n3'/`n0')
display as result "Final regflag share: " %6.3f (`n5'/`n0')

display as text "============================================================"
display as text "Done."
display as text "============================================================"

