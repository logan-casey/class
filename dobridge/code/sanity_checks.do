clear all
set more off

use "../data/main_data.dta", clear

capture confirm variable regflag_0203
if _rc {
    display as error "WARNING: regflag_0203 not found. Re-run clean_data.do."
    exit 111
}
capture confirm variable regflag_2010
if _rc {
    display as error "WARNING: regflag_2010 not found. Re-run clean_data.do."
    exit 111
}

display as text "=============================="
display as text "SANITY CHECKS: DATA INTEGRITY"
display as text "=============================="

* 1) Core key checks
isid gvkey fyear
assert !missing(gvkey, fyear)

* 2) Internal consistency of profit/loss decomposition
capture assert abs((profit - loss) - taxinc) < 1e-8 if !missing(profit, loss, taxinc)
if _rc {
    display as error "WARNING: profit - loss != taxinc for some observations."
}
capture assert profit >= 0 if !missing(profit)
capture assert loss   >= 0 if !missing(loss)

* 3) Refund variable should only be populated in policy years
capture assert missing(potential_refund) if !inlist(fyear, 2001, 2002, 2008, 2009)
if _rc {
    display as error "WARNING: potential_refund is populated outside policy years."
}

display as text "=============================="
display as text "SANITY CHECKS: COVERAGE"
display as text "=============================="
count if regflag_0203 == 1
display as result "Obs in regression sample (2002 policy outcomes: 2002/2003): " r(N)

count if regflag_2010 == 1
display as result "Obs in regression sample (2009 policy outcomes: 2010/2011): " r(N)

tab fyear if regflag_0203 == 1 | regflag_2010 == 1, missing

* 4) FF48 checks (if present)
capture confirm variable ffi48
if _rc == 0 {
    count if !missing(sic)
    display as result "Obs with nonmissing sic: " r(N)

    count if !missing(ffi48)
    display as result "Obs with nonmissing ffi48: " r(N)

    tab ffi48 if regflag_0203 == 1 | regflag_2010 == 1, missing
}
else {
    display as error "WARNING: ffi48 not found in dataset."
}

display as text "=============================="
display as text "TABLE: POLICY-PERIOD MEANS/SD"
display as text "=============================="

* Periods for table comparison: use actual regression samples
gen byte reg_period = .
replace reg_period = 2002 if regflag_0203 == 1
replace reg_period = 2009 if regflag_2010 == 1
label define policy_period 2002 "2002 policy period (FY2001-2002)" ///
                           2009 "2009 policy period (FY2010-2011 outcomes)"
label values reg_period policy_period

* Paper variables in display units:
* ($M) variables are already Compustat USD millions.
* Change in employment is in thousands in Compustat EMP.
xtset gvkey fyear

* Match regression timing for assignment variable and tax refund
gen v_assign_m = .
replace v_assign_m = L.assignment_v if regflag_0203 == 1
bysort gvkey: egen v_2009_val = max(cond(fyear == 2009, assignment_v, .))
replace v_assign_m = v_2009_val if regflag_2010 == 1

gen tax_refund_m = .
replace tax_refund_m = L.potential_refund if regflag_0203 == 1
bysort gvkey: egen refund_2010_val = max(cond(fyear == 2009, potential_refund, .))
replace tax_refund_m = refund_2010_val if regflag_2010 == 1

gen investment_m          = investment
gen change_cash_m         = d_cash
gen payout_m              = payout
gen change_debt_m         = d_totdebt
gen change_stinv_m        = d_stinv
gen change_ltinv_m        = d_inv
gen change_emp_thou       = d_emp

label var v_assign_m            "V (Assignment Variable, $M)"
label var tax_refund_m          "Tax Refund ($M)"
label var investment_m          "Investment ($M)"
label var change_cash_m         "Change in Cash ($M)"
label var payout_m              "Payout ($M)"
label var change_debt_m         "Change in Debt ($M)"
label var change_stinv_m        "Change in Short-term Investments ($M)"
label var change_ltinv_m        "Change in Long-term Investments ($M)"
label var change_emp_thou       "Change in Employment (Thousands)"
label var zscore                "Altman's z-score"
label var oscore                "Ohlson's o-score"

local table_vars ///
    v_assign_m ///
    tax_refund_m ///
    investment_m ///
    change_cash_m ///
    payout_m ///
    change_debt_m ///
    change_stinv_m ///
    change_ltinv_m ///
    change_emp_thou ///
    zscore ///
    oscore

* Mean / SD / N by regression sample period
foreach p in 2002 2009 {
    display as text "----------------------------------------"
    display as text "Regression sample period: `p'"
    foreach v of local table_vars {
        quietly summarize `v' if reg_period == `p' & has_crsp_link_obs == 1, detail
        display as text "`: variable label `v''"
        display as result "  mean = " %12.4f r(mean) "   sd = " %12.4f r(sd) "   N = " %9.0f r(N)
    }
}

display as text "========================================"
display as text "Done."
display as text "========================================"

