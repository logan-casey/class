clear all
set more off

use "../data/main_data.dta", clear

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
count if inlist(fyear, 2001, 2002, 2008, 2009)
display as result "Obs in policy years (2001, 2002, 2008, 2009): " r(N)

count if inlist(fyear, 2001, 2002, 2008, 2009) & !missing(potential_refund)
display as result "Obs with nonmissing potential_refund in policy years: " r(N)

tab fyear if inlist(fyear, 2001, 2002, 2008, 2009), missing

* 4) FF48 checks (if present)
capture confirm variable ffi48
if _rc == 0 {
    count if !missing(sic)
    display as result "Obs with nonmissing sic: " r(N)

    count if !missing(ffi48)
    display as result "Obs with nonmissing ffi48: " r(N)

    tab ffi48 if inlist(fyear, 2001, 2002, 2008, 2009), missing
}
else {
    display as error "WARNING: ffi48 not found in dataset."
}

display as text "=============================="
display as text "TABLE: POLICY-PERIOD MEANS/SD"
display as text "=============================="

* Periods for table comparison
gen byte policy_period = .
replace policy_period = 2002 if inlist(fyear, 2001, 2002)
replace policy_period = 2009 if inlist(fyear, 2008, 2009)
label define policy_period 2002 "2002 policy period (FY2001-2002)" ///
                           2009 "2009 policy period (FY2008-2009)"
label values policy_period policy_period

* Paper variables in display units:
* ($M) variables are already Compustat USD millions.
* Change in employment is in thousands in Compustat EMP.
gen v_profit_minus_loss_m = profit - loss
gen tax_refund_m          = potential_refund
gen investment_m          = investment
gen change_cash_m         = d_cash
gen payout_m              = payout
gen change_debt_m         = d_totdebt
gen change_stinv_m        = d_stinv
gen change_ltinv_m        = d_inv
gen change_emp_thou       = d_emp

label var v_profit_minus_loss_m "V (Profits - Losses) ($M)"
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
    v_profit_minus_loss_m ///
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

* Mean / SD / N by policy period
foreach p in 2002 2009 {
    display as text "----------------------------------------"
    display as text "Policy period: `p'"
    foreach v of local table_vars {
        quietly summarize `v' if policy_period == `p', detail
        display as text "`: variable label `v''"
        display as result "  mean = " %12.4f r(mean) "   sd = " %12.4f r(sd) "   N = " %9.0f r(N)
    }
}

display as text "========================================"
display as text "Done."
display as text "========================================"
