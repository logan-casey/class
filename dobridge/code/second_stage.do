clear all
set more off

use "../data/main_data.dta", clear

capture confirm variable ffi48
if _rc {
    di as err "Missing industry variable ffi48. Run clean_data.do first."
    exit 111
}

xtset gvkey fyear

* Assignment variable V in levels ($M): profits minus losses
gen v_level = assignment_v

* Outcomes from specification (first policy only)
* "other" = change in short-term investments + change in long-term investments + acquisitions
gen other = d_stinv + d_inv + acquisitions
local outcomes investment payout d_cash d_totdebt other d_emp

* ---------------------------------------------------------------------
* Second stage: 2002 policy period (outcome years 2002 and 2003)
* ---------------------------------------------------------------------
gen refund_0203 = L.potential_refund if inlist(fyear, 2002, 2003)
gen v_0203      = L.v_level          if inlist(fyear, 2002, 2003)
gen d_0203      = (v_0203 < 0)       if !missing(v_0203) & inlist(fyear, 2002, 2003)
gen v1_0203     = v_0203
gen v2_0203     = v_0203^2
gen zv1_0203    = d_0203 * v1_0203
gen zv2_0203    = d_0203 * v2_0203

foreach c in tobinq roa cf_assets sales_assets leverage ln_assets {
    gen pre_0203_`c' = L.`c' if inlist(fyear, 2002, 2003)
}
gen pre_0203_loss  = L.loss    if inlist(fyear, 2002, 2003)
gen pre_0203_loss2 = pre_0203_loss^2

foreach y of local outcomes {
    ivregress 2sls `y' ///
        v1_0203 v2_0203 ///
        pre_0203_tobinq pre_0203_roa pre_0203_cf_assets pre_0203_sales_assets ///
        pre_0203_leverage pre_0203_ln_assets ///
        pre_0203_loss pre_0203_loss2 ///
        i.ffi48 ///
        (refund_0203 = zv1_0203 zv2_0203) ///
        if regflag_0203, vce(cluster ffi48)
    estimates store ss_0203_`y'
}

display as text "Stored estimates created for first-policy (2002/2003) regressions."


