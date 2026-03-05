clear all
set more off

use "../data/main_data.dta", clear

capture confirm variable ffi48
if _rc {
    di as err "Missing industry variable ffi48. Run clean_data.do first."
    exit 111
}
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

* Assignment variable V in levels ($M): paper definition from clean_data.do
gen v_level = assignment_v

* ---------------------------------------------------------------------
* First stage: 2002 policy period (outcome years 2002 and 2003)
* TaxRefund in year t is mapped from potential_refund in t-1.
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

reg refund_0203 ///
    v1_0203 v2_0203 ///
    zv1_0203 zv2_0203 ///
    pre_0203_tobinq pre_0203_roa pre_0203_cf_assets pre_0203_sales_assets ///
    pre_0203_leverage pre_0203_ln_assets ///
    pre_0203_loss pre_0203_loss2 ///
    i.ffi48 ///
    if regflag_0203 == 1, vce(cluster ffi48)

test zv1_0203 zv2_0203
estimates store fs_0203

break
* ---------------------------------------------------------------------
* First stage: 2009 policy period
* Run separate cross-sections for 2010 and 2011 outcomes.
* ---------------------------------------------------------------------
bysort gvkey: egen refund_2010_val = max(cond(fyear == 2009, potential_refund, .))
bysort gvkey: egen v_2009_val      = max(cond(fyear == 2009, v_level, .))
bysort gvkey: egen loss_2009_val   = max(cond(fyear == 2009, loss, .))
foreach c in tobinq roa cf_assets sales_assets leverage ln_assets {
    bysort gvkey: egen `c'_2009_val = max(cond(fyear == 2009, `c', .))
}

foreach yy in 2010 2011 {
    gen in_`yy' = (fyear == `yy') & (regflag_2010 == 1)

    gen refund_`yy'_reg = refund_2010_val if in_`yy'
    gen v_`yy'      = v_2009_val      if in_`yy'
    gen d_`yy'      = (v_`yy' < 0)    if in_`yy' & !missing(v_`yy')
    gen v1_`yy'     = v_`yy'
    gen v2_`yy'     = v_`yy'^2
    gen zv1_`yy'    = d_`yy' * v1_`yy'
    gen zv2_`yy'    = d_`yy' * v2_`yy'

    gen pre_`yy'_tobinq       = tobinq_2009_val       if in_`yy'
    gen pre_`yy'_roa          = roa_2009_val          if in_`yy'
    gen pre_`yy'_cf_assets    = cf_assets_2009_val    if in_`yy'
    gen pre_`yy'_sales_assets = sales_assets_2009_val if in_`yy'
    gen pre_`yy'_leverage     = leverage_2009_val     if in_`yy'
    gen pre_`yy'_ln_assets    = ln_assets_2009_val    if in_`yy'
    gen pre_`yy'_loss         = loss_2009_val         if in_`yy'
    gen pre_`yy'_loss2        = pre_`yy'_loss^2

    reg refund_`yy'_reg ///
        v1_`yy' v2_`yy' ///
        zv1_`yy' zv2_`yy' ///
        pre_`yy'_tobinq pre_`yy'_roa pre_`yy'_cf_assets pre_`yy'_sales_assets ///
        pre_`yy'_leverage pre_`yy'_ln_assets ///
        pre_`yy'_loss pre_`yy'_loss2 ///
        i.ffi48 ///
        if in_`yy', vce(cluster ffi48)

    test zv1_`yy' zv2_`yy'
    estimates store fs_`yy'
}

display as text "Stored estimates: fs_0203 fs_2010 fs_2011"
