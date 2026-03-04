clear all
set more off

use "../data/main_data.dta", clear

capture confirm variable regflag_0203
if _rc {
    di as err "Missing regflag_0203. Re-run clean_data.do."
    exit 111
}

xtset gvkey fyear
cap mkdir "../output"

* ---------------------------------------------------------------------
* Build 2002-period regression sample variables (same timing as first_stage.do)
* ---------------------------------------------------------------------
gen v = L.assignment_v if regflag_0203 == 1
gen tax_refund = L.potential_refund if regflag_0203 == 1
gen inv = investment if regflag_0203 == 1

foreach c in tobinq roa cf_assets sales_assets leverage ln_assets {
    gen pre_`c' = L.`c' if regflag_0203 == 1
}
gen pre_loss  = L.loss if regflag_0203 == 1
gen pre_loss2 = pre_loss^2

* Residualize outcomes on controls + industry FE (recentered residuals)
reg tax_refund pre_tobinq pre_roa pre_cf_assets pre_sales_assets ///
    pre_leverage pre_ln_assets pre_loss pre_loss2 i.ffi48 ///
    if regflag_0203 == 1
predict tax_refund_resid if e(sample), resid

reg inv pre_tobinq pre_roa pre_cf_assets pre_sales_assets ///
    pre_leverage pre_ln_assets pre_loss pre_loss2 i.ffi48 ///
    if regflag_0203 == 1
predict inv_resid if e(sample), resid

* ---------------------------------------------------------------------
* 1M bins of V with 99% CI for means
* ---------------------------------------------------------------------
preserve
keep if regflag_0203 == 1 & !missing(v, tax_refund_resid, inv_resid)

gen bin = floor(v)
gen bin_mid = bin + 0.5

collapse ///
    (mean) tax_mean=tax_refund_resid inv_mean=inv_resid ///
    (sd)   tax_sd=tax_refund_resid inv_sd=inv_resid ///
    (count) n_tax=tax_refund_resid n_inv=inv_resid, by(bin bin_mid)

gen tax_se = tax_sd / sqrt(n_tax)
gen inv_se = inv_sd / sqrt(n_inv)
local z99 = invnormal(0.995)
gen tax_lo99 = tax_mean - `z99' * tax_se
gen tax_hi99 = tax_mean + `z99' * tax_se
gen inv_lo99 = inv_mean - `z99' * inv_se
gen inv_hi99 = inv_mean + `z99' * inv_se

tempfile binned
save `binned'
restore

* ---------------------------------------------------------------------
* Left panel: Tax refund kink
* ---------------------------------------------------------------------
twoway ///
    (lfit tax_refund_resid v if regflag_0203 == 1 & v < 0, lcolor(blue) lwidth(medthick)) ///
    (lfit tax_refund_resid v if regflag_0203 == 1 & v >= 0, lcolor(blue) lwidth(medthick)) ///
    (rarea tax_lo99 tax_hi99 bin_mid using `binned', color(gs12%45) lcolor(none)) ///
    (scatter tax_mean bin_mid using `binned', mcolor(blue) msymbol(Oh) msize(small)), ///
    xline(0, lcolor(red) lpattern(dash)) ///
    xtitle("Assignment variable V ($M)") ///
    ytitle("Tax Refund Residual ($M)") ///
    title("Tax Refund (2002 Policy Period)") ///
    legend(off) ///
    name(fig2_left, replace)

* ---------------------------------------------------------------------
* Right panel: Investment kink
* ---------------------------------------------------------------------
twoway ///
    (lfit inv_resid v if regflag_0203 == 1 & v < 0, lcolor(navy) lwidth(medthick)) ///
    (lfit inv_resid v if regflag_0203 == 1 & v >= 0, lcolor(navy) lwidth(medthick)) ///
    (rarea inv_lo99 inv_hi99 bin_mid using `binned', color(gs12%45) lcolor(none)) ///
    (scatter inv_mean bin_mid using `binned', mcolor(navy) msymbol(Oh) msize(small)), ///
    xline(0, lcolor(red) lpattern(dash)) ///
    xtitle("Assignment variable V ($M)") ///
    ytitle("Investment Residual ($M)") ///
    title("Investment (2002 Policy Period)") ///
    legend(off) ///
    name(fig2_right, replace)

graph combine fig2_left fig2_right, col(2) ///
    title("Figure 2: Examples of Kinks in Variables")

graph export "../output/figure2_kinks_2002.png", width(2400) replace
