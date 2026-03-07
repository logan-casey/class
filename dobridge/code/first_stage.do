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

foreach c in tobinq roa cf_assets sales_assets leverage ln_assets mtr {
    gen pre_0203_`c' = `c' if inlist(fyear, 2002, 2003)
}
gen pre_0203_loss  = loss    if inlist(fyear, 2002, 2003)
gen pre_0203_loss2 = pre_0203_loss^2

reg refund_0203 ///
     v1_0203 v2_0203 ///
     zv1_0203 zv2_0203 ///
     pre_0203_tobinq pre_0203_roa pre_0203_cf_assets pre_0203_sales_assets ///
     pre_0203_leverage pre_0203_ln_assets pre_0203_mtr ///
     pre_0203_loss pre_0203_loss2 ///
     i.ffi48 ///
     if regflag_0203 == 1, vce(cluster ffi48)

scalar b_0203  = _b[zv1_0203]
scalar se_0203 = _se[zv1_0203]
scalar p_0203  = 2*ttail(e(df_r), abs(_b[zv1_0203]/_se[zv1_0203]))
scalar N_0203  = e(N)
scalar r2_0203 = e(r2)

test zv1_0203 zv2_0203
scalar pF_0203 = r(p)
estimates store fs_0203

* ---------------------------------------------------------------------
* First stage: 2009 policy period
* Run 2010 outcome only.
* ---------------------------------------------------------------------
bysort gvkey: egen refund_2010_val = max(cond(fyear == 2009, potential_refund, .))
bysort gvkey: egen v_2009_val      = max(cond(fyear == 2009, v_level, .))

local yy = 2010
gen in_`yy' = (fyear == `yy') & (regflag_2010 == 1)

gen refund_`yy'_reg = refund_2010_val if in_`yy'
gen v_`yy'      = v_2009_val      if in_`yy'
gen d_`yy'      = (v_`yy' < 0)    if in_`yy' & !missing(v_`yy')
gen v1_`yy'     = v_`yy'
gen v2_`yy'     = v_`yy'^2
gen zv1_`yy'    = d_`yy' * v1_`yy'
gen zv2_`yy'    = d_`yy' * v2_`yy'

gen pre_`yy'_tobinq       = tobinq       if in_`yy'
gen pre_`yy'_roa          = roa          if in_`yy'
gen pre_`yy'_cf_assets    = cf_assets    if in_`yy'
gen pre_`yy'_sales_assets = sales_assets if in_`yy'
gen pre_`yy'_leverage     = leverage     if in_`yy'
gen pre_`yy'_ln_assets    = ln_assets    if in_`yy'
gen pre_`yy'_mtr         = mtr          if in_`yy'
gen pre_`yy'_loss         = loss         if in_`yy'
gen pre_`yy'_loss2        = pre_`yy'_loss^2

reg refund_`yy'_reg ///
    v1_`yy' v2_`yy' ///
    zv1_`yy' zv2_`yy' ///
    pre_`yy'_tobinq pre_`yy'_roa pre_`yy'_cf_assets pre_`yy'_sales_assets ///
    pre_`yy'_leverage pre_`yy'_ln_assets pre_`yy'_mtr ///
    pre_`yy'_loss pre_`yy'_loss2 ///
    i.ffi48 ///
    if in_`yy', vce(cluster ffi48)

scalar b_2010  = _b[zv1_`yy']
scalar se_2010 = _se[zv1_`yy']
scalar p_2010  = 2*ttail(e(df_r), abs(_b[zv1_`yy']/_se[zv1_`yy']))
scalar N_2010  = e(N)
scalar r2_2010 = e(r2)

test zv1_`yy' zv2_`yy'
scalar pF_2010 = r(p)
estimates store fs_`yy'

* ---------------------------------------------------------------------
* Export LaTeX first-stage table (zv1 coefficient = change in slope at kink)
* ---------------------------------------------------------------------
cap mkdir "../output"

local star_0203 ""
if (scalar(p_0203) < 0.10) local star_0203 "*"
if (scalar(p_0203) < 0.05) local star_0203 "**"
if (scalar(p_0203) < 0.01) local star_0203 "***"

local star_2010 ""
if (scalar(p_2010) < 0.10) local star_2010 "*"
if (scalar(p_2010) < 0.05) local star_2010 "**"
if (scalar(p_2010) < 0.01) local star_2010 "***"

local b0203   : display %6.3f scalar(b_0203)
local se0203  : display %7.4f scalar(se_0203)
local pf0203  : display %4.2f scalar(pF_0203)
local n0203   : display %12.0fc scalar(N_0203)
local r20203  : display %5.3f scalar(r2_0203)

local b2010   : display %6.3f scalar(b_2010)
local se2010  : display %7.4f scalar(se_2010)
local pf2010  : display %4.2f scalar(pF_2010)
local n2010   : display %12.0fc scalar(N_2010)
local r22010  : display %5.3f scalar(r2_2010)

file open fh using "../output/first_stage_table.tex", write replace
file write fh "\begin{table}[htbp]" _n
file write fh "\centering" _n
file write fh "\caption{First-stage regression: tax refund kink coefficients}" _n
file write fh "\begin{tabular}{lcc}" _n
file write fh "\hline\hline" _n
file write fh "Dependent variable = Tax Refund & 2002 Policy & 2009 Policy \\" _n
file write fh " & (1) & (2) \\" _n
file write fh "\hline" _n
file write fh "Change in Slope & `b0203'`star_0203' & `b2010'`star_2010' \\" _n
file write fh " & [`se0203'] & [`se2010'] \\" _n
file write fh "Controls & + & + \\" _n
file write fh "Industry F.E. & + & + \\" _n
file write fh "F-test (p-value) & `pf0203' & `pf2010' \\" _n
file write fh "Chi-squared test for coefficient differences (p-value) &  &  \\" _n
file write fh "Observations & `n0203' & `n2010' \\" _n
file write fh "R-squared & `r20203' & `r22010' \\" _n
file write fh "\hline\hline" _n
file write fh "\end{tabular}" _n
file write fh "\end{table}" _n
file close fh

display as text "Stored estimates: fs_0203 fs_2010"
display as text "LaTeX table written: ../output/first_stage_table.tex"
