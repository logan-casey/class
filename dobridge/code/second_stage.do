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

xtset gvkey fyear

* Assignment variable V in levels ($M): profits minus losses
gen v_level = assignment_v

* Outcomes from specification (first policy only)
* "other" = change in short-term investments + change in long-term investments + acquisitions
gen other = d_stinv + d_inv + acquisitions

* Winsorize outcomes at 1/99 within the 2002/2003 regression sample
* before constructing total uses.
foreach y in investment d_cash d_totdebt payout other d_emp {
    quietly _pctile `y' if regflag_0203 == 1 & inlist(fyear, 2002, 2003) & !missing(`y'), p(1 99)
    local p1_`y' = r(r1)
    local p99_`y' = r(r2)
    replace `y' = `p1_`y'' if regflag_0203 == 1 & inlist(fyear, 2002, 2003) & `y' < `p1_`y''
    replace `y' = `p99_`y'' if regflag_0203 == 1 & inlist(fyear, 2002, 2003) & `y' > `p99_`y''
}

* "total" adds uses and nets out debt issuance (so debt enters with opposite sign)
gen total = investment + d_cash - d_totdebt + payout + other
local outcomes investment d_cash d_totdebt payout other total d_emp

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

* Control timing switch:
* 0 = contemporaneous controls (year of refund receipt)
* 1 = lagged controls (year before refund receipt)
local use_lagged_controls 0

foreach c in tobinq roa cf_assets sales_assets leverage ln_assets mtr {
    gen pre_0203_`c' = `c' if inlist(fyear, 2002, 2003)
    if `use_lagged_controls' {
        replace pre_0203_`c' = L.`c' if inlist(fyear, 2002, 2003)
    }
}
gen pre_0203_loss  = loss if inlist(fyear, 2002, 2003)
if `use_lagged_controls' {
    replace pre_0203_loss = L.loss if inlist(fyear, 2002, 2003)
}
gen pre_0203_loss2 = pre_0203_loss^2

foreach y of local outcomes {
     ivregress 2sls `y' ///
         v1_0203 v2_0203 ///
         pre_0203_tobinq pre_0203_roa pre_0203_cf_assets pre_0203_sales_assets ///
         pre_0203_leverage pre_0203_ln_assets pre_0203_mtr ///
         pre_0203_loss pre_0203_loss2 ///
         i.ffi48 ///
         (refund_0203 = zv1_0203 zv2_0203) ///
        if regflag_0203 == 1, vce(cluster ffi48)
    estimates store ss_0203_`y'

    scalar b_`y'  = _b[refund_0203]
    scalar se_`y' = _se[refund_0203]
    scalar p_`y'  = 2*normal(-abs(_b[refund_0203]/_se[refund_0203]))
    scalar N_`y'  = e(N)
    scalar r2_`y' = e(r2)
}

* ---------------------------------------------------------------------
* Log summary (Table-5 style)
* ---------------------------------------------------------------------
display as text " "
display as text "Table 5-style summary (2002 policy period, refund years 2002/2003)"
display as text "Outcome                      Coef (Tax Refund)     SE          N        R2"
foreach y in investment d_cash d_totdebt payout other total d_emp {
    local b_`y'  : display %9.4f scalar(b_`y')
    local se_`y' : display %9.4f scalar(se_`y')
    local n_`y'  : display %10.0fc scalar(N_`y')
    local r_`y'  : display %7.3f scalar(r2_`y')
    display as text %24s "`y'" "   " as result "`b_`y''" "   " as text "[" as result "`se_`y''" as text "]   " as result "`n_`y''" "   " "`r_`y''"
}

* ---------------------------------------------------------------------
* Export LaTeX table
* ---------------------------------------------------------------------
cap mkdir "../output"

foreach y in investment d_cash d_totdebt payout other total d_emp {
    local star_`y' ""
    if (scalar(p_`y') < 0.10) local star_`y' "*"
    if (scalar(p_`y') < 0.05) local star_`y' "**"
    if (scalar(p_`y') < 0.01) local star_`y' "***"

    local bfmt_`y'  : display %9.4f scalar(b_`y')
    local sefmt_`y' : display %9.4f scalar(se_`y')
    local nfmt_`y'  : display %12.0fc scalar(N_`y')
    local rfmt_`y'  : display %6.3f scalar(r2_`y')
}

file open fh using "../output/table5_second_stage_2002.tex", write replace
file write fh "\begin{table}[htbp]" _n
file write fh "\centering" _n
file write fh "\caption{Tax Refund Allocation in the 2002 Policy Period (Years of Receipt: 2002 and 2003)}" _n
file write fh "\begin{tabular}{lccccccc}" _n
file write fh "\hline\hline" _n
file write fh "Dependent variable & Investment & Change in Cash & Change in Total Debt & Payout & Other Uses & Total & Change in Employment \\" _n
file write fh " & (1) & (2) & (3) & (4) & (5) & (6) & (7) \\" _n
file write fh "\hline" _n
file write fh "Tax Refund & `bfmt_investment'`star_investment' & `bfmt_d_cash'`star_d_cash' & `bfmt_d_totdebt'`star_d_totdebt' & `bfmt_payout'`star_payout' & `bfmt_other'`star_other' & `bfmt_total'`star_total' & `bfmt_d_emp'`star_d_emp' \\" _n
file write fh " & [`sefmt_investment'] & [`sefmt_d_cash'] & [`sefmt_d_totdebt'] & [`sefmt_payout'] & [`sefmt_other'] & [`sefmt_total'] & [`sefmt_d_emp'] \\" _n
file write fh "Controls & + & + & + & + & + & + & + \\" _n
file write fh "Industry F.E. & + & + & + & + & + & + & + \\" _n
file write fh "Observations & `nfmt_investment' & `nfmt_d_cash' & `nfmt_d_totdebt' & `nfmt_payout' & `nfmt_other' & `nfmt_total' & `nfmt_d_emp' \\" _n
file write fh "R-squared & `rfmt_investment' & `rfmt_d_cash' & `rfmt_d_totdebt' & `rfmt_payout' & `rfmt_other' & `rfmt_total' & `rfmt_d_emp' \\" _n
file write fh "\hline\hline" _n
file write fh "\multicolumn{8}{l}{\footnotesize Standard errors in brackets; clustered by FF48 industry.} \\" _n
file write fh "\multicolumn{8}{l}{\footnotesize *** p<0.01, ** p<0.05, * p<0.10.} \\" _n
file write fh "\end{tabular}" _n
file write fh "\end{table}" _n
file close fh

display as text "Stored estimates created for first-policy (2002/2003) regressions."
display as text "LaTeX table written: ../output/table5_second_stage_2002.tex"


