clear all
set more off

use "../data/main_data.dta", clear

capture confirm variable ffi48
if _rc {
    di as err "Missing ffi48. Run clean_data.do first."
    exit 111
}
capture confirm variable regflag_0203
if _rc {
    di as err "Missing regflag_0203. Run clean_data.do first."
    exit 111
}

xtset gvkey fyear
cap mkdir "../output"

* ---------------------------------------------------------------------
* Core second-stage timing for 2002 policy period (refund years 2002/2003)
* ---------------------------------------------------------------------
gen v_level = assignment_v
gen refund_0203 = L.potential_refund if inlist(fyear, 2002, 2003)
gen v_0203      = L.v_level          if inlist(fyear, 2002, 2003)
gen d_0203      = (v_0203 < 0)       if inlist(fyear, 2002, 2003) & !missing(v_0203)
gen v1_0203     = v_0203
gen v2_0203     = v_0203^2
gen zv1_0203    = d_0203 * v1_0203
gen zv2_0203    = d_0203 * v2_0203

* Pre-treatment controls (year prior to refund receipt).
foreach c in tobinq roa cf_assets sales_assets leverage ln_assets mtr {
    gen pre_0203_`c' = L.`c' if inlist(fyear, 2002, 2003)
}
gen pre_0203_loss  = L.loss if inlist(fyear, 2002, 2003)
gen pre_0203_loss2 = pre_0203_loss^2

* ---------------------------------------------------------------------
* Subgroup variables
* ---------------------------------------------------------------------
* Prior-year payout and DFF-style payout-to-operating-income ratio.
gen pre_payout = L.dvc + L.prstkc if inlist(fyear, 2002, 2003)
gen pre_dff = pre_payout / abs(L.oibdp) if inlist(fyear, 2002, 2003) & !missing(pre_payout, L.oibdp) & L.oibdp != 0

* Prior-year multinational indicator: |PIFO/PI| > 0.05
gen pre_multi = (abs(L.pifo / L.pi) > 0.05) if inlist(fyear, 2002, 2003) & !missing(L.pifo, L.pi) & L.pi != 0

* Prior-year KZ index. Use existing variable if present; otherwise construct.
capture confirm variable kz_index
if _rc == 0 {
    gen pre_kz = L.kz_index if inlist(fyear, 2002, 2003)
}
else {
    capture confirm variable kz
    if _rc == 0 {
        gen pre_kz = L.kz if inlist(fyear, 2002, 2003)
    }
    else {
        * Construct KZ only if required inputs exist; otherwise leave missing.
        capture confirm variable ppent
        local have_ppent = (_rc == 0)
        capture confirm variable dvp
        local have_dvp = (_rc == 0)

        if `have_ppent' & `have_dvp' {
            gen ppent_lag = L.ppent
            gen prccf_use = prcc_f
            capture confirm variable prcc_c
            if _rc == 0 replace prccf_use = prcc_c if missing(prccf_use)

            capture confirm variable txdb
            if _rc == 0 {
                gen txdb_use = txdb
                replace txdb_use = 0 if missing(txdb_use)
            }
            else {
                gen txdb_use = 0
            }

            gen kz_calc = .
            replace kz_calc = ///
                -1.002 * ((ib + dp) / ppent_lag) + ///
                 0.283 * ((at + prccf_use*csho - ceq - txdb_use) / at) + ///
                 3.139 * ((dltt + dlc) / (dltt + seq)) - ///
                39.368 * ((dvc + dvp) / ppent_lag) - ///
                 1.315 * (che / ppent_lag) ///
                if !missing(ib, dp, ppent_lag, at, prccf_use, csho, ceq, dltt, dlc, seq, dvc, dvp, che) & ///
                   ppent_lag > 0 & at != 0 & (dltt + seq) != 0
            gen pre_kz = L.kz_calc if inlist(fyear, 2002, 2003)
        }
        else {
            gen pre_kz = . if inlist(fyear, 2002, 2003)
        }
    }
}

* Base sample for Table 6 (same structural sample as second stage, investment outcome).
gen byte base_ok = regflag_0203 == 1 & inlist(fyear, 2002, 2003) & ///
    !missing(investment, refund_0203, v1_0203, v2_0203, zv1_0203, zv2_0203, ///
             pre_0203_tobinq, pre_0203_roa, pre_0203_cf_assets, pre_0203_sales_assets, ///
             pre_0203_leverage, pre_0203_ln_assets, pre_0203_mtr, pre_0203_loss, pre_0203_loss2, ffi48)

quietly count if base_ok
di as text "Table 6 base sample (investment spec): " %12.0fc r(N)

* Availability diagnostics.
quietly count if base_ok & !missing(pre_kz)
local n_kz = r(N)
quietly count if base_ok & !missing(pre_payout)
local n_payout = r(N)
quietly count if base_ok & !missing(pre_dff)
local n_dff = r(N)
quietly count if base_ok & !missing(pre_multi)
local n_multi = r(N)
quietly count if base_ok & !missing(pre_0203_tobinq)
local n_q = r(N)

local have_kz = (`n_kz' > 0)
local have_q = (`n_q' > 0)

* Quartiles for KZ and Tobin's Q in base sample (if available).
if `have_kz' {
    quietly _pctile pre_kz if base_ok & !missing(pre_kz), p(25 75)
    scalar kz_p25 = r(r1)
    scalar kz_p75 = r(r2)
}
else {
    scalar kz_p25 = .
    scalar kz_p75 = .
}

if `have_q' {
    quietly _pctile pre_0203_tobinq if base_ok & !missing(pre_0203_tobinq), p(25 75)
    scalar q_p25 = r(r1)
    scalar q_p75 = r(r2)
}
else {
    scalar q_p25 = .
    scalar q_p75 = .
}

di as text "Availability in base sample:"
di as text "  pre_kz:      " %12.0fc `n_kz'
di as text "  pre_payout:  " %12.0fc `n_payout'
di as text "  pre_dff:     " %12.0fc `n_dff'
di as text "  pre_multi:   " %12.0fc `n_multi'
di as text "  pre_tobinq:  " %12.0fc `n_q'

tempfile t6res
postfile P str1 panel int col str40 label double b se p N r2 using `t6res', replace

program define _run_t6, rclass
    syntax , Panel(string) Col(integer) Label(string) Cond(string)
    capture noisily ivregress 2sls investment ///
        v1_0203 v2_0203 ///
        pre_0203_tobinq pre_0203_roa pre_0203_cf_assets pre_0203_sales_assets ///
        pre_0203_leverage pre_0203_ln_assets pre_0203_mtr ///
        pre_0203_loss pre_0203_loss2 ///
        i.ffi48 ///
        (refund_0203 = zv1_0203 zv2_0203) ///
        if base_ok & (`cond'), vce(cluster ffi48)

    if _rc == 0 {
        return scalar b  = _b[refund_0203]
        return scalar se = _se[refund_0203]
        return scalar p  = 2*normal(-abs(_b[refund_0203]/_se[refund_0203]))
        return scalar N  = e(N)
        return scalar r2 = e(r2)
    }
    else {
        return scalar b  = .
        return scalar se = .
        return scalar p  = .
        return scalar N  = .
        return scalar r2 = .
    }
end

* Panel A: Financial constraints
if `have_kz' {
    quietly _run_t6, panel("A") col(1) label("KZ: constrained (top quartile)") cond("pre_kz >= kz_p75 & !missing(pre_kz)")
    post P ("A") (1) ("KZ: constrained (top quartile)") (r(b)) (r(se)) (r(p)) (r(N)) (r(r2))

    quietly _run_t6, panel("A") col(2) label("KZ: unconstrained (bottom quartile)") cond("pre_kz <= kz_p25 & !missing(pre_kz)")
    post P ("A") (2) ("KZ: unconstrained (bottom quartile)") (r(b)) (r(se)) (r(p)) (r(N)) (r(r2))
}
else {
    post P ("A") (1) ("KZ: constrained (top quartile)") (.) (.) (.) (.) (.)
    post P ("A") (2) ("KZ: unconstrained (bottom quartile)") (.) (.) (.) (.) (.)
}

quietly _run_t6, panel("A") col(3) label("Payout=0 (constrained)") cond("pre_payout == 0 & !missing(pre_payout)")
post P ("A") (3) ("Payout=0 (constrained)") (r(b)) (r(se)) (r(p)) (r(N)) (r(r2))

quietly _run_t6, panel("A") col(4) label("Payout>0 (unconstrained)") cond("pre_payout > 0 & !missing(pre_payout)")
post P ("A") (4) ("Payout>0 (unconstrained)") (r(b)) (r(se)) (r(p)) (r(N)) (r(r2))

quietly _run_t6, panel("A") col(5) label("DFF<=0 (constrained)") cond("pre_dff <= 0 & !missing(pre_dff)")
post P ("A") (5) ("DFF<=0 (constrained)") (r(b)) (r(se)) (r(p)) (r(N)) (r(r2))

quietly _run_t6, panel("A") col(6) label("DFF>0 (unconstrained)") cond("pre_dff > 0 & !missing(pre_dff)")
post P ("A") (6) ("DFF>0 (unconstrained)") (r(b)) (r(se)) (r(p)) (r(N)) (r(r2))

* Panel B: Investment opportunities
if `have_q' {
    quietly _run_t6, panel("B") col(1) label("High Tobin's Q (top quartile)") cond("pre_0203_tobinq >= q_p75 & !missing(pre_0203_tobinq)")
    post P ("B") (1) ("High Tobin's Q (top quartile)") (r(b)) (r(se)) (r(p)) (r(N)) (r(r2))

    quietly _run_t6, panel("B") col(2) label("Low Tobin's Q (bottom quartile)") cond("pre_0203_tobinq <= q_p25 & !missing(pre_0203_tobinq)")
    post P ("B") (2) ("Low Tobin's Q (bottom quartile)") (r(b)) (r(se)) (r(p)) (r(N)) (r(r2))
}
else {
    post P ("B") (1) ("High Tobin's Q (top quartile)") (.) (.) (.) (.) (.)
    post P ("B") (2) ("Low Tobin's Q (bottom quartile)") (.) (.) (.) (.) (.)
}

quietly _run_t6, panel("B") col(3) label("Refund year 2002") cond("fyear == 2002")
post P ("B") (3) ("Refund year 2002") (r(b)) (r(se)) (r(p)) (r(N)) (r(r2))

quietly _run_t6, panel("B") col(4) label("Refund year 2003") cond("fyear == 2003")
post P ("B") (4) ("Refund year 2003") (r(b)) (r(se)) (r(p)) (r(N)) (r(r2))

quietly _run_t6, panel("B") col(5) label("Multinational") cond("pre_multi == 1 & !missing(pre_multi)")
post P ("B") (5) ("Multinational") (r(b)) (r(se)) (r(p)) (r(N)) (r(r2))

quietly _run_t6, panel("B") col(6) label("Domestic") cond("pre_multi == 0 & !missing(pre_multi)")
post P ("B") (6) ("Domestic") (r(b)) (r(se)) (r(p)) (r(N)) (r(r2))

postclose P

preserve
use `t6res', clear
sort panel col
export delimited using "../output/table6_subgroups_results.csv", replace
di as text "CSV written: ../output/table6_subgroups_results.csv"
di as text " "
di as text "Table 6 regression summary (Tax Refund coefficient):"
list panel col label b se p N r2, noobs sepby(panel)

forvalues j = 1/6 {
    quietly summarize b if panel == "A" & col == `j', meanonly
    local A_b`j' : display %9.3f r(mean)
    quietly summarize se if panel == "A" & col == `j', meanonly
    local A_se`j' : display %9.3f r(mean)
    quietly summarize p if panel == "A" & col == `j', meanonly
    local A_p`j' = r(mean)
    local A_star`j' ""
    if (`A_p`j'' < 0.10) local A_star`j' "*"
    if (`A_p`j'' < 0.05) local A_star`j' "**"
    if (`A_p`j'' < 0.01) local A_star`j' "***"
    quietly summarize N if panel == "A" & col == `j', meanonly
    local A_N`j' : display %12.0fc r(mean)
    quietly summarize r2 if panel == "A" & col == `j', meanonly
    local A_r2`j' : display %6.3f r(mean)
}

forvalues j = 1/6 {
    quietly summarize b if panel == "B" & col == `j', meanonly
    local B_b`j' : display %9.3f r(mean)
    quietly summarize se if panel == "B" & col == `j', meanonly
    local B_se`j' : display %9.3f r(mean)
    quietly summarize p if panel == "B" & col == `j', meanonly
    local B_p`j' = r(mean)
    local B_star`j' ""
    if (`B_p`j'' < 0.10) local B_star`j' "*"
    if (`B_p`j'' < 0.05) local B_star`j' "**"
    if (`B_p`j'' < 0.01) local B_star`j' "***"
    quietly summarize N if panel == "B" & col == `j', meanonly
    local B_N`j' : display %12.0fc r(mean)
    quietly summarize r2 if panel == "B" & col == `j', meanonly
    local B_r2`j' : display %6.3f r(mean)
}

file open fh using "../output/table6_subgroups.tex", write replace
file write fh "\begin{table}[htbp]" _n
file write fh "\centering" _n
file write fh "\caption{Investment, Financial Constraints, and Investment Opportunities: 2002 Policy Period}" _n
file write fh "\begin{tabular}{lcccccc}" _n
file write fh "\hline\hline" _n
file write fh "\multicolumn{7}{l}{\textbf{Panel A: Financial Constraints}} \\" _n
file write fh " & \multicolumn{2}{c}{KZ} & \multicolumn{2}{c}{Payout} & \multicolumn{2}{c}{DFF} \\" _n
file write fh "Financially constrained? & Yes & No & Yes & No & Yes & No \\" _n
file write fh " & (1) & (2) & (3) & (4) & (5) & (6) \\" _n
file write fh "Tax Refund & `A_b1'`A_star1' & `A_b2'`A_star2' & `A_b3'`A_star3' & `A_b4'`A_star4' & `A_b5'`A_star5' & `A_b6'`A_star6' \\" _n
file write fh " & [`A_se1'] & [`A_se2'] & [`A_se3'] & [`A_se4'] & [`A_se5'] & [`A_se6'] \\" _n
file write fh "Controls & + & + & + & + & + & + \\" _n
file write fh "Industry F.E. & + & + & + & + & + & + \\" _n
file write fh "Observations & `A_N1' & `A_N2' & `A_N3' & `A_N4' & `A_N5' & `A_N6' \\" _n
file write fh "R-squared & `A_r21' & `A_r22' & `A_r23' & `A_r24' & `A_r25' & `A_r26' \\" _n
file write fh "\hline" _n
file write fh "\multicolumn{7}{l}{\textbf{Panel B: Investment Opportunities}} \\" _n
file write fh " & \multicolumn{2}{c}{Tobin's Q} & \multicolumn{2}{c}{Year of Tax Refund} & \multicolumn{2}{c}{Multinational} \\" _n
file write fh "High investment opportunities? & Yes & No & No & Yes & Yes & No \\" _n
file write fh " & (1) & (2) & (3) & (4) & (5) & (6) \\" _n
file write fh "Tax Refund & `B_b1'`B_star1' & `B_b2'`B_star2' & `B_b3'`B_star3' & `B_b4'`B_star4' & `B_b5'`B_star5' & `B_b6'`B_star6' \\" _n
file write fh " & [`B_se1'] & [`B_se2'] & [`B_se3'] & [`B_se4'] & [`B_se5'] & [`B_se6'] \\" _n
file write fh "Controls & + & + & + & + & + & + \\" _n
file write fh "Industry F.E. & + & + & + & + & + & + \\" _n
file write fh "Observations & `B_N1' & `B_N2' & `B_N3' & `B_N4' & `B_N5' & `B_N6' \\" _n
file write fh "R-squared & `B_r21' & `B_r22' & `B_r23' & `B_r24' & `B_r25' & `B_r26' \\" _n
file write fh "\hline\hline" _n
file write fh "\multicolumn{7}{l}{\footnotesize Standard errors in brackets; clustered by FF48 industry.} \\" _n
file write fh "\multicolumn{7}{l}{\footnotesize *** p<0.01, ** p<0.05, * p<0.10.} \\" _n
file write fh "\end{tabular}" _n
file write fh "\end{table}" _n
file close fh
di as text "LaTeX written: ../output/table6_subgroups.tex"

if !`have_kz' {
    di as err "Cannot estimate KZ quartile columns with current data: pre_kz unavailable."
}
if !`have_q' {
    di as err "Cannot estimate Tobin's Q quartile columns with current data: pre_0203_tobinq unavailable."
}
if (`n_multi' == 0) {
    di as err "Cannot estimate multinational/domestic columns with current data: pre_multi unavailable."
}
if (`n_dff' == 0) {
    di as err "Cannot estimate DFF columns with current data: pre_dff unavailable."
}
restore
