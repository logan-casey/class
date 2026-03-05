
use "../data/compustat_output1.dta", clear

********************************************************************************
* Sample Restrictions
********************************************************************************

destring(gvkey), replace

* Keep firms with total assets greater than $1 million
keep if at > 1 & !missing(at)

* Sort for panel operations
gsort gvkey fyear -datadate
duplicates drop gvkey fyear, force

* Generate indicators for policy-specific year coverage
* 2002 refund period coverage: 1996-2001
* 2003 refund period coverage: 1997-2002
* 2010 refund period coverage: 2003-2008
gen pre2002 = (fyear >= 1996 & fyear <= 2001)
gen pre2003 = (fyear >= 1997 & fyear <= 2002)
gen pre2009 = (fyear >= 2003 & fyear <= 2008)

* Count number of years each firm appears in each window
bysort gvkey: egen years_pre2002 = total(pre2002)
bysort gvkey: egen years_pre2003 = total(pre2003)
bysort gvkey: egen years_pre2009 = total(pre2009)

* Keep firms present in at least one policy-relevant window
keep if (years_pre2002 == 6) | (years_pre2003 == 6) | (years_pre2009 == 6)

* Regression eligibility flags based on complete preceding-year panels
* For 2002/2003 refund regressions: require full 1996-2002 history
gen pre_0203_full = inrange(fyear, 1996, 2002)
bysort gvkey: egen n_pre_0203_full = total(pre_0203_full)
gen firm_complete_0203 = (n_pre_0203_full == 7)

* For 2010/2011 refund regressions: require full 2003-2009 history
gen pre_2010_full = inrange(fyear, 2003, 2009)
bysort gvkey: egen n_pre_2010_full = total(pre_2010_full)
gen firm_complete_2010 = (n_pre_2010_full == 7)

* Row-level flags are finalized after all period-specific restrictions (Step 7)
gen byte regflag_0203 = .
gen byte regflag_2010 = .

drop pre2002 pre2003 pre2009 years_pre2002 years_pre2003 years_pre2009
drop pre_0203_full pre_2010_full n_pre_0203_full n_pre_2010_full


********************************************************************************
* Taxable Income and Tax Rate
********************************************************************************

* --- Step 1: Construct Pretax Income when PIDOM is missing ---
* PI = OIADP - XINT + SPI + NOPI
gen pi_calc = oiadp - xint + spi + nopi

* Fill missing pi with pi_calc
replace pi = pi_calc if missing(pi)

* --- Step 2: Construct taxable income ---
local tau = 0.35

* Replace missing XIDO with zero
replace xido = 0 if missing(xido)

* Primary: use PIDOM and TXDFED
gen taxinc = pidom - txdfed / `tau' + xido / (1 - `tau') ///
    if !missing(pidom) & !missing(txdfed)

* Secondary: use PI and TXDI when PIDOM or TXDFED are missing
replace taxinc = pi - txdi / `tau' + xido / (1 - `tau') ///
    if missing(taxinc) & !missing(pi) & !missing(txdi)

* If TXDI is also missing, use PI with no deferred tax adjustment
replace taxinc = pi + xido / (1 - `tau') ///
    if missing(taxinc) & !missing(pi) & missing(txdi)

* --- Step 3: Construct Federal Income Taxes ---
* Replace missing components with zero before computing fallback
foreach v in txfo txs txdi txo {
    replace `v' = 0 if missing(`v')
}

* Fallback federal taxes = TXT - TXFO - TXS - TXDI - TXO
gen txfed_calc = txt - txfo - txs - txdi - txo

* Use TXFED if available, otherwise use calculated value
gen txfed_use = txfed
replace txfed_use = txfed_calc if missing(txfed_use)

* --- Step 4: Tax Rate ---
gen taxrate = txfed_use / taxinc
* Winsorize or bound as needed in later steps; raw rate constructed here

********************************************************************************
* Lagged Variables (needed for changes)
********************************************************************************

* Declare panel
xtset gvkey fyear

* Lags
gen ch_lag      = L.ch
gen dlc_lag     = L.dlc
gen dltt_lag    = L.dltt
gen ivst_lag    = L.ivst
gen emp_lag     = L.emp
gen ni_lag1     = L.ni       // NI(t-1)
gen ni_lag2     = L2.ni      // NI(t-2)

********************************************************************************
* Outcome Variables
********************************************************************************

* Investment
gen investment = capxv - sppe

* Change in cash
gen d_cash = ch - ch_lag

* Change in total debt
gen d_totdebt = (dlc + dltt) - (dlc_lag + dltt_lag)

* Change in long-term debt
gen d_ltdebt = dltt - dltt_lag

* Change in short-term debt
gen d_stdebt = dlc - dlc_lag

* Payout
gen payout = dvc + prstkc

* Change in short-term investments
gen d_stinv = ivst - ivst_lag

* Change in investments
gen d_inv = ivch

* Acquisitions
gen acquisitions = aqc

* Change in employment
gen d_emp = emp - emp_lag

********************************************************************************
* Firm Characteristics
********************************************************************************
ffi48_macro sic, gen(ffi48) desc(ffi48_desc) replace

* Drop financial firms: FF48 = Banks (44), Insur (45), RlEst (46), Fin (47)
// drop if inrange(ffi48, 44, 47)

* Altman's Z-score
* EBIT approximated as OIADP (operating income after depreciation) or use ebit if available
* The paper uses EBIT; Compustat does not have a direct EBIT field so use OIADP
gen zscore = 3.3*(ebit/at) + 1.0*(sale/at) + 1.4*(re/at) + 1.2*(wcap/at)

* Ohlson's O-score
* Indicator: net income negative in both t-1 and t-2
gen ni_neg2 = (ni_lag1 < 0 & ni_lag2 < 0 & !missing(ni_lag1) & !missing(ni_lag2))

* Indicator: liabilities exceed assets
gen insolvent = (lt > at & !missing(lt) & !missing(at))

* Change in NI scaled by sum of absolute values
gen d_ni_scaled = (ni - ni_lag1) / (abs(ni) + abs(ni_lag1)) ///
    if !missing(ni) & !missing(ni_lag1) & (abs(ni) + abs(ni_lag1)) != 0
replace d_ni_scaled = 0 if missing(d_ni_scaled)

gen oscore = -1.32 ///
    - 0.407  * ln(at) ///
    + 6.03   * (lt/at) ///
    - 1.43   * (act - lct)/at ///
    + 0.0757 * (lct/act) ///
    - 2.37   * (ni/at) ///
    - 1.83   * ((pi + dp)/lt) ///
    + 0.285  * ni_neg2 ///
    - 1.72   * insolvent ///
    - 0.521  * d_ni_scaled

* ROA
gen roa = oibdp / at

* Tobin's Q
* SEQ + TXDITC - PSTK = book equity; handle missing TXDITC and PSTK
replace txditc = 0 if missing(txditc)
replace pstk   = 0 if missing(pstk)
gen tobinq = (at + prcc_f*csho - (seq + txditc - pstk)) / at

* Cash Flow / Assets
gen cf_assets = (ib + dp) / at

* Ln(Assets)
gen ln_assets = ln(at)

* Leverage
gen leverage = (dlc + dltt) / (dlc + dltt + prcc_f*csho)

* Sales / Assets
gen sales_assets = sale / at

********************************************************************************
* Label Variables
********************************************************************************

label var taxinc      "Taxable Income"
label var taxrate     "Tax Rate (Federal Taxes / Taxable Income)"
label var investment  "Investment (CAPXV - SPPE)"
label var d_cash      "Change in Cash"
label var d_totdebt   "Change in Total Debt"
label var d_ltdebt    "Change in Long-Term Debt"
label var d_stdebt    "Change in Short-Term Debt"
label var payout      "Payout (DVC + PRSTKC)"
label var d_stinv     "Change in Short-Term Investments"
label var d_inv       "Change in Investments (IVCH)"
label var acquisitions "Acquisitions (AQC)"
label var d_emp       "Change in Employment"
label var zscore      "Altman's Z-Score"
label var oscore      "Ohlson's O-Score"
label var roa         "ROA (OIBDP/AT)"
label var tobinq      "Tobin's Q"
label var cf_assets   "Cash Flow / Assets"
label var ln_assets   "Ln(Assets)"
label var leverage    "Leverage"
label var sales_assets "Sales / Assets"



********************************************************************************
* STEP 1: Define profit/loss variables
********************************************************************************

* Profits = positive taxable income, Losses = negative taxable income (as positive number)
gen profit = max(taxinc, 0)
gen loss   = max(-taxinc, 0)

* We'll need taxrate by year too — already constructed as taxrate

********************************************************************************
* STEP 2: Two-year carryback adjustments
* Implemented as helper functions and applied within each policy window
********************************************************************************

sort gvkey fyear

mata:
real colvector adjust_2yr_profits(real colvector profits, real colvector losses) {
    real colvector padj
    real scalar n, t, loss_rem
    n = rows(profits)
    padj = profits

    for (t = 3; t <= n; t++) {
        if (losses[t] <= 0) continue
        loss_rem = losses[t]

        // Apply to t-2 first
        if (padj[t-2] > 0) {
            if (padj[t-2] >= loss_rem) {
                padj[t-2] = padj[t-2] - loss_rem
                loss_rem = 0
            }
            else {
                loss_rem = loss_rem - padj[t-2]
                padj[t-2] = 0
            }
        }

        // Then apply remaining to t-1
        if (loss_rem > 0 & padj[t-1] > 0) {
            if (padj[t-1] >= loss_rem) {
                padj[t-1] = padj[t-1] - loss_rem
                loss_rem = 0
            }
            else {
                loss_rem = loss_rem - padj[t-1]
                padj[t-1] = 0
            }
        }
    }
    return(padj)
}

real colvector apply_loss_oldest(real scalar loss_amt, real colvector profits) {
    real colvector prem
    real scalar j, loss_rem, p
    prem = profits
    loss_rem = loss_amt

    for (j = 1; j <= rows(prem); j++) {
        if (loss_rem <= 0) break
        p = prem[j]
        if (p <= 0) continue
        if (loss_rem >= p) {
            prem[j] = 0
            loss_rem = loss_rem - p
        }
        else {
            prem[j] = p - loss_rem
            loss_rem = 0
        }
    }
    return(prem)
}

real scalar calc_refund(real scalar loss_amt,
                        real colvector profits,
                        real colvector rates,
                        real scalar half_fifth) {
    real scalar refund, loss_rem, j, p, r
    refund = 0
    loss_rem = loss_amt

    for (j = 1; j <= rows(profits); j++) {
        p = profits[j]
        if (j == 1 & half_fifth == 1) p = p * 0.5
        r = rates[j]

        if (loss_rem <= 0) break
        if (p <= 0) continue

        if (loss_rem >= p) {
            refund   = refund + r * p
            loss_rem = loss_rem - p
        }
        else {
            refund   = refund + r * loss_rem
            loss_rem = 0
        }
    }
    return(refund)
}
end

********************************************************************************
* STEP 4: 2002 POLICY - Refunds based on 2001 and 2002 losses
* Window for 2001 loss: profits from 1996-2000 (after 2yr carryback adj)
* Window for 2002 loss: profits from 1997-2001 (after 2yr carryback adj)
********************************************************************************

gen refund_2002 = .
gen refund_2003 = .
gen loss_2002   = .
gen loss_2003   = .
gen assign_v_2002 = .
gen assign_v_2003 = .

mata:
void policy2002() {
    st_view(gv,     ., "gvkey")
    st_view(fy,     ., "fyear")
    st_view(lv,     ., "loss")
    st_view(pv,     ., "profit")
    st_view(tr,     ., "taxrate")
    st_view(r2002,  ., "refund_2002")
    st_view(r2003,  ., "refund_2003")
    st_view(l2002,  ., "loss_2002")
    st_view(l2003,  ., "loss_2003")
    st_view(v2002,  ., "assign_v_2002")
    st_view(v2003,  ., "assign_v_2003")

    n = rows(gv)

    for (i = 1; i <= n; i++) {
        g = gv[i]
        y = fy[i]

        // 2002 refund: based on fiscal year 2001 loss
        if (y == 2001) {
            loss_amt = lv[i]
            if (loss_amt <= 0) {
                r2002[i] = 0
                l2002[i] = 0
                continue
            }
            l2002[i] = loss_amt

            profits = J(5, 1, 0)
            losses  = J(5, 1, 0)
            rates   = J(5, 1, 0)
            years   = (1996, 1997, 1998, 1999, 2000)
            for (j = 1; j <= n; j++) {
                if (gv[j] != g) continue
                for (k = 1; k <= 5; k++) {
                    if (fy[j] == years[k]) {
                        profits[k] = pv[j]
                        losses[k]  = lv[j]
                        rates[k]   = tr[j]
                    }
                }
            }
            profits_adj = adjust_2yr_profits(profits, losses)
            // Assignment V for 2002 refund year:
            // available profits in 1996-2000 minus 2001 policy loss
            v2002[i] = sum(profits_adj) - loss_amt
            r2002[i] = calc_refund(loss_amt, profits_adj, rates, 0)
        }

        // 2003 refund: based on fiscal year 2002 loss
        if (y == 2002) {
            loss_amt = lv[i]
            if (loss_amt <= 0) {
                r2003[i] = 0
                l2003[i] = 0
                continue
            }
            l2003[i] = loss_amt

            profits_9601 = J(6, 1, 0)
            losses_9601  = J(6, 1, 0)
            rates_9601   = J(6, 1, 0)
            years = (1996, 1997, 1998, 1999, 2000, 2001)
            for (j = 1; j <= n; j++) {
                if (gv[j] != g) continue
                for (k = 1; k <= 6; k++) {
                    if (fy[j] == years[k]) {
                        profits_9601[k] = pv[j]
                        losses_9601[k]  = lv[j]
                        rates_9601[k]   = tr[j]
                    }
                }
            }

            // Exclude 2001 policy loss from the generic 2yr adjustment pass;
            // it is consumed separately via the 2002-policy oldest-year rule.
            losses_for_adj_9601 = losses_9601
            losses_for_adj_9601[6] = 0
            profits_adj_9601 = adjust_2yr_profits(profits_9601, losses_for_adj_9601)
            loss_2001 = losses_9601[6]
            profits_remaining = apply_loss_oldest(loss_2001, profits_adj_9601)

            profits_2003 = profits_remaining[2::6]
            rates_2003   = rates_9601[2::6]
            // Assignment V for 2003 refund year:
            // available profits in 1997-2001 minus 2002 policy loss
            v2003[i] = sum(profits_2003) - loss_amt
            r2003[i] = calc_refund(loss_amt, profits_2003, rates_2003, 0)
        }
    }
}
end

mata: policy2002()

********************************************************************************
* STEP 5: 2009 POLICY - Refunds based on 2008 and/or 2009 losses
* Window for 2008 loss: profits from 2003-2007 (after 2yr carryback adj)
* Window for 2009 loss: profits from 2004-2008 (after 2yr carryback adj)
* 50% cap on 5th year (oldest year) profits
* Firms choose whichever loss generates higher refund
********************************************************************************

gen refund_2010       = .
gen loss_applied_2010 = .
gen used_2008_loss    = .
gen assign_v_2010     = .

mata:
void policy2009() {
    st_view(gv,      ., "gvkey")
    st_view(fy,      ., "fyear")
    st_view(lv,      ., "loss")
    st_view(pv,      ., "profit")
    st_view(tr,      ., "taxrate")
    st_view(r2010,   ., "refund_2010")
    st_view(lappl,   ., "loss_applied_2010")
    st_view(u08,     ., "used_2008_loss")
    st_view(v2010,   ., "assign_v_2010")

    n = rows(gv)

    for (i = 1; i <= n; i++) {
        g = gv[i]
        y = fy[i]

        if (y != 2008 & y != 2009) continue

        years_all = (2003,2004,2005,2006,2007,2008,2009)
        profits_all = J(7,1,0)
        rates_all   = J(7,1,0)
        loss_all    = J(7,1,0)

        for (j = 1; j <= n; j++) {
            if (gv[j] != g) continue
            for (k = 1; k <= 7; k++) {
                if (fy[j] == years_all[k]) {
                    profits_all[k] = pv[j]
                    rates_all[k]   = tr[j]
                    loss_all[k]    = lv[j]
                }
            }
        }

        // Only compute once per firm (at 2008 obs)
        if (y != 2008) continue

        loss_2008 = loss_all[6]
        loss_2009 = loss_all[7]

        if (loss_2008 > 0) {
            profits_A = adjust_2yr_profits(profits_all[1::5], loss_all[1::5])
            rates_A   = rates_all[1::5]
            refund_A  = calc_refund(loss_2008, profits_A, rates_A, 1)
        }
        else {
            refund_A = 0
        }

        if (loss_2009 > 0) {
            profits_B = adjust_2yr_profits(profits_all[2::6], loss_all[2::6])
            rates_B   = rates_all[2::6]
            refund_B  = calc_refund(loss_2009, profits_B, rates_B, 1)
        }
        else {
            refund_B = 0
        }

        if (refund_A >= refund_B) {
            // 2010 refund amount is based on the chosen 5yr-loss year (2008 here).
            loss_for_2010 = loss_2008
            r2010[i]  = refund_A
            lappl[i]  = loss_for_2010
            u08[i]    = 1
            // Assignment V for 2010 refund year when 2008 option is chosen:
            // available profits in 2003-2007 minus (2008 + 2009) policy losses
            v2010[i]  = sum(profits_A) - (loss_2008 + loss_2009)
        }
        else {
            loss_for_2010 = loss_2009
            r2010[i]  = refund_B
            lappl[i]  = loss_for_2010
            u08[i]    = 0
            // Assignment V for 2010 refund year when 2009 option is chosen:
            // available profits in 2004-2008 minus 2009 policy loss
            v2010[i]  = sum(profits_B) - loss_for_2010
        }
    }
}
end

mata: policy2009()

* Propagate 2010 refund info to 2009 obs for the same firm
bysort gvkey (fyear): replace refund_2010    = refund_2010[_n-1]    ///
    if fyear == 2009 & missing(refund_2010)
bysort gvkey (fyear): replace used_2008_loss = used_2008_loss[_n-1] ///
    if fyear == 2009 & missing(used_2008_loss)
bysort gvkey (fyear): replace assign_v_2010 = assign_v_2010[_n-1] ///
    if fyear == 2009 & missing(assign_v_2010)

********************************************************************************
* STEP 6: Consolidate refund variables
********************************************************************************

* Create a single refund variable relevant to each observation year
gen potential_refund = .
replace potential_refund = refund_2002 if fyear == 2001
replace potential_refund = refund_2003 if fyear == 2002
replace potential_refund = refund_2010 if inlist(fyear, 2008, 2009)

* Assignment variable V = carryback-window profits available - policy losses
* 2002 refund (fyear==2001): sum profits 1996-2000 - 2001 loss
* 2003 refund (fyear==2002): sum profits 1997-2001 - 2002 loss
* 2010 refund (fyear==2008/2009): either
*   sum profits 2003-2007 - (2008+2009 losses), or
*   sum profits 2004-2008 - 2009 loss,
* depending on the chosen policy-loss year (higher refund option)
gen assignment_v = .
replace assignment_v = assign_v_2002 if fyear == 2001
replace assignment_v = assign_v_2003 if fyear == 2002
replace assignment_v = assign_v_2010 if inlist(fyear, 2008, 2009)

* Scale by assets for use as a regressor (common in the paper)
gen refund_assets = potential_refund / L.at

label var refund_2002      "Estimated 2002 tax refund (2001 loss)"
label var refund_2003      "Estimated 2003 tax refund (2002 loss)"
label var refund_2010      "Estimated 2010 tax refund (2008 or 2009 loss)"
label var used_2008_loss   "=1 if firm applies 2008 loss to 5yr carryback"
label var potential_refund "Potential refund from carryback policy"
label var refund_assets    "Potential refund / lagged assets"
label var assignment_v     "Assignment variable V (available carryback profits - policy losses)"

********************************************************************************
* STEP 7: Period-Specific Regression Sample Flags
********************************************************************************

* V = assignment variable in levels ($M)
gen v = assignment_v

* Scale flows by lagged assets (fallback to current assets if lag is missing)
// gen flow_scale_at = L.at
// replace flow_scale_at = at if missing(flow_scale_at) | flow_scale_at <= 0
// gen d_cash_scaled = d_cash / flow_scale_at
// gen d_debt_scaled = d_totdebt / flow_scale_at
// gen d_inv_scaled  = d_inv / flow_scale_at
// gen payout_scaled = payout / flow_scale_at

* Policy-period windows (for restriction calculations)
gen byte policy_0203 = inlist(fyear, 2001, 2002)
gen byte policy_2010 = inlist(fyear, 2008, 2009)

* Loss eligibility by policy period
gen byte loss_elig_0203_obs = (loss_2002 > 0 | loss_2003 > 0)
bysort gvkey: egen byte firm_has_policy_loss_0203 = max(loss_elig_0203_obs)

gen byte loss_elig_2010_obs = (loss_applied_2010 > 0)
bysort gvkey: egen byte firm_has_policy_loss_2010 = max(loss_elig_2010_obs)

* ---------- Outlier rules for 2002/2003 sample ----------
quietly _pctile v if policy_0203 & !missing(v), p(1 99)
scalar v_p1_0203  = r(r1)
scalar v_p99_0203 = r(r2)

quietly _pctile investment if policy_0203 & !missing(investment), p(1 99)
scalar inv_p1_0203  = r(r1)
scalar inv_p99_0203 = r(r2)

// quietly _pctile d_cash_scaled if policy_0203 & !missing(d_cash_scaled), p(1 99)
// scalar cash_p1_0203  = r(r1)
// scalar cash_p99_0203 = r(r2)

// quietly _pctile d_debt_scaled if policy_0203 & !missing(d_debt_scaled), p(1 99)
// scalar debt_p1_0203  = r(r1)
// scalar debt_p99_0203 = r(r2)

// quietly _pctile d_inv_scaled if policy_0203 & !missing(d_inv_scaled), p(1 99)
// scalar dinv_p1_0203  = r(r1)
// scalar dinv_p99_0203 = r(r2)

// quietly _pctile payout_scaled if policy_0203 & !missing(payout_scaled), p(1 99)
// scalar payout_p1_0203  = r(r1)
// scalar payout_p99_0203 = r(r2)

quietly _pctile potential_refund if policy_0203 & !missing(potential_refund), p(99.5)
scalar refund_p995_0203 = r(r1)

gen byte drop_0203_obs = policy_0203 & ///
    ((v < v_p1_0203 | v > v_p99_0203) | ///
    (investment < inv_p1_0203 | investment > inv_p99_0203) | ///
    // (d_cash_scaled < cash_p1_0203 | d_cash_scaled > cash_p99_0203) | ///
    // (d_debt_scaled < debt_p1_0203 | d_debt_scaled > debt_p99_0203) | ///
    // (d_inv_scaled < dinv_p1_0203 | d_inv_scaled > dinv_p99_0203) | ///
    // (payout_scaled < payout_p1_0203 | payout_scaled > payout_p99_0203) | ///
    (potential_refund > refund_p995_0203))
bysort gvkey: egen byte drop_firm_0203 = max(drop_0203_obs)

* ---------- Outlier rules for 2010/2011 sample ----------
quietly _pctile v if policy_2010 & !missing(v), p(1 99)
scalar v_p1_2010  = r(r1)
scalar v_p99_2010 = r(r2)

quietly _pctile investment if policy_2010 & !missing(investment), p(1 99)
scalar inv_p1_2010  = r(r1)
scalar inv_p99_2010 = r(r2)

// quietly _pctile d_cash_scaled if policy_2010 & !missing(d_cash_scaled), p(1 99)
// scalar cash_p1_2010  = r(r1)
// scalar cash_p99_2010 = r(r2)

// quietly _pctile d_debt_scaled if policy_2010 & !missing(d_debt_scaled), p(1 99)
// scalar debt_p1_2010  = r(r1)
// scalar debt_p99_2010 = r(r2)

// quietly _pctile d_inv_scaled if policy_2010 & !missing(d_inv_scaled), p(1 99)
// scalar dinv_p1_2010  = r(r1)
// scalar dinv_p99_2010 = r(r2)

// quietly _pctile payout_scaled if policy_2010 & !missing(payout_scaled), p(1 99)
// scalar payout_p1_2010  = r(r1)
// scalar payout_p99_2010 = r(r2)

quietly _pctile potential_refund if policy_2010 & !missing(potential_refund), p(99.5)
scalar refund_p995_2010 = r(r1)

gen byte drop_2010_obs = policy_2010 & ///
    ((v < v_p1_2010 | v > v_p99_2010) | ///
    (investment < inv_p1_2010 | investment > inv_p99_2010) | ///
    // (d_cash_scaled < cash_p1_2010 | d_cash_scaled > cash_p99_2010) | ///
    // (d_debt_scaled < debt_p1_2010 | d_debt_scaled > debt_p99_2010) | ///
    // (d_inv_scaled < dinv_p1_2010 | d_inv_scaled > dinv_p99_2010) | ///
    // (payout_scaled < payout_p1_2010 | payout_scaled > payout_p99_2010) | ///
    (potential_refund > refund_p995_2010))
bysort gvkey: egen byte drop_firm_2010 = max(drop_2010_obs)

* Final firm-level sample flags by policy period
gen byte sample_firm_0203 = firm_complete_0203 & firm_has_policy_loss_0203 & ///
    (drop_firm_0203 == 0)
gen byte sample_firm_2010 = firm_complete_2010 & firm_has_policy_loss_2010 & ///
    (drop_firm_2010 == 0)

* Row-level regression flags (use these in first/second-stage scripts)
replace regflag_0203 = inlist(fyear, 2002, 2003) & sample_firm_0203
replace regflag_2010 = inlist(fyear, 2010, 2011) & sample_firm_2010

* Cleanup temporary variables and scalars
// drop flow_scale_at d_cash_scaled d_debt_scaled d_inv_scaled payout_scaled
drop policy_0203 policy_2010
drop loss_elig_0203_obs loss_elig_2010_obs
drop at_ok_0203_obs at_ok_2010_obs n_at_ok_0203 n_at_ok_2010
drop drop_0203_obs drop_2010_obs drop_firm_0203 drop_firm_2010
scalar drop v_p1_0203 v_p99_0203 inv_p1_0203 inv_p99_0203
// scalar drop debt_p1_0203 debt_p99_0203 dinv_p1_0203 dinv_p99_0203 payout_p1_0203 payout_p99_0203  cash_p1_0203 cash_p99_0203
scalar drop refund_p995_0203
scalar drop v_p1_2010 v_p99_2010 inv_p1_2010 inv_p99_2010
// scalar drop debt_p1_2010 debt_p99_2010 dinv_p1_2010 dinv_p99_2010 payout_p1_2010 payout_p99_2010 cash_p1_2010 cash_p99_2010
scalar drop refund_p995_2010

save "../data/main_data.dta", replace
