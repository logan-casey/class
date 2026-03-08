
use "../data/compustat_output1.dta", clear

********************************************************************************
* Sample Restrictions
********************************************************************************

destring(gvkey), replace

* Keep firms with total assets greater than $1 million
keep if at > 1 & !missing(at)

* Exclude financials/utilities/international/non-operating establishments by SIC
* Financials: 6000-6999, Utilities: 4900-4999, Intl/non-operating: 9000-9999
gen double sic_num = real(trim(sic))
drop if !missing(sic_num) & ///
    (inrange(sic_num, 4900, 4999) | inrange(sic_num, 6000, 6999) | inrange(sic_num, 9000, 9999))
drop sic_num

* Sort for panel operations
gen strL datafmt_l = lower(trim(datafmt))
gen byte datafmt_pref = (datafmt_l != "indl")
gen strL consol_l = lower(trim(consol))
gen byte consol_pref = (consol_l != "c")
gen double at_keep = at
replace at_keep = 0 if missing(at_keep)

* Prioritize consolidated industrial data, then by largest assets.
sort gvkey fyear datafmt_pref consol_pref -at_keep
by gvkey fyear: keep if _n == 1

drop datafmt_l datafmt_pref consol_l consol_pref at_keep

* Generate indicators for policy-relevant five-year carryback windows
* 2002-period assignment windows: 1996-2000 (for 2001 loss) or 1997-2001 (for 2002 loss)
* 2010-period assignment windows: 2003-2007 (for 2008 loss) or 2004-2008 (for 2009 loss)
gen win_2001_5 = inrange(fyear, 1996, 2000)
gen win_2002_5 = inrange(fyear, 1997, 2001)
gen win_2008_5 = inrange(fyear, 2003, 2007)
gen win_2009_5 = inrange(fyear, 2004, 2008)

* Count years observed per five-year window
bysort gvkey: egen n_win_2001_5 = total(win_2001_5)
bysort gvkey: egen n_win_2002_5 = total(win_2002_5)
bysort gvkey: egen n_win_2008_5 = total(win_2008_5)
bysort gvkey: egen n_win_2009_5 = total(win_2009_5)

* Keep firms present in at least one required five-year window
keep if (n_win_2001_5 == 5) | (n_win_2002_5 == 5) | ///
        (n_win_2008_5 == 5) | (n_win_2009_5 == 5)

* Regression eligibility flags based on required five-year carryback presence
gen firm_complete_0203 = (n_win_2001_5 == 5 | n_win_2002_5 == 5)
gen firm_complete_2010 = (n_win_2008_5 == 5 | n_win_2009_5 == 5)

* Row-level flags are finalized after all period-specific restrictions (Step 7)
gen byte regflag_0203 = .
gen byte regflag_2010 = .

drop win_2001_5 win_2002_5 win_2008_5 win_2009_5
drop n_win_2001_5 n_win_2002_5 n_win_2008_5 n_win_2009_5


********************************************************************************
* Taxable Income and Tax Rate
********************************************************************************

* --- Step 1: Construct Pretax Income when PIDOM is missing ---
* PI = OIADP - XINT + SPI + NOPI
gen pi_calc = oiadp - xint + spi + nopi

* Fill missing pi with pi_calc
replace pi = pi_calc if missing(pi)

* --- Step 2: Construct taxable income (as in dobridge_na.R) ---
local tau = 0.35

* Replace missing XIDO with zero
replace xido = 0 if missing(xido)

* Fallback PIDOM definition: PIDOM, else PI - PIFO
gen pidom_imp = pidom
replace pidom_imp = pi - cond(missing(pifo), 0, pifo) if missing(pidom_imp) & !missing(pi)
* MTR-style txfed fallback for TI: use TXDFED if available, else TXDI
gen txdfed_for_ti = txdfed
replace txdfed_for_ti = txdi if missing(txdfed_for_ti)
replace txdfed_for_ti = 0 if missing(txdfed_for_ti)

* Taxable income used in refunds
gen taxinc = pidom_imp - txdfed_for_ti / `tau' + xido / (1 - `tau')

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
gen taxrate = .
replace taxrate = txfed_use / taxinc if !missing(txfed_use) & txfed_use > 0 & !missing(taxinc) & taxinc > 0
replace taxrate = . if taxrate < 0.01 | taxrate > 0.52

********************************************************************************
* Fallback Marginal Tax Rate (Graham and Mills, 2008)
********************************************************************************

* Reuse pidom_imp for MTR fallback, when PIDOM is missing.

* MTR ratio: TXFED/PIDOM, with fallback TXT/PI when needed.
gen double tax_ratio_gm = .
replace tax_ratio_gm = txfed_use / pidom_imp if !missing(txfed_use) & !missing(pidom_imp) & pidom_imp != 0
replace tax_ratio_gm = txt / pi if missing(tax_ratio_gm) & !missing(txt) & !missing(pi) & pi != 0

* Graham-Mills fallback.
gen double mtr = 0.331 ///
    - 0.075 * (tax_ratio_gm < 0.1) ///
    - 0.012 * (tlcf > 0) ///
    - 0.106 * (pi < 0) ///
    + 0.037 * (abs(pifo / pi) > 0.05)

replace mtr = . if missing(tax_ratio_gm) & missing(txt) & missing(txfed_use)

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
label var mtr         "Marginal tax rate (Graham-Mills fallback)"



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
gen v2010_profit_window = .
gen v2010_loss_2008     = .
gen v2010_loss_2009     = .
gen v2010_chosen_loss   = .

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
    st_view(vpwin,   ., "v2010_profit_window")
    st_view(vl08,    ., "v2010_loss_2008")
    st_view(vl09,    ., "v2010_loss_2009")
    st_view(vch,     ., "v2010_chosen_loss")

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
        vl08[i] = loss_2008
        vl09[i] = loss_2009

        // Initialize option-specific objects each loop to avoid stale values
        profits_A = J(5,1,0)
        rates_A   = J(5,1,.)
        profits_B = J(5,1,0)
        rates_B   = J(5,1,.)

        if (loss_2008 > 0) {
            profits_A = adjust_2yr_profits(profits_all[1::5], loss_all[1::5])
            rates_A   = rates_all[1::5]
            refund_A  = calc_refund(loss_2008, profits_A, rates_A, 1)
        }
        else {
            refund_A = .
        }

        if (loss_2009 > 0) {
            profits_B = adjust_2yr_profits(profits_all[2::6], loss_all[2::6])
            rates_B   = rates_all[2::6]
            refund_B  = calc_refund(loss_2009, profits_B, rates_B, 1)
        }
        else {
            refund_B = .
        }

        // No eligible policy loss year: leave outputs missing
        if (loss_2008 <= 0 & loss_2009 <= 0) continue

        // Choose among eligible options only
        if (loss_2008 > 0 & (loss_2009 <= 0 | refund_A >= refund_B)) {
            // 2010 refund amount is based on the chosen 5yr-loss year (2008 here).
            loss_for_2010 = loss_2008
            r2010[i]  = refund_A
            lappl[i]  = loss_for_2010
            u08[i]    = 1
            vpwin[i]  = sum(profits_A)
            vch[i]    = loss_for_2010
            // Assignment V for 2010 refund year when 2008 option is chosen:
            // available profits in 2003-2007 minus chosen policy loss (2008)
            v2010[i]  = sum(profits_A) - loss_for_2010
        }
        else {
            loss_for_2010 = loss_2009
            r2010[i]  = refund_B
            lappl[i]  = loss_for_2010
            u08[i]    = 0
            vpwin[i]  = sum(profits_B)
            vch[i]    = loss_for_2010
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
bysort gvkey (fyear): replace v2010_profit_window = v2010_profit_window[_n-1] ///
    if fyear == 2009 & missing(v2010_profit_window)
bysort gvkey (fyear): replace v2010_loss_2008 = v2010_loss_2008[_n-1] ///
    if fyear == 2009 & missing(v2010_loss_2008)
bysort gvkey (fyear): replace v2010_loss_2009 = v2010_loss_2009[_n-1] ///
    if fyear == 2009 & missing(v2010_loss_2009)
bysort gvkey (fyear): replace v2010_chosen_loss = v2010_chosen_loss[_n-1] ///
    if fyear == 2009 & missing(v2010_chosen_loss)

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
label var v2010_profit_window "2010-policy V component: sum of adjusted carryback-window profits"
label var v2010_loss_2008     "2010-policy V component: 2008 loss"
label var v2010_loss_2009     "2010-policy V component: 2009 loss"
label var v2010_chosen_loss   "2010-policy V component: chosen policy loss (2008 or 2009)"
label var potential_refund "Potential refund from carryback policy"
label var refund_assets    "Potential refund / lagged assets"
label var assignment_v     "Assignment variable V (available carryback profits - policy losses)"

********************************************************************************
* STEP 7: Period-Specific Regression Sample Flags
********************************************************************************

* V = assignment variable in levels ($M)
gen v = assignment_v

* Policy-period windows (for restriction calculations)
gen byte policy_0203 = inlist(fyear, 2001, 2002)
gen byte policy_2010 = inlist(fyear, 2008, 2009)

* Loss eligibility by policy period
gen byte loss_elig_0203_obs = ///
    ((!missing(loss_2002) & loss_2002 > 0) | ///
     (!missing(loss_2003) & loss_2003 > 0))
bysort gvkey: egen byte firm_has_policy_loss_0203 = max(loss_elig_0203_obs)

gen byte loss_elig_2010_obs = (!missing(loss_applied_2010) & loss_applied_2010 > 0)
bysort gvkey: egen byte firm_has_policy_loss_2010 = max(loss_elig_2010_obs)

* ---------- Outlier rules for first-stage samples (row-level only) ----------
* 2002/2003 outcomes use lagged policy-year assignment/refund.
gen v_fs_0203 = L.assignment_v if inlist(fyear, 2002, 2003)
gen refund_fs_0203 = L.potential_refund if inlist(fyear, 2002, 2003)

quietly _pctile v_fs_0203 if !missing(v_fs_0203), p(1 99)
scalar v_p1_0203  = r(r1)
scalar v_p99_0203 = r(r2)

quietly _pctile refund_fs_0203 if !missing(refund_fs_0203), p(99.9)
scalar refund_p999_0203 = r(r1)

gen byte outlier_fs_0203 = inlist(fyear, 2002, 2003) & !missing(v_fs_0203, refund_fs_0203) & ///
    ((v_fs_0203 < v_p1_0203 | v_fs_0203 > v_p99_0203) | ///
     (refund_fs_0203 > refund_p999_0203))

* 2010/2011 outcomes use 2009 policy-year assignment/refund carried by firm.
bysort gvkey: egen v_2009_for_2010 = max(cond(fyear == 2009, assignment_v, .))
bysort gvkey: egen refund_2009_for_2010 = max(cond(fyear == 2009, potential_refund, .))
gen v_fs_2010 = v_2009_for_2010 if inlist(fyear, 2010, 2011)
gen refund_fs_2010 = refund_2009_for_2010 if inlist(fyear, 2010, 2011)

quietly _pctile v_fs_2010 if !missing(v_fs_2010), p(1 99)
scalar v_p1_2010  = r(r1)
scalar v_p99_2010 = r(r2)

quietly _pctile refund_fs_2010 if !missing(refund_fs_2010), p(99.9)
scalar refund_p999_2010 = r(r1)

gen byte outlier_fs_2010 = inlist(fyear, 2010, 2011) & !missing(v_fs_2010, refund_fs_2010) & ///
    ((v_fs_2010 < v_p1_2010 | v_fs_2010 > v_p99_2010) | ///
     (refund_fs_2010 > refund_p999_2010))

* Align with R: trim investment by 1st and 99th percentiles by period sample
quietly _pctile investment if inlist(fyear, 2002, 2003), p(1 99)
scalar invest_p1_0203  = r(r1)
scalar invest_p99_0203 = r(r2)
quietly _pctile investment if inlist(fyear, 2010, 2011), p(1 99)
scalar invest_p1_2010  = r(r1)
scalar invest_p99_2010 = r(r2)

gen byte invest_trim_0203 = inlist(fyear, 2002, 2003) & ///
    !missing(investment) & (investment >= invest_p1_0203 & investment <= invest_p99_0203)
gen byte invest_trim_2010 = inlist(fyear, 2010, 2011) & ///
    !missing(investment) & (investment >= invest_p1_2010 & investment <= invest_p99_2010)

* Final firm-level sample flags by policy period
gen byte sample_firm_0203 = firm_complete_0203 & firm_has_policy_loss_0203
gen byte sample_firm_2010 = firm_complete_2010 & firm_has_policy_loss_2010

* Row-level regression flags (use these in first/second-stage scripts)
replace regflag_0203 = inlist(fyear, 2002, 2003) & sample_firm_0203 & ///
    !missing(v_fs_0203, refund_fs_0203) & (outlier_fs_0203 == 0) & ///
    (invest_trim_0203 == 1)
replace regflag_2010 = inlist(fyear, 2010, 2011) & sample_firm_2010 & ///
    !missing(v_fs_2010, refund_fs_2010) & (outlier_fs_2010 == 0) & ///
    (invest_trim_2010 == 1)

* Cleanup temporary variables and scalars
drop policy_0203 policy_2010
drop loss_elig_0203_obs loss_elig_2010_obs
drop v_fs_0203 refund_fs_0203 v_2009_for_2010 refund_2009_for_2010 v_fs_2010 refund_fs_2010
drop outlier_fs_0203 outlier_fs_2010
drop invest_trim_0203 invest_trim_2010
scalar drop v_p1_0203 v_p99_0203
scalar drop refund_p999_0203
scalar drop v_p1_2010 v_p99_2010
scalar drop refund_p999_2010
scalar drop invest_p1_0203 invest_p99_0203
scalar drop invest_p1_2010 invest_p99_2010

********************************************************************************
* STEP 8: Flag firms/observations with a valid CRSP link (CCM)
********************************************************************************

* Observation date used for link validity: datadate if available, else year-end.
gen double obsdate = datadate
replace obsdate = mdy(12, 31, fyear) if missing(obsdate) & !missing(fyear)
format obsdate %td

* Stable row id for merge-back after many-to-many join.
gen long obs_id = _n

preserve
keep obs_id gvkey obsdate
tempfile base_obs
save `base_obs'
restore

preserve
use "../data/ccm_link.dta", clear
rename GVKEY gvkey
rename LPERMNO permno
rename LINKDT linkdt
rename LINKENDDT linkend
capture confirm variable LINKTYPE
if (_rc == 0) rename LINKTYPE linktype
capture confirm variable LINKPRIM
if (_rc == 0) rename LINKPRIM linkprim

capture confirm variable linktype
if (_rc == 0) {
    tostring linktype, replace
    replace linktype = trim(upper(linktype))
}

capture confirm variable linkprim
if (_rc == 0) {
    tostring linkprim, replace
    replace linkprim = trim(upper(linkprim))
}

capture confirm numeric variable linkdt
if (_rc) {
    gen double linkdt_num = date(linkdt, "YMD")
    replace linkdt_num = date(linkdt, "MDY") if missing(linkdt_num) & !missing(linkdt)
}
else {
    gen double linkdt_num = linkdt
}

capture confirm numeric variable linkend
if (_rc) {
    gen double linkend_num = date(linkend, "YMD")
    replace linkend_num = date(linkend, "MDY") if missing(linkend_num) & !missing(linkend)
}
else {
    gen double linkend_num = linkend
}

replace linkend_num = td(31dec9999) if missing(linkend_num)
format linkdt_num linkend_num %td
destring gvkey, replace

* Keep CCM link types used in R: LC/LU/LS and P/C
capture confirm variable linktype
if (_rc == 0) keep if inlist(linktype, "LC", "LU", "LS")
capture confirm variable linkprim
if (_rc == 0) keep if inlist(linkprim, "P", "C")

keep gvkey permno linkdt_num linkend_num
drop if missing(gvkey) | missing(permno)
tempfile ccm
save `ccm'
restore

preserve
use `base_obs', clear
joinby gvkey using `ccm'
keep if inrange(obsdate, linkdt_num, linkend_num)
bysort obs_id: keep if _n == 1
keep obs_id
gen byte has_crsp_link_obs = 1
tempfile crsp_hits
save `crsp_hits'
restore

merge 1:1 obs_id using `crsp_hits', nogen
replace has_crsp_link_obs = 0 if missing(has_crsp_link_obs)
bysort gvkey: egen byte has_crsp_link_firm = max(has_crsp_link_obs)

label var has_crsp_link_obs  "Obs has valid CCM gvkey-permno link at obsdate"
label var has_crsp_link_firm "Firm has any valid CCM gvkey-permno link in sample"

* Optional: enforce CRSP-link requirement inside regression flags.
* Default off to align with paper-based sample that does not require CRSP coverage.
local require_crsp_link_for_reg 0
if `require_crsp_link_for_reg' {
    replace regflag_0203 = 0 if inlist(fyear, 2002, 2003) & has_crsp_link_obs == 0
    replace regflag_2010 = 0 if fyear == 2010 & has_crsp_link_obs == 0
}

drop obs_id obsdate

save "../data/main_data.dta", replace
