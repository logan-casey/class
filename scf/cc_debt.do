
clear all
set maxvar 6000
*---------------------------------------------------------------*
* 1. Load raw SCF
*---------------------------------------------------------------*
cd C:\Users\logan\Documents\GitHub\class\scf
use p22i6.dta, clear

foreach v of varlist _all {
    rename `v' `=upper("`v'")'
}

*---------------------------------------------------------------*
* 2. Construct variables
*---------------------------------------------------------------*

* Total revolving credit card debt
gen ccbal = X413 + X421 + X427

* High-interest indicator
gen highint = (X7132 > 5) if !missing(X7132)

* Card ownership
gen has_card = (X411 > 0 | X419 > 0 | X425 > 0)

* Banked indicator
gen banked = (X3501 == 1)

*---------------------------------------------------------------*
* 3. Sample restrictions
*---------------------------------------------------------------*

* Drop top 5% total income (NETWORTH summary var is not in codebk2022 main file)
xtile nw_pct = X5729 [pw=X42001], nq(100)
drop if nw_pct > 95

* Drop retirees
drop if X14 >= 68
drop if X14 <= 22

* Drop unbanked / underbanked
drop if banked==0
drop if has_card==0

*---------------------------------------------------------------*
* 4. Revolving high-interest debt
*---------------------------------------------------------------*
* Any revolving credit card debt (balance > 0)
gen revolve_any = (ccbal > 0) if !missing(ccbal)

* High-interest revolving debt indicator (main-card rate > 5%)
gen revolve = (ccbal > 0 & highint==1)

svyset [pw=X42001]

* Share with any revolving debt among households with a credit card
svy, subpop(if has_card==1): mean revolve_any

* Existing outcome: share with high-interest revolving debt
svy: mean revolve

* X7132 is PERCENT*100; include -1 ("no interest") as 0%, exclude 0 (inap.)
gen apr_main_pct = .
replace apr_main_pct = 0 if X7132 == -1
replace apr_main_pct = X7132/100 if X7132 > 0

preserve
keep if has_card==1 & revolve_any==1 & X7132!=0 & !missing(apr_main_pct, X42001)
sort apr_main_pct
gen w_cum = sum(X42001)
quietly summarize X42001, meanonly
scalar w_half = r(sum)/2
quietly summarize apr_main_pct if w_cum >= w_half, meanonly
display "Weighted median APR (percent) = " %6.2f r(min)
restore

