use "../data/compustat_output1.dta", clear

********************************************************************************
* Sample Restrictions
********************************************************************************

* Keep firms with total assets greater than $1 million
keep if at > 1 & !missing(at)

* Sort for panel operations
sort gvkey fyear

* Generate indicator for presence in the five years preceding 2001 (1996-2000) and 2001
* and five years preceding 2008 (2003-2007) and 2008
gen pre2002 = (fyear >= 1996 & fyear <= 2001)
gen pre2009 = (fyear >= 2003 & fyear <= 2008)

* Count number of years each firm appears in each window
bysort gvkey: egen years_pre2002 = total(pre2002)
bysort gvkey: egen years_pre2009 = total(pre2009)

* Keep firms present in at least one year in either window
* (Adjust threshold if "present in the five years" means all five years)
keep if (years_pre2002 == 6) | (years_pre2009 == 6)

drop pre2002 pre2009 years_pre2002 years_pre2009

tempfile compustat
save `compustat', replace

* ============================================================
* Step 1: Generate fiscal year-end date in CRSP monthly data
* ============================================================

use "../data/crsp.dta", clear

* Convert date to monthly (mofd) for merging
gen modate = mofd(date)
format modate %tm

rename PERMNO permno

tempfile crsp
save `crsp', replace

* ============================================================
* Step 2: Prepare Compustat file
* ============================================================

use "../data/compustat_output1.dta", clear

* Generate year-month of fiscal year end from datadate
gen modate = mofd(datadate)
format modate %tm

destring gvkey, replace

tempfile compustat
save `compustat', replace

* ============================================================
* Step 3: Prepare CCM link table
* ============================================================

preserve
use "../data/ccm_link.dta", clear

gen linkend = LINKENDDT
replace linkend = td(31dec9999) if missing(linkend)

rename LPERMNO permno
rename LINKDT  linkdt
rename GVKEY   gvkey

keep permno linkend linkdt gvkey

* Convert to monthly dates for annual comparison
gen linkdt_m  = mofd(linkdt)
gen linkend_m = mofd(linkend)
format linkdt_m linkend_m %tm

* Keep only valid link types if link type variable is present
* keep if inlist(LINKTYPE, "LC", "LU")  // uncomment if available

destring gvkey, replace
drop linkdt linkend

tempfile ccm
save `ccm', replace
restore

* ============================================================
* Step 4: Merge Compustat with CCM link table via CUSIP/PERMNO
* ============================================================

use `compustat', clear

* Merge on cusip to get permno candidates, then validate with link dates
* Use joinby on permno after merging link table on gvkey

* First join Compustat to link table on gvkey
joinby gvkey using `ccm'

* Keep only obs where datadate falls within link validity window
keep if modate >= linkdt_m & modate <= linkend_m
drop linkdt_m linkend_m

* ============================================================
* Step 5: Merge with CRSP monthly to get PRC and SHROUT
* ============================================================

joinby permno using `crsp'

* Keep only the CRSP observation matching the fiscal year-end month
keep if modate == mofd(datadate)

* Clean up price (CRSP uses negative for bid-ask midpoint)
replace PRC = abs(PRC)

* Convert SHROUT from thousands to millions to match Compustat CSHO units
gen csho_crsp  = SHROUT / 1000
gen prcc_f     = PRC

* ============================================================
* Step 6: Deduplicate and finalize
* ============================================================

duplicates drop gvkey datadate, force
isid gvkey datadate
sort gvkey datadate

save "../data/ccm_merged.dta"