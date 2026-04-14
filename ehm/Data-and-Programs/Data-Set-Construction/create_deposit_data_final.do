*************************************************************************
*****     Financial Distress, Retail Banking and FDIC Insurance     *****
*****                     Create Deposit Data File                  *****
*****                            09-10-2016                         *****
*************************************************************************

* Import, clean and combine data for deposit analysis
* 1) Import deposit level and bank balance sheet data from SDI
* 2) Import CDS to hazard rate data
* 3) Import deposit rate data from ratewatch files
* 4) Import CDS Data 

cd "/mnt/ide0/home/egan/Documents/Research/Bank Discipline/Analysis/Finalized Analysis"
set matsize 1000
set more off

********************************************************************************

***************************************
****    (1) Import Deposit Data    ****
***************************************
cd "\\ad.umn.edu\csom\users\eganx076\My Documents\Research\Deposit Competition\AER Version"
*Import bank deposit and balance/income data from the Statistics on Depository Institutions
use "\\ad.umn.edu\csom\users\eganx076\My Documents\Research\Data Sets\Statistics on Depository Institutions\q12002-q32013\Cleaned Data\SDI_Total_Deposits.dta",clear
merge 1:1 repdte date cert using "\\ad.umn.edu\csom\users\eganx076\My Documents\Research\Data Sets\Statistics on Depository Institutions\q12002-q32013\Cleaned Data\SDI_Income_and_Expense.dta", keep(1 3) nogen
merge 1:1 repdte date cert using "\\ad.umn.edu\csom\users\eganx076\My Documents\Research\Data Sets\Statistics on Depository Institutions\q12002-q32013\Cleaned Data\SDI_Assets_and_Liabilitites.dta", keep(1 3) nogen
merge 1:1 repdte date cert using "\\ad.umn.edu\csom\users\eganx076\My Documents\Research\Data Sets\Statistics on Depository Institutions\q12002-q32013\Cleaned Data\SDI_Loan_Charge-Offs.dta", keep(1 3) nogen
merge 1:1 repdte date cert using "\\ad.umn.edu\csom\users\eganx076\My Documents\Research\Data Sets\Statistics on Depository Institutions\q12002-q32013\Cleaned Data\SDI_Noncurrent_Loans.dta", keep(1 3) nogen
merge 1:1 repdte date cert using "\\ad.umn.edu\csom\users\eganx076\My Documents\Research\Data Sets\Statistics on Depository Institutions\q12002-q32013\Cleaned Data\SDI_Securities.dta", keep(1 3) nogen
*Calculate Deposits and Shares by Quarter
gen month= mofd(date)
format month %tm
gen depunins = dep - depins
bysort month: egen t_depins = sum(depins)
bysort month: egen t_depunins = sum(depunins)
gen share_unins = depunins/t_depunins
gen share_ins = depins/t_depins
* Keep relevant sample banks only
rename cert certno
keep if (certno==3510)|(certno==628)|(certno==33869)|(certno==3511)|(certno==6548)|(certno==7213)|(certno==6384)|(certno==57957)|(certno==867)|(certno==18409)|(certno==6672)|(certno==12368)|(certno==57890)|(certno==9846)|(certno==588)|(certno==17534)|(certno==29950)
bysort month: egen g_depins = sum(depins)
bysort month: egen g_depunins = sum(depunins)
gen o_ins_share = 1- g_depins/t_depins
gen o_unins_share = 1- g_depunins/t_depunins
gen shortname=""
replace shortname="BofA" if certno==3510
replace shortname="JPM" if certno==628
replace shortname="Wachovia" if certno==33869
replace shortname="Wells" if certno==3511
replace shortname="US" if certno==6548
replace shortname="Citi" if certno==7213
replace shortname="PNC" if certno==6384
replace shortname="RBS" if certno==57957
replace shortname="Suntrust" if certno==867
replace shortname="TD" if certno==18409
replace shortname="Fifth" if certno==6672
replace shortname="Regions" if certno==12368
replace shortname="HSBC" if certno==57890
replace shortname="BB&T" if certno==9846
replace shortname="M&T" if certno==588
replace shortname="NatCity" if certno==6557
replace shortname="Key" if certno==17534
replace shortname="Santander" if certno==29950
save "deposit_data_final.dta", replace

********************************************************************************

*********************************************
****    (2) Import CDS to Hazard Data    ****
*********************************************
*Note cds to hazard.csv is created in the corresponding .R file based on Hull 2012.
insheet using "cds to hazard.csv", clear comma
save "cds to hazard.dta", replace



*******************************************
****    (3) Import Deposit Rate Data    ***
*******************************************
use "/share/AEH/Data Files/Rate-Watch Data/Deposit Data/deposit_rate_data_2001_12.dta",clear 
keep if (charter_nbr==3510)|(charter_nbr==628)|(charter_nbr==33869)|(charter_nbr==3511)|(charter_nbr==6548)|(charter_nbr==7213)|(charter_nbr==6384)|(charter_nbr==57957)|(charter_nbr==867)|(charter_nbr==18409)|(charter_nbr==6672)|(charter_nbr==12368)|(charter_nbr==57890)|(charter_nbr==9846)|(charter_nbr==588)|(charter_nbr==6557)|(charter_nbr==17534)|(charter_nbr==29950)
keep if (productname=="01MCD100K")|(productname=="01MCD10K")|(productname=="03MCD100K")|(productname=="03MCD10K")|(productname=="06MCD100K")|(productname=="06MCD10K")|(productname=="09MCD100K")|(productname=="09MCD10K")|(productname=="12MCD100K")|(productname=="12MCD10K")|(productname=="18MCD100K")|(productname=="18MCD10K")|(productname=="24MCD100K")|(productname=="24MCD10K")|(productname=="36MCD100K")|(productname=="36MCD10K")|(productname=="48MCD100K")|(productname=="48MCD10K")|(productname=="60MCD100K")|(productname=="60MCD10K")
gen month= mofd(date(datesurveyed,"YMD"))
collapse (p50) rate, by(charter_nbr month productname)
forvalues x = 2002(1)2013{
	forvalues y = 1(1)12{
		if (`y'>10 & `x'==2013){
		}
		else if (`y'<10){
			append using "/share/AEH/Data Files/Rate-Watch Data/Deposit Data/deposit_rate_data_`x'_0`y'.dta"
			keep if (charter_nbr==3510)|(charter_nbr==628)|(charter_nbr==33869)|(charter_nbr==3511)|(charter_nbr==6548)|(charter_nbr==7213)|(charter_nbr==6384)|(charter_nbr==57957)|(charter_nbr==867)|(charter_nbr==18409)|(charter_nbr==6672)|(charter_nbr==12368)|(charter_nbr==57890)|(charter_nbr==9846)|(charter_nbr==588)|(charter_nbr==6557)|(charter_nbr==17534)|(charter_nbr==29950)
			keep if (productname=="01MCD100K")|(productname=="01MCD10K")|(productname=="03MCD100K")|(productname=="03MCD10K")|(productname=="06MCD100K")|(productname=="06MCD10K")|(productname=="09MCD100K")|(productname=="09MCD10K")|(productname=="12MCD100K")|(productname=="12MCD10K")|(productname=="18MCD100K")|(productname=="18MCD10K")|(productname=="24MCD100K")|(productname=="24MCD10K")|(productname=="36MCD100K")|(productname=="36MCD10K")|(productname=="48MCD100K")|(productname=="48MCD10K")|(productname=="60MCD100K")|(productname=="60MCD10K")
			replace month= mofd(date(datesurveyed,"YMD")) if month==.
			collapse (p50) rate, by(charter_nbr month productname)
		}
		else if (`y'>9){
			append using "/share/AEH/Data Files/Rate-Watch Data/Deposit Data/deposit_rate_data_`x'_`y'.dta"
			keep if (charter_nbr==3510)|(charter_nbr==628)|(charter_nbr==33869)|(charter_nbr==3511)|(charter_nbr==6548)|(charter_nbr==7213)|(charter_nbr==6384)|(charter_nbr==57957)|(charter_nbr==867)|(charter_nbr==18409)|(charter_nbr==6672)|(charter_nbr==12368)|(charter_nbr==57890)|(charter_nbr==9846)|(charter_nbr==588)|(charter_nbr==6557)|(charter_nbr==17534)|(charter_nbr==29950)
			keep if (productname=="01MCD100K")|(productname=="01MCD10K")|(productname=="03MCD100K")|(productname=="03MCD10K")|(productname=="06MCD100K")|(productname=="06MCD10K")|(productname=="09MCD100K")|(productname=="09MCD10K")|(productname=="12MCD100K")|(productname=="12MCD10K")|(productname=="18MCD100K")|(productname=="18MCD10K")|(productname=="24MCD100K")|(productname=="24MCD10K")|(productname=="36MCD100K")|(productname=="36MCD10K")|(productname=="48MCD100K")|(productname=="48MCD10K")|(productname=="60MCD100K")|(productname=="60MCD10K")
			replace month= mofd(date(datesurveyed,"YMD")) if month==.
			collapse (p50) rate, by(charter_nbr month productname)
		}
	}
}
format month %tm
gen shortname=""
replace shortname="BofA" if charter_nbr==3510
replace shortname="JPM" if charter_nbr==628
replace shortname="Wachovia" if charter_nbr==33869
replace shortname="Wells" if charter_nbr==3511
replace shortname="US" if charter_nbr==6548
replace shortname="Citi" if charter_nbr==7213
replace shortname="PNC" if charter_nbr==6384
replace shortname="RBS" if charter_nbr==57957
replace shortname="Suntrust" if charter_nbr==867
replace shortname="TD" if charter_nbr==18409
replace shortname="Fifth" if charter_nbr==6672
replace shortname="Regions" if charter_nbr==12368
replace shortname="HSBC" if charter_nbr==57890
replace shortname="BB&T" if charter_nbr==9846
replace shortname="M&T" if charter_nbr==588
replace shortname="NatCity" if charter_nbr==6557
replace shortname="Key" if charter_nbr==17534
replace shortname="Santander" if charter_nbr==29950
preserve
gen large = strpos(productname,"100K")>0
gen length = real(substr(productname,1,2))
replace rate =rate/100
drop if length!=12
gen size = large+1
keep rate charter_nbr shortname month size
reshape wide rate, i(shortname charter_nbr month) j(size)
rename rate1 ins_rate
rename rate2 unins_rate
save "depost_rate_data_1yr_cd_wide.dta",replace 
restore
save "deposit_rate_data.dta",replace



********************************************************************************


********************************************************************************

***********************************
****    (4) Import CDS Data    ****
***********************************
use "bankcds.dta", clear
drop if spread1y ==. & spread5y ==.
*Create bank identifier variables
gen certno =0
replace certno=3510 if redcode=="0G656B" /*Bank of America*/
replace certno=628 if redcode=="4C933D" /*JPM*/
replace certno=628 if redcode=="4C93HG" /*JPM*/
replace certno=33869 if redcode=="9BBGDB" /*Wachovia*/
replace certno=3511 if redcode=="9DDH8V" /*Wells*/
replace certno=7213 if redcode=="PR9D89" /*Citi*/
replace certno=6384 if redcode=="6FC7BC" /*PNC*/
replace certno=57957 if redcode=="NUD88R" /*RBS*/
replace certno=867 if redcode=="8EDGA5" /*Suntrust*/
replace certno=18409 if redcode=="8HA277" /*TD*/
replace certno=6672 if redcode=="347DEA" /*Fifth Third*/
replace certno=12368 if redcode=="7CDHD4" /*Regions*/
replace certno=12368 if redcode=="7CEAF3" /*Regions*/
replace certno=57890 if redcode=="4E46A8" /*HSBC*/
replace certno=9846 if redcode=="059DCA" /*BB&T*/
replace certno=4297 if redcode=="1F444N" /*CapOne*/
replace certno=4297 if redcode=="1F445B" /*CapOne*/
replace certno=6557 if redcode=="6J8945" /*National City*/
replace certno=17534 if redcode=="4DC58D" /*KeyBank*/
replace certno=29950 if redcode=="8C9E95" /*Santander*/
replace certno =6548 if redcode== "992BGA" /*US Bank*/
drop if certno==0
replace shortname="BofA" if certno==3510
replace shortname="JPM" if certno==628
replace shortname="Wachovia" if certno==33869
replace shortname="Wells" if certno==3511
replace shortname="US" if certno==6548
replace shortname="Citi" if certno==7213
replace shortname="PNC" if certno==6384
replace shortname="RBS" if certno==57957
replace shortname="Suntrust" if certno==867
replace shortname="TD" if certno==18409
replace shortname="Fifth" if certno==6672
replace shortname="Regions" if certno==12368
replace shortname="HSBC" if certno==57890
replace shortname="BB&T" if certno==9846
replace shortname="M&T" if certno==588
replace shortname="CapOne" if certno==4297
replace shortname="NatCity" if certno==6557
replace shortname="Key" if certno==17534
replace shortname="Santander" if certno==29950
*Collapse at the monthly level
gen month =mofd(date)
collapse spread1y spread5y, by(certno shortname month)
save "monthly cds.dta",replace
