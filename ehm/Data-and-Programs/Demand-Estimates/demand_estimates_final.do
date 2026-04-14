*************************************************************************
*****     Financial Distress, Retail Banking and FDIC Insurance     *****
*****                        Deposit Analaysis                      *****
*****                            09-10-2016                         *****
*************************************************************************

* Import, clean and combine data for deposit analysis
* 1) Import data
* 2) Reduced form analysis Diff & Diff
* 3) Demand Estimates
* 4) Create data for calibration exercise

cd "\\ad.umn.edu\csom\users\eganx076\My Documents\Research\Deposit Competition\AER Version"
set matsize 1000
set more off

********************************************************************************

******************************
****    (1) Import Data    ***
******************************
use "deposit_data_final.dta",clear
gen eq_ratio = eq/dep
*Combine with CDS data 
merge m:1 shortname month using "\\ad.umn.edu\csom\users\eganx076\My Documents\Research\Deposit Competition\CDS Data\monthly cds.dta", keep(1 3) nogen
*Merge with CMS and CMT data
gen date2= date
replace date = dofm(month)
merge m:1 date using  "\\ad.umn.edu\csom\users\eganx076\My Documents\Research\Data Sets\FREDS Interest Rate Data\2013\monthly_swaps_and_treasuries_rates.dta", keep(1 3) nogen
drop date
rename date2 date
*Merge with CD rate data
gen charter_nbr=certno
merge 1:1 charter_nbr month using "depost_rate_data_1yr_cd_wide.dta", keep(1 3) nogen
*Convert cds to hazard
gen cds = round(spread5y,0.0001)
merge m:1 cds using "cds to hazard.dta", keep(matched) nogen 
*Create relevant variables
gen d_share_unins = ln(share_unins)-ln(o_unins_share)
gen d_share_ins = ln(share_ins)-ln(o_ins_share)
gen loan_exp = lnlsnet/asset
gen iv_ln1 = ntre/asset
gen iv_ln2= ntlnls/asset-ntre/asset
gen iv_cmo = sccol/sc
gen iv_cmo2= idsccmo/sc
gen inter_3510= cmt1y*(certno==3510)
gen inter_628= cmt1y*(certno==628)
gen inter_33869= cmt1y*(certno==33869)
gen inter_3511= cmt1y*(certno==3511)
gen inter_6548= cmt1y*(certno==6548)
gen inter_7213= cmt1y*(certno==7213)
gen inter_6384= cmt1y*(certno==6384)
gen inter_57957= cmt1y*(certno==57957)
gen inter_867= cmt1y*(certno==867)
gen inter_18409= cmt1y*(certno==18409)
gen inter_6672= cmt1y*(certno==6672)
gen inter_12368= cmt1y*(certno==12368)
gen inter_57890= cmt1y*(certno==57890)
gen inter_9846= cmt1y*(certno==9846)
gen inter_588= cmt1y*(certno==588)
gen inter_4297= cmt1y*(certno==4297)
gen inter_6557= cmt1y*(certno==6557)
gen inter_17534= cmt1y*(certno==17534)
gen inter_29950= cmt1y*(certno==29950)
gen q = qofd(date)
gen ins_rt_spd= ins_rate-cmt1y
gen unins_rt_spd= unins_rate-cmt1y
su depins depunins cds unins_rt_spd ins_rt_spd

********************************************************************************

****************************************************
****    (2) Reduced Form Diff-Diff Estimates    ****
****************************************************

* Create Table 2
gen share_diff = ln(share_unins)-ln(share_ins)
xi: reg share_diff hazard i.q i.certno offdom,robust
outreg2 using reduced_form, ctitle(Ins1) bdec(2) sdec(2) tdec(2) alpha(0.01,0.05,0.10) drop(_I*) tex replace
xi: reg d_share_unins hazard i.q i.certno offdom,robust
outreg2 using reduced_form, ctitle(Ins1) bdec(2) sdec(2) tdec(2) alpha(0.01,0.05,0.10) drop(_I*) tex append
xi: reg d_share_ins hazard i.q i.certno offdom,robust
outreg2 using reduced_form, ctitle(Ins1) bdec(2) sdec(2) tdec(2) alpha(0.01,0.05,0.10) drop(_I*) tex append

********************************************************************************

************************************
****    (3) Demand Estimates    ****
************************************

*Calculate Table 3

******************************
*** 3.a Uninsured Deposits ***
******************************
*Calculate passthrough IV 
gen rt= unins_rate-cmt1y
xi: reg unins_rate inter* i.q i.certno offdom [weight = sqrt(t_depunins)], robust
predict p_iv_unins_fe
*Demand estimates
xi, noomit: ivregress 2sls d_share_unins (rt hazard = p_iv_unins_fe iv_cmo) i.certno  i.q  offdom  [weight = sqrt(t_depunins)], noconstant robust
outreg2 using demand_estimates, ctitle(Ins1) bdec(2) sdec(2) tdec(2) alpha(0.01,0.05,0.10) drop(_I*) tex replace
xi, noomit: ivregress 2sls d_share_unins (rt hazard = p_iv_unins_fe iv_ln1 iv_ln2) i.certno  i.q  offdom  [weight = sqrt(t_depunins)], noconstant robust
outreg2 using demand_estimates, ctitle(Ins1) bdec(2) sdec(2) tdec(2) alpha(0.01,0.05,0.10) drop(_I*) tex append
xi, noomit: ivregress 2sls d_share_unins (rt hazard = p_iv_unins_fe iv_ln1 iv_ln2 iv_cmo) i.certno  i.q  offdom  [weight = sqrt(t_depunins)], noconstant robust
outreg2 using demand_estimates, ctitle(Ins1) bdec(2) sdec(2) tdec(2) alpha(0.01,0.05,0.10) drop(_I*) tex append
matrix coeff = e(b)
predict unins_error, resid

****************************
*** 3.b Insured Deposits ***
****************************
*Calculate passthrough IV 
replace rt= (ins_rate - cmt1y)
xi: reg ins_rate  inter* i.q i.certno offdom  [weight = sqrt(t_depins)], robust
predict p_iv_ins_fe
*Demand estimates
xi, noomit: ivregress 2sls d_share_ins (rt = p_iv_ins_fe) i.certno  i.q offdom  [weight = sqrt(t_depins)], noconstant robust 
outreg2 using demand_estimates, ctitle(Ins1) bdec(2) sdec(2) tdec(2) alpha(0.01,0.05,0.10) drop(_I*) tex append
matrix coeff2 = e(b)
predict ins_error, resid
xi, noomit: ivregress 2sls d_share_ins (rt hazard = p_iv_ins_fe iv_cmo) i.certno  i.q  offdom  [weight = sqrt(t_depins)], noconstant robust
outreg2 using demand_estimates, ctitle(Ins1) bdec(2) sdec(2) tdec(2) alpha(0.01,0.05,0.10) drop(_I*) tex append
xi, noomit: ivregress 2sls d_share_ins (rt hazard = p_iv_ins_fe iv_ln1 iv_ln2) i.certno  i.q  offdom  [weight = sqrt(t_depins)], noconstant robust
outreg2 using demand_estimates, ctitle(Ins1) bdec(2) sdec(2) tdec(2) alpha(0.01,0.05,0.10) drop(_I*) tex append
xi, noomit: ivregress 2sls d_share_ins (rt hazard = p_iv_ins_fe iv_ln1 iv_ln2 iv_cmo) i.certno  i.q  offdom  [weight = sqrt(t_depins)], noconstant robust
outreg2 using demand_estimates, ctitle(Ins1) bdec(2) sdec(2) tdec(2) alpha(0.01,0.05,0.10) drop(_I*) tex append


***********************************
*** 3.c Store Demand Parameters ***
***********************************
*Uninsured brand effects
gen date_u = coeff[1,47]
gen delta_u = 0
replace delta_u= coeff[1,3] if certno==628
replace delta_u= coeff[1,4] if certno==867
replace delta_u= coeff[1,5] if certno==3510
replace delta_u= coeff[1,6] if certno==3511
replace delta_u= coeff[1,7] if certno==6384
replace delta_u= coeff[1,8] if certno==6548
replace delta_u= coeff[1,9] if certno==6672
replace delta_u= coeff[1,10] if certno==7213
replace delta_u= coeff[1,11] if certno==9846
replace delta_u= coeff[1,12] if certno==12368
replace delta_u= coeff[1,13] if certno==17534
replace delta_u= coeff[1,14] if certno==18409
replace delta_u= coeff[1,15] if certno==29950
replace delta_u= coeff[1,16] if certno==33869
replace delta_u= coeff[1,17] if certno==57890
replace delta_u= coeff[1,18] if certno==57957
gen office_u = coeff[1,19]
*Insured brand effects
gen date_i = coeff2[1,46]
gen delta_i = 0
replace delta_i = coeff2[1,2] if certno==628
replace delta_i = coeff2[1,3] if certno==867
replace delta_i = coeff2[1,4] if certno==3510
replace delta_i = coeff2[1,5] if certno==3511
replace delta_i = coeff2[1,6] if certno==6384
replace delta_i = coeff2[1,7] if certno==6548
replace delta_i = coeff2[1,8] if certno==6672
replace delta_i = coeff2[1,9] if certno==7213
replace delta_i = coeff2[1,10] if certno==9846
replace delta_i = coeff2[1,11] if certno==12368
replace delta_i = coeff2[1,12] if certno==17534
replace delta_i = coeff2[1,13] if certno==18409
replace delta_i = coeff2[1,14] if certno==29950
replace delta_i = coeff2[1,15] if certno==33869
replace delta_i = coeff2[1,16] if certno==57890
replace delta_i = coeff2[1,17] if certno==57957
gen office_i = coeff2[1,18]
*Elasticities
gen alpha_unins = coeff[1,1]
gen alpha_ins = coeff2[1,1]
gen gamma = coeff[1,2]

********************************************************************************

*******************************************
****    (4) Create Calibration Data    ****
*******************************************
gen unins_share = share_unins
gen ins_share = share_ins
gen ins_dep= depins	
gen unins_dep= depunins
gen credit_spread_t= spread5y+cmt10y if date==date("March 31, 2007","MDY")
bysort certno: egen credit_spread= max(credit_spread_t)
gen b = (liab-dep)*(credit_spread)
gen int_rate=0.05
gen sd = (1+int_rate)*(((b+ins_dep*(unins_rate+1/(alpha_unins*(1-unins_share))-1/(alpha_ins*(1-ins_share)))+unins_dep*unins_rate)/(ins_dep+unins_dep))-1/(alpha_unins*(1-unins_share))-unins_rate)/(normalden(invnormal(hazard))+invnormal(hazard)*(hazard+int_rate)-(1+int_rate)*normalden(invnormal(hazard))/(1-hazard))
gen mu = 1/(alpha_unins*(1-unins_share))+unins_rate-sd*normalden(invnormal(hazard))/(1-hazard)
gen c = -1/((1-ins_share)*alpha_ins) + mu + sd*normalden(invnormal(hazard))/(1-hazard)-ins_rate
gen r_bar = mu+invnormal(hazard)*sd
gen lhs = (1+int_rate)*(b-ins_dep*(r_bar-c-ins_rate)-unins_dep*(r_bar-unins_rate))
gen rhs = (ins_dep+unins_dep)*((1-hazard)*mu+sd*normalden(invnormal(hazard))-(1-hazard)*r_bar)
gen double M_i  = t_depins
gen  double M_u = t_depunins
*Export data
preserve
format date %td
keep shortname int_rate b hazard mu sd r_bar  delta_u* delta_i* ins_dep unins_dep ins_share unins_share alpha_unins alpha_ins ins_rate unins_rate ins_error unins_error c gamma o_unins_share o_ins_share date_i date_u office_u office_i offdom cmt1y M_i M_u date eq_ratio
outsheet using "calibration params final.csv",replace comma
restore
