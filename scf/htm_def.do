gen chk_sav_lt_1_24 = (checking + saving + mma + call + cds + prepaid) / income < (1/24) if income > 0 & !missing(checking, saving, income)
mean chk_sav_lt_1_24 [pw=wgt]
drop chk_sav_lt_1_24

gen all_liq = (checking + saving + mma + call + cds + prepaid + stocks) / income < (1/24) if income > 0 & !missing(checking, saving, income)
mean all_liq [pw=wgt]
drop all_liq
