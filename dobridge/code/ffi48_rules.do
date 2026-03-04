* ---- BEGIN GENERATED RULES ----
* Assumes tempvar sicnum exists and output vars `gen' and `desc' exist.
* Guarded with missing(`gen') to mimic SAS else-if (first match wins).

replace `gen'  = 1     if `touse' & missing(`gen') & inrange(`sicnum', 100, 199)
replace `desc' = "Agric" if `touse' & `gen'==1 & `desc'==""

replace `gen'  = 1     if `touse' & missing(`gen') & inrange(`sicnum', 200, 299)
replace `desc' = "Agric" if `touse' & `gen'==1 & `desc'==""

replace `gen'  = 1     if `touse' & missing(`gen') & inrange(`sicnum', 700, 799)
replace `desc' = "Agric" if `touse' & `gen'==1 & `desc'==""

replace `gen'  = 1     if `touse' & missing(`gen') & inrange(`sicnum', 910, 919)
replace `desc' = "Agric" if `touse' & `gen'==1 & `desc'==""

replace `gen'  = 1     if `touse' & missing(`gen') & `sicnum'==2048
replace `desc' = "Agric" if `touse' & `gen'==1 & `desc'==""

replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2000, 2009)
replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2010, 2019)
replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2020, 2029)
replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2030, 2039)
replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2040, 2046)
replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2050, 2059)
replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2060, 2063)
replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2070, 2079)
replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2090, 2092)
replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

replace `gen'  = 2     if `touse' & missing(`gen') & `sicnum'==2095
replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2098, 2099)
replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

replace `gen'  = 3     if `touse' & missing(`gen') & inrange(`sicnum', 2064, 2068)
replace `desc' = "Soda" if `touse' & `gen'==3 & `desc'==""

replace `gen'  = 3     if `touse' & missing(`gen') & `sicnum'==2086
replace `desc' = "Soda" if `touse' & `gen'==3 & `desc'==""

replace `gen'  = 3     if `touse' & missing(`gen') & `sicnum'==2087
replace `desc' = "Soda" if `touse' & `gen'==3 & `desc'==""

replace `gen'  = 3     if `touse' & missing(`gen') & `sicnum'==2096
replace `desc' = "Soda" if `touse' & `gen'==3 & `desc'==""

replace `gen'  = 3     if `touse' & missing(`gen') & `sicnum'==2097
replace `desc' = "Soda" if `touse' & `gen'==3 & `desc'==""

replace `gen'  = 4     if `touse' & missing(`gen') & `sicnum'==2080
replace `desc' = "Beer" if `touse' & `gen'==4 & `desc'==""

replace `gen'  = 4     if `touse' & missing(`gen') & `sicnum'==2082
replace `desc' = "Beer" if `touse' & `gen'==4 & `desc'==""

replace `gen'  = 4     if `touse' & missing(`gen') & `sicnum'==2083
replace `desc' = "Beer" if `touse' & `gen'==4 & `desc'==""

replace `gen'  = 4     if `touse' & missing(`gen') & `sicnum'==2084
replace `desc' = "Beer" if `touse' & `gen'==4 & `desc'==""

replace `gen'  = 4     if `touse' & missing(`gen') & `sicnum'==2085
replace `desc' = "Beer" if `touse' & `gen'==4 & `desc'==""

replace `gen'  = 5     if `touse' & missing(`gen') & inrange(`sicnum', 2100, 2199)
replace `desc' = "Smoke" if `touse' & `gen'==5 & `desc'==""

replace `gen'  = 6     if `touse' & missing(`gen') & inrange(`sicnum', 920, 999)
replace `desc' = "Toys" if `touse' & `gen'==6 & `desc'==""

replace `gen'  = 6     if `touse' & missing(`gen') & inrange(`sicnum', 3650, 3651)
replace `desc' = "Toys" if `touse' & `gen'==6 & `desc'==""

replace `gen'  = 6     if `touse' & missing(`gen') & `sicnum'==3652
replace `desc' = "Toys" if `touse' & `gen'==6 & `desc'==""

replace `gen'  = 6     if `touse' & missing(`gen') & `sicnum'==3732
replace `desc' = "Toys" if `touse' & `gen'==6 & `desc'==""

replace `gen'  = 6     if `touse' & missing(`gen') & inrange(`sicnum', 3930, 3931)
replace `desc' = "Toys" if `touse' & `gen'==6 & `desc'==""

replace `gen'  = 6     if `touse' & missing(`gen') & inrange(`sicnum', 3940, 3949)
replace `desc' = "Toys" if `touse' & `gen'==6 & `desc'==""

replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7800, 7829)
replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7830, 7833)
replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7840, 7841)
replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

replace `gen'  = 7     if `touse' & missing(`gen') & `sicnum'==7900
replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7910, 7911)
replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7920, 7929)
replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7930, 7933)
replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7940, 7949)
replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

replace `gen'  = 7     if `touse' & missing(`gen') & `sicnum'==7980
replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7990, 7999)
replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2700, 2709)
replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2710, 2719)
replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2720, 2729)
replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2730, 2739)
replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2740, 2749)
replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2770, 2771)
replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2780, 2789)
replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2790, 2799)
replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==2047
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 2391, 2392)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 2510, 2519)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 2590, 2599)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 2840, 2843)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==2844
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3160, 3161)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3170, 3171)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3172
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3190, 3199)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3229
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3260
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3262, 3263)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3269
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3230, 3231)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3630, 3639)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3750, 3751)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3800
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3860, 3861)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3870, 3873)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3910, 3911)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3914
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3915
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3960, 3962)
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3991
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3995
replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

replace `gen'  = 10     if `touse' & missing(`gen') & inrange(`sicnum', 2300, 2390)
replace `desc' = "Clths" if `touse' & `gen'==10 & `desc'==""

replace `gen'  = 10     if `touse' & missing(`gen') & inrange(`sicnum', 3020, 3021)
replace `desc' = "Clths" if `touse' & `gen'==10 & `desc'==""

replace `gen'  = 10     if `touse' & missing(`gen') & inrange(`sicnum', 3100, 3111)
replace `desc' = "Clths" if `touse' & `gen'==10 & `desc'==""

replace `gen'  = 10     if `touse' & missing(`gen') & inrange(`sicnum', 3130, 3131)
replace `desc' = "Clths" if `touse' & `gen'==10 & `desc'==""

replace `gen'  = 10     if `touse' & missing(`gen') & inrange(`sicnum', 3140, 3149)
replace `desc' = "Clths" if `touse' & `gen'==10 & `desc'==""

replace `gen'  = 10     if `touse' & missing(`gen') & inrange(`sicnum', 3150, 3151)
replace `desc' = "Clths" if `touse' & `gen'==10 & `desc'==""

replace `gen'  = 10     if `touse' & missing(`gen') & inrange(`sicnum', 3963, 3965)
replace `desc' = "Clths" if `touse' & `gen'==10 & `desc'==""

replace `gen'  = 11     if `touse' & missing(`gen') & inrange(`sicnum', 8000, 8099)
replace `desc' = "Hlth" if `touse' & `gen'==11 & `desc'==""

replace `gen'  = 12     if `touse' & missing(`gen') & `sicnum'==3693
replace `desc' = "MedEq" if `touse' & `gen'==12 & `desc'==""

replace `gen'  = 12     if `touse' & missing(`gen') & inrange(`sicnum', 3840, 3849)
replace `desc' = "MedEq" if `touse' & `gen'==12 & `desc'==""

replace `gen'  = 12     if `touse' & missing(`gen') & inrange(`sicnum', 3850, 3851)
replace `desc' = "MedEq" if `touse' & `gen'==12 & `desc'==""

replace `gen'  = 13     if `touse' & missing(`gen') & `sicnum'==2830
replace `desc' = "Drugs" if `touse' & `gen'==13 & `desc'==""

replace `gen'  = 13     if `touse' & missing(`gen') & `sicnum'==2831
replace `desc' = "Drugs" if `touse' & `gen'==13 & `desc'==""

replace `gen'  = 13     if `touse' & missing(`gen') & `sicnum'==2833
replace `desc' = "Drugs" if `touse' & `gen'==13 & `desc'==""

replace `gen'  = 13     if `touse' & missing(`gen') & `sicnum'==2834
replace `desc' = "Drugs" if `touse' & `gen'==13 & `desc'==""

replace `gen'  = 13     if `touse' & missing(`gen') & `sicnum'==2835
replace `desc' = "Drugs" if `touse' & `gen'==13 & `desc'==""

replace `gen'  = 13     if `touse' & missing(`gen') & `sicnum'==2836
replace `desc' = "Drugs" if `touse' & `gen'==13 & `desc'==""

replace `gen'  = 14     if `touse' & missing(`gen') & inrange(`sicnum', 2800, 2809)
replace `desc' = "Chems" if `touse' & `gen'==14 & `desc'==""

replace `gen'  = 14     if `touse' & missing(`gen') & inrange(`sicnum', 2810, 2819)
replace `desc' = "Chems" if `touse' & `gen'==14 & `desc'==""

replace `gen'  = 14     if `touse' & missing(`gen') & inrange(`sicnum', 2820, 2829)
replace `desc' = "Chems" if `touse' & `gen'==14 & `desc'==""

replace `gen'  = 14     if `touse' & missing(`gen') & inrange(`sicnum', 2850, 2859)
replace `desc' = "Chems" if `touse' & `gen'==14 & `desc'==""

replace `gen'  = 14     if `touse' & missing(`gen') & inrange(`sicnum', 2860, 2869)
replace `desc' = "Chems" if `touse' & `gen'==14 & `desc'==""

replace `gen'  = 14     if `touse' & missing(`gen') & inrange(`sicnum', 2870, 2879)
replace `desc' = "Chems" if `touse' & `gen'==14 & `desc'==""

replace `gen'  = 14     if `touse' & missing(`gen') & inrange(`sicnum', 2890, 2899)
replace `desc' = "Chems" if `touse' & `gen'==14 & `desc'==""

replace `gen'  = 15     if `touse' & missing(`gen') & `sicnum'==3031
replace `desc' = "Rubbr" if `touse' & `gen'==15 & `desc'==""

replace `gen'  = 15     if `touse' & missing(`gen') & `sicnum'==3041
replace `desc' = "Rubbr" if `touse' & `gen'==15 & `desc'==""

replace `gen'  = 15     if `touse' & missing(`gen') & inrange(`sicnum', 3050, 3053)
replace `desc' = "Rubbr" if `touse' & `gen'==15 & `desc'==""

replace `gen'  = 15     if `touse' & missing(`gen') & inrange(`sicnum', 3060, 3069)
replace `desc' = "Rubbr" if `touse' & `gen'==15 & `desc'==""

replace `gen'  = 15     if `touse' & missing(`gen') & inrange(`sicnum', 3070, 3079)
replace `desc' = "Rubbr" if `touse' & `gen'==15 & `desc'==""

replace `gen'  = 15     if `touse' & missing(`gen') & inrange(`sicnum', 3080, 3089)
replace `desc' = "Rubbr" if `touse' & `gen'==15 & `desc'==""

replace `gen'  = 15     if `touse' & missing(`gen') & inrange(`sicnum', 3090, 3099)
replace `desc' = "Rubbr" if `touse' & `gen'==15 & `desc'==""

replace `gen'  = 16     if `touse' & missing(`gen') & inrange(`sicnum', 2200, 2269)
replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

replace `gen'  = 16     if `touse' & missing(`gen') & inrange(`sicnum', 2270, 2279)
replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

replace `gen'  = 16     if `touse' & missing(`gen') & inrange(`sicnum', 2280, 2284)
replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

replace `gen'  = 16     if `touse' & missing(`gen') & inrange(`sicnum', 2290, 2295)
replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

replace `gen'  = 16     if `touse' & missing(`gen') & `sicnum'==2297
replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

replace `gen'  = 16     if `touse' & missing(`gen') & `sicnum'==2298
replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

replace `gen'  = 16     if `touse' & missing(`gen') & `sicnum'==2299
replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

replace `gen'  = 16     if `touse' & missing(`gen') & inrange(`sicnum', 2393, 2395)
replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

replace `gen'  = 16     if `touse' & missing(`gen') & inrange(`sicnum', 2397, 2399)
replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 800, 899)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 2400, 2439)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 2450, 2459)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 2490, 2499)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 2660, 2661)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 2950, 2952)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3200
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3210, 3211)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3240, 3241)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3250, 3259)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3261
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3264
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3270, 3275)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3280, 3281)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3290, 3293)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3295, 3299)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3420, 3429)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3430, 3433)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3440, 3441)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3442
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3446
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3448
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3449
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3450, 3451)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3452
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3490, 3499)
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3996
replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

replace `gen'  = 18     if `touse' & missing(`gen') & inrange(`sicnum', 1500, 1511)
replace `desc' = "Cnstr" if `touse' & `gen'==18 & `desc'==""

replace `gen'  = 18     if `touse' & missing(`gen') & inrange(`sicnum', 1520, 1529)
replace `desc' = "Cnstr" if `touse' & `gen'==18 & `desc'==""

replace `gen'  = 18     if `touse' & missing(`gen') & inrange(`sicnum', 1530, 1539)
replace `desc' = "Cnstr" if `touse' & `gen'==18 & `desc'==""

replace `gen'  = 18     if `touse' & missing(`gen') & inrange(`sicnum', 1540, 1549)
replace `desc' = "Cnstr" if `touse' & `gen'==18 & `desc'==""

replace `gen'  = 18     if `touse' & missing(`gen') & inrange(`sicnum', 1600, 1699)
replace `desc' = "Cnstr" if `touse' & `gen'==18 & `desc'==""

replace `gen'  = 18     if `touse' & missing(`gen') & inrange(`sicnum', 1700, 1799)
replace `desc' = "Cnstr" if `touse' & `gen'==18 & `desc'==""

replace `gen'  = 19     if `touse' & missing(`gen') & `sicnum'==3300
replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3310, 3317)
replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3320, 3325)
replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3330, 3339)
replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3340, 3341)
replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3350, 3357)
replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3360, 3369)
replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3370, 3379)
replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3390, 3399)
replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

replace `gen'  = 20     if `touse' & missing(`gen') & `sicnum'==3400
replace `desc' = "FabPr" if `touse' & `gen'==20 & `desc'==""

replace `gen'  = 20     if `touse' & missing(`gen') & `sicnum'==3443
replace `desc' = "FabPr" if `touse' & `gen'==20 & `desc'==""

replace `gen'  = 20     if `touse' & missing(`gen') & `sicnum'==3444
replace `desc' = "FabPr" if `touse' & `gen'==20 & `desc'==""

replace `gen'  = 20     if `touse' & missing(`gen') & inrange(`sicnum', 3460, 3469)
replace `desc' = "FabPr" if `touse' & `gen'==20 & `desc'==""

replace `gen'  = 20     if `touse' & missing(`gen') & inrange(`sicnum', 3470, 3479)
replace `desc' = "FabPr" if `touse' & `gen'==20 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & inrange(`sicnum', 3510, 3519)
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & inrange(`sicnum', 3520, 3529)
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3530
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3531
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3532
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3533
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3534
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3535
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3536
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3538
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & inrange(`sicnum', 3540, 3549)
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & inrange(`sicnum', 3550, 3559)
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & inrange(`sicnum', 3560, 3569)
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3580
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3581
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3582
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3585
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3586
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3589
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 21     if `touse' & missing(`gen') & inrange(`sicnum', 3590, 3599)
replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

replace `gen'  = 22     if `touse' & missing(`gen') & `sicnum'==3600
replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

replace `gen'  = 22     if `touse' & missing(`gen') & inrange(`sicnum', 3610, 3613)
replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

replace `gen'  = 22     if `touse' & missing(`gen') & inrange(`sicnum', 3620, 3621)
replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

replace `gen'  = 22     if `touse' & missing(`gen') & inrange(`sicnum', 3623, 3629)
replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

replace `gen'  = 22     if `touse' & missing(`gen') & inrange(`sicnum', 3640, 3644)
replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

replace `gen'  = 22     if `touse' & missing(`gen') & `sicnum'==3645
replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

replace `gen'  = 22     if `touse' & missing(`gen') & `sicnum'==3646
replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

replace `gen'  = 22     if `touse' & missing(`gen') & inrange(`sicnum', 3648, 3649)
replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

replace `gen'  = 22     if `touse' & missing(`gen') & `sicnum'==3660
replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

replace `gen'  = 22     if `touse' & missing(`gen') & `sicnum'==3690
replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

replace `gen'  = 22     if `touse' & missing(`gen') & inrange(`sicnum', 3691, 3692)
replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

replace `gen'  = 22     if `touse' & missing(`gen') & `sicnum'==3699
replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==2296
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==2396
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & inrange(`sicnum', 3010, 3011)
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3537
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3647
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3694
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3700
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3710
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3711
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3713
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3714
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3715
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3716
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3792
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & inrange(`sicnum', 3790, 3791)
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3799
replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

replace `gen'  = 24     if `touse' & missing(`gen') & `sicnum'==3720
replace `desc' = "Aero" if `touse' & `gen'==24 & `desc'==""

replace `gen'  = 24     if `touse' & missing(`gen') & `sicnum'==3721
replace `desc' = "Aero" if `touse' & `gen'==24 & `desc'==""

replace `gen'  = 24     if `touse' & missing(`gen') & inrange(`sicnum', 3723, 3724)
replace `desc' = "Aero" if `touse' & `gen'==24 & `desc'==""

replace `gen'  = 24     if `touse' & missing(`gen') & `sicnum'==3725
replace `desc' = "Aero" if `touse' & `gen'==24 & `desc'==""

replace `gen'  = 24     if `touse' & missing(`gen') & inrange(`sicnum', 3728, 3729)
replace `desc' = "Aero" if `touse' & `gen'==24 & `desc'==""

replace `gen'  = 25     if `touse' & missing(`gen') & inrange(`sicnum', 3730, 3731)
replace `desc' = "Ships" if `touse' & `gen'==25 & `desc'==""

replace `gen'  = 25     if `touse' & missing(`gen') & inrange(`sicnum', 3740, 3743)
replace `desc' = "Ships" if `touse' & `gen'==25 & `desc'==""

replace `gen'  = 26     if `touse' & missing(`gen') & inrange(`sicnum', 3760, 3769)
replace `desc' = "Guns" if `touse' & `gen'==26 & `desc'==""

replace `gen'  = 26     if `touse' & missing(`gen') & `sicnum'==3795
replace `desc' = "Guns" if `touse' & `gen'==26 & `desc'==""

replace `gen'  = 26     if `touse' & missing(`gen') & inrange(`sicnum', 3480, 3489)
replace `desc' = "Guns" if `touse' & `gen'==26 & `desc'==""

replace `gen'  = 27     if `touse' & missing(`gen') & inrange(`sicnum', 1040, 1049)
replace `desc' = "Gold" if `touse' & `gen'==27 & `desc'==""

replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1000, 1009)
replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1010, 1019)
replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1020, 1029)
replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1030, 1039)
replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1050, 1059)
replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1060, 1069)
replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1070, 1079)
replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1080, 1089)
replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1090, 1099)
replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1100, 1119)
replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1400, 1499)
replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

replace `gen'  = 29     if `touse' & missing(`gen') & inrange(`sicnum', 1200, 1299)
replace `desc' = "Coal" if `touse' & `gen'==29 & `desc'==""

replace `gen'  = 30     if `touse' & missing(`gen') & `sicnum'==1300
replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

replace `gen'  = 30     if `touse' & missing(`gen') & inrange(`sicnum', 1310, 1319)
replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

replace `gen'  = 30     if `touse' & missing(`gen') & inrange(`sicnum', 1320, 1329)
replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

replace `gen'  = 30     if `touse' & missing(`gen') & inrange(`sicnum', 1330, 1339)
replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

replace `gen'  = 30     if `touse' & missing(`gen') & inrange(`sicnum', 1370, 1379)
replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

replace `gen'  = 30     if `touse' & missing(`gen') & `sicnum'==1380
replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

replace `gen'  = 30     if `touse' & missing(`gen') & `sicnum'==1381
replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

replace `gen'  = 30     if `touse' & missing(`gen') & `sicnum'==1382
replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

replace `gen'  = 30     if `touse' & missing(`gen') & `sicnum'==1389
replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

replace `gen'  = 30     if `touse' & missing(`gen') & inrange(`sicnum', 2900, 2912)
replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

replace `gen'  = 30     if `touse' & missing(`gen') & inrange(`sicnum', 2990, 2999)
replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

replace `gen'  = 31     if `touse' & missing(`gen') & `sicnum'==4900
replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

replace `gen'  = 31     if `touse' & missing(`gen') & inrange(`sicnum', 4910, 4911)
replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

replace `gen'  = 31     if `touse' & missing(`gen') & inrange(`sicnum', 4920, 4922)
replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

replace `gen'  = 31     if `touse' & missing(`gen') & `sicnum'==4923
replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

replace `gen'  = 31     if `touse' & missing(`gen') & inrange(`sicnum', 4924, 4925)
replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

replace `gen'  = 31     if `touse' & missing(`gen') & inrange(`sicnum', 4930, 4931)
replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

replace `gen'  = 31     if `touse' & missing(`gen') & `sicnum'==4932
replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

replace `gen'  = 31     if `touse' & missing(`gen') & `sicnum'==4939
replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

replace `gen'  = 31     if `touse' & missing(`gen') & inrange(`sicnum', 4940, 4942)
replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

replace `gen'  = 32     if `touse' & missing(`gen') & `sicnum'==4800
replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

replace `gen'  = 32     if `touse' & missing(`gen') & inrange(`sicnum', 4810, 4813)
replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

replace `gen'  = 32     if `touse' & missing(`gen') & inrange(`sicnum', 4820, 4822)
replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

replace `gen'  = 32     if `touse' & missing(`gen') & inrange(`sicnum', 4830, 4839)
replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

replace `gen'  = 32     if `touse' & missing(`gen') & inrange(`sicnum', 4840, 4841)
replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

replace `gen'  = 32     if `touse' & missing(`gen') & inrange(`sicnum', 4880, 4889)
replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

replace `gen'  = 32     if `touse' & missing(`gen') & `sicnum'==4890
replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

replace `gen'  = 32     if `touse' & missing(`gen') & `sicnum'==4891
replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

replace `gen'  = 32     if `touse' & missing(`gen') & `sicnum'==4892
replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

replace `gen'  = 32     if `touse' & missing(`gen') & `sicnum'==4899
replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7020, 7021)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7030, 7033)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7200
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7210, 7212)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7214
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7215, 7216)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7217
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7219
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7220, 7221)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7230, 7231)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7240, 7241)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7250, 7251)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7260, 7269)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7270, 7290)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7291
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7292, 7299)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7395
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7500
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7520, 7529)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7530, 7539)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7540, 7549)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7600
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7620
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7622
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7623
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7629
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7630, 7631)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7640, 7641)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7690, 7699)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 8100, 8199)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 8200, 8299)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 8300, 8399)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 8400, 8499)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 8600, 8699)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 8800, 8899)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7510, 7515)
replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 2750, 2759)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==3993
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7218
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7300
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7310, 7319)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7320, 7329)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7330, 7339)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7340, 7342)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7349
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7350, 7351)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7352
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7353
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7359
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7360, 7369)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7370, 7372)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7374
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7375
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7376
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7377
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7378
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7379
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7380
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7381, 7382)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7383
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7384
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7385
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7389, 7390)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7391
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7392
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7393
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7394
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7396
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7397
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7399
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7519
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==8700
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 8710, 8713)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 8720, 8721)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 8730, 8734)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 8740, 8748)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 8900, 8910)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==8911
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 8920, 8999)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 4220, 4229)
replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

replace `gen'  = 35     if `touse' & missing(`gen') & inrange(`sicnum', 3570, 3579)
replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3680
replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3681
replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3682
replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3683
replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3684
replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3685
replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3686
replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3687
replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3688
replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3689
replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3695
replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==7373
replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3622
replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3661
replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3662
replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3663
replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3664
replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3665
replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3666
replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3669
replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

replace `gen'  = 36     if `touse' & missing(`gen') & inrange(`sicnum', 3670, 3679)
replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3810
replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3812
replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3811
replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3820
replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3821
replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3822
replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3823
replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3824
replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3825
replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3826
replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3827
replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3829
replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

replace `gen'  = 37     if `touse' & missing(`gen') & inrange(`sicnum', 3830, 3839)
replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

replace `gen'  = 38     if `touse' & missing(`gen') & inrange(`sicnum', 2520, 2549)
replace `desc' = "Paper" if `touse' & `gen'==38 & `desc'==""

replace `gen'  = 38     if `touse' & missing(`gen') & inrange(`sicnum', 2600, 2639)
replace `desc' = "Paper" if `touse' & `gen'==38 & `desc'==""

replace `gen'  = 38     if `touse' & missing(`gen') & inrange(`sicnum', 2670, 2699)
replace `desc' = "Paper" if `touse' & `gen'==38 & `desc'==""

replace `gen'  = 38     if `touse' & missing(`gen') & inrange(`sicnum', 2760, 2761)
replace `desc' = "Paper" if `touse' & `gen'==38 & `desc'==""

replace `gen'  = 38     if `touse' & missing(`gen') & inrange(`sicnum', 3950, 3955)
replace `desc' = "Paper" if `touse' & `gen'==38 & `desc'==""

replace `gen'  = 39     if `touse' & missing(`gen') & inrange(`sicnum', 2440, 2449)
replace `desc' = "Boxes" if `touse' & `gen'==39 & `desc'==""

replace `gen'  = 39     if `touse' & missing(`gen') & inrange(`sicnum', 2640, 2659)
replace `desc' = "Boxes" if `touse' & `gen'==39 & `desc'==""

replace `gen'  = 39     if `touse' & missing(`gen') & inrange(`sicnum', 3220, 3221)
replace `desc' = "Boxes" if `touse' & `gen'==39 & `desc'==""

replace `gen'  = 39     if `touse' & missing(`gen') & inrange(`sicnum', 3410, 3412)
replace `desc' = "Boxes" if `touse' & `gen'==39 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4000, 4013)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4040, 4049)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4100
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4110, 4119)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4120, 4121)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4130, 4131)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4140, 4142)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4150, 4151)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4170, 4173)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4190, 4199)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4200
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4210, 4219)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4230, 4231)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4240, 4249)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4400, 4499)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4500, 4599)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4600, 4699)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4700
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4710, 4712)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4720, 4729)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4730, 4739)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4740, 4749)
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4780
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4782
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4783
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4784
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4785
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4789
replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5000
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5010, 5015)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5020, 5023)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5030, 5039)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5040, 5042)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5043
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5044
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5045
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5046
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5047
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5048
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5049
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5050, 5059)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5060
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5063
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5064
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5065
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5070, 5078)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5080
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5081
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5082
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5083
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5084
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5085
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5086, 5087)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5088
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5090
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5091, 5092)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5093
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5094
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5099
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5100
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5110, 5113)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5120, 5122)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5130, 5139)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5140, 5149)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5150, 5159)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5160, 5169)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5170, 5172)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5180, 5182)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5190, 5199)
replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5200
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5210, 5219)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5220, 5229)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5230, 5231)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5250, 5251)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5260, 5261)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5270, 5271)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5300
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5310, 5311)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5320
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5330, 5331)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5334
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5340, 5349)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5390, 5399)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5400
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5410, 5411)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5412
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5420, 5429)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5430, 5439)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5440, 5449)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5450, 5459)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5460, 5469)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5490, 5499)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5500
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5510, 5529)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5530, 5539)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5540, 5549)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5550, 5559)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5560, 5569)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5570, 5579)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5590, 5599)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5600, 5699)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5700
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5710, 5719)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5720, 5722)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5730, 5733)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5734
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5735
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5736
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5750, 5799)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5900
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5910, 5912)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5920, 5929)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5930, 5932)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5940
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5941
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5942
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5943
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5944
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5945
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5946
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5947
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5948
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5949
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5950, 5959)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5960, 5969)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5970, 5979)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5980, 5989)
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5990
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5992
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5993
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5994
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5995
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5999
replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

replace `gen'  = 43     if `touse' & missing(`gen') & inrange(`sicnum', 5800, 5819)
replace `desc' = "Meals" if `touse' & `gen'==43 & `desc'==""

replace `gen'  = 43     if `touse' & missing(`gen') & inrange(`sicnum', 5820, 5829)
replace `desc' = "Meals" if `touse' & `gen'==43 & `desc'==""

replace `gen'  = 43     if `touse' & missing(`gen') & inrange(`sicnum', 5890, 5899)
replace `desc' = "Meals" if `touse' & `gen'==43 & `desc'==""

replace `gen'  = 43     if `touse' & missing(`gen') & `sicnum'==7000
replace `desc' = "Meals" if `touse' & `gen'==43 & `desc'==""

replace `gen'  = 43     if `touse' & missing(`gen') & inrange(`sicnum', 7010, 7019)
replace `desc' = "Meals" if `touse' & `gen'==43 & `desc'==""

replace `gen'  = 43     if `touse' & missing(`gen') & inrange(`sicnum', 7040, 7049)
replace `desc' = "Meals" if `touse' & `gen'==43 & `desc'==""

replace `gen'  = 43     if `touse' & missing(`gen') & `sicnum'==7213
replace `desc' = "Meals" if `touse' & `gen'==43 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6000
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6010, 6019)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6020
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6021
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6022
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6023, 6024)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6025
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6026
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6027
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6028, 6029)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6030, 6036)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6040, 6059)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6060, 6062)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6080, 6082)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6090, 6099)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6100
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6110, 6111)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6112, 6113)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6120, 6129)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6130, 6139)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6140, 6149)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6150, 6159)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6160, 6169)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6170, 6179)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6190, 6199)
replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

replace `gen'  = 45     if `touse' & missing(`gen') & `sicnum'==6300
replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6310, 6319)
replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6320, 6329)
replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6330, 6331)
replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6350, 6351)
replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6360, 6361)
replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6370, 6379)
replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6390, 6399)
replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6400, 6411)
replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

replace `gen'  = 46     if `touse' & missing(`gen') & `sicnum'==6500
replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

replace `gen'  = 46     if `touse' & missing(`gen') & `sicnum'==6510
replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

replace `gen'  = 46     if `touse' & missing(`gen') & `sicnum'==6512
replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

replace `gen'  = 46     if `touse' & missing(`gen') & `sicnum'==6513
replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

replace `gen'  = 46     if `touse' & missing(`gen') & `sicnum'==6514
replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

replace `gen'  = 46     if `touse' & missing(`gen') & `sicnum'==6515
replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

replace `gen'  = 46     if `touse' & missing(`gen') & inrange(`sicnum', 6517, 6519)
replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

replace `gen'  = 46     if `touse' & missing(`gen') & inrange(`sicnum', 6520, 6529)
replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

replace `gen'  = 46     if `touse' & missing(`gen') & inrange(`sicnum', 6530, 6531)
replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

replace `gen'  = 46     if `touse' & missing(`gen') & `sicnum'==6532
replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

replace `gen'  = 46     if `touse' & missing(`gen') & inrange(`sicnum', 6540, 6541)
replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

replace `gen'  = 46     if `touse' & missing(`gen') & inrange(`sicnum', 6550, 6553)
replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

replace `gen'  = 46     if `touse' & missing(`gen') & inrange(`sicnum', 6590, 6599)
replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

replace `gen'  = 46     if `touse' & missing(`gen') & inrange(`sicnum', 6610, 6611)
replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & inrange(`sicnum', 6200, 6299)
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6700
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & inrange(`sicnum', 6710, 6719)
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & inrange(`sicnum', 6720, 6722)
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6723
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6724
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6725
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6726
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & inrange(`sicnum', 6730, 6733)
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & inrange(`sicnum', 6740, 6779)
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & inrange(`sicnum', 6790, 6791)
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6792
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6793
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6794
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6795
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6798
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6799
replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

replace `gen'  = 48     if `touse' & missing(`gen') & inrange(`sicnum', 4950, 4959)
replace `desc' = "Other" if `touse' & `gen'==48 & `desc'==""

replace `gen'  = 48     if `touse' & missing(`gen') & inrange(`sicnum', 4960, 4961)
replace `desc' = "Other" if `touse' & `gen'==48 & `desc'==""

replace `gen'  = 48     if `touse' & missing(`gen') & inrange(`sicnum', 4970, 4971)
replace `desc' = "Other" if `touse' & `gen'==48 & `desc'==""

replace `gen'  = 48     if `touse' & missing(`gen') & inrange(`sicnum', 4990, 4991)
replace `desc' = "Other" if `touse' & `gen'==48 & `desc'==""

* ---- END GENERATED RULES ----
