program define ffi48_macro
    version 13.0
    /*
      Create Fama-French 48 industry classification from SIC.

      Usage:
        ffi48 sicvar, gen(newvar) desc(newstrvar)

      Options:
        gen()   : name of numeric output variable for FF48 code
        desc()  : name of string output variable for 5-char description
        replace : allow overwriting existing gen()/desc() variables
    */

    syntax varname [if] [in], GEN(name) DESC(name) [REPLACE]

    marksample touse, novarlist

    // Validate/handle existing vars
    capture confirm new variable `gen'
    if (_rc) {
        if ("`replace'" == "") {
            di as err "Variable `gen' already exists. Specify option replace to overwrite."
            exit 110
        }
        else {
            capture drop `gen'
        }
    }

    capture confirm new variable `desc'
    if (_rc) {
        if ("`replace'" == "") {
            di as err "Variable `desc' already exists. Specify option replace to overwrite."
            exit 110
        }
        else {
            capture drop `desc'
        }
    }

    // Coerce SIC to numeric safely (handles string SIC like "0100")
    tempvar sicnum
    gen double `sicnum' = .

    capture confirm numeric variable `varlist'
    if !_rc {
        replace `sicnum' = `varlist' if `touse'
    }
    else {
        replace `sicnum' = real(trim(`varlist')) if `touse'
    }



    // Initialize outputs
    gen int `gen' = . if `touse'
    gen str5 `desc' = "" if `touse'

    /*
      IMPORTANT:
      SAS uses "else if" so the FIRST matching condition wins.
      In Stata, you must prevent later rules from overwriting earlier ones.
      We do that by guarding every rule with: if missing(`gen') & <condition>
    */

    * ---- BEGIN GENERATED RULES ----
* Assumes tempvar sicnum exists and output vars `gen' and `desc' exist.
* Guarded with missing(`gen') to mimic SAS else-if (first match wins).

qui replace `gen'  = 1     if `touse' & missing(`gen') & inrange(`sicnum', 100, 199)
qui replace `desc' = "Agric" if `touse' & `gen'==1 & `desc'==""

qui replace `gen'  = 1     if `touse' & missing(`gen') & inrange(`sicnum', 200, 299)
qui replace `desc' = "Agric" if `touse' & `gen'==1 & `desc'==""

qui replace `gen'  = 1     if `touse' & missing(`gen') & inrange(`sicnum', 700, 799)
qui replace `desc' = "Agric" if `touse' & `gen'==1 & `desc'==""

qui replace `gen'  = 1     if `touse' & missing(`gen') & inrange(`sicnum', 910, 919)
qui replace `desc' = "Agric" if `touse' & `gen'==1 & `desc'==""

qui replace `gen'  = 1     if `touse' & missing(`gen') & `sicnum'==2048
qui replace `desc' = "Agric" if `touse' & `gen'==1 & `desc'==""

qui replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2000, 2009)
qui replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

qui replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2010, 2019)
qui replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

qui replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2020, 2029)
qui replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

qui replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2030, 2039)
qui replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

qui replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2040, 2046)
qui replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

qui replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2050, 2059)
qui replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

qui replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2060, 2063)
qui replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

qui replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2070, 2079)
qui replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

qui replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2090, 2092)
qui replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

qui replace `gen'  = 2     if `touse' & missing(`gen') & `sicnum'==2095
qui replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

qui replace `gen'  = 2     if `touse' & missing(`gen') & inrange(`sicnum', 2098, 2099)
qui replace `desc' = "Food" if `touse' & `gen'==2 & `desc'==""

qui replace `gen'  = 3     if `touse' & missing(`gen') & inrange(`sicnum', 2064, 2068)
qui replace `desc' = "Soda" if `touse' & `gen'==3 & `desc'==""

qui replace `gen'  = 3     if `touse' & missing(`gen') & `sicnum'==2086
qui replace `desc' = "Soda" if `touse' & `gen'==3 & `desc'==""

qui replace `gen'  = 3     if `touse' & missing(`gen') & `sicnum'==2087
qui replace `desc' = "Soda" if `touse' & `gen'==3 & `desc'==""

qui replace `gen'  = 3     if `touse' & missing(`gen') & `sicnum'==2096
qui replace `desc' = "Soda" if `touse' & `gen'==3 & `desc'==""

qui replace `gen'  = 3     if `touse' & missing(`gen') & `sicnum'==2097
qui replace `desc' = "Soda" if `touse' & `gen'==3 & `desc'==""

qui replace `gen'  = 4     if `touse' & missing(`gen') & `sicnum'==2080
qui replace `desc' = "Beer" if `touse' & `gen'==4 & `desc'==""

qui replace `gen'  = 4     if `touse' & missing(`gen') & `sicnum'==2082
qui replace `desc' = "Beer" if `touse' & `gen'==4 & `desc'==""

qui replace `gen'  = 4     if `touse' & missing(`gen') & `sicnum'==2083
qui replace `desc' = "Beer" if `touse' & `gen'==4 & `desc'==""

qui replace `gen'  = 4     if `touse' & missing(`gen') & `sicnum'==2084
qui replace `desc' = "Beer" if `touse' & `gen'==4 & `desc'==""

qui replace `gen'  = 4     if `touse' & missing(`gen') & `sicnum'==2085
qui replace `desc' = "Beer" if `touse' & `gen'==4 & `desc'==""

qui replace `gen'  = 5     if `touse' & missing(`gen') & inrange(`sicnum', 2100, 2199)
qui replace `desc' = "Smoke" if `touse' & `gen'==5 & `desc'==""

qui replace `gen'  = 6     if `touse' & missing(`gen') & inrange(`sicnum', 920, 999)
qui replace `desc' = "Toys" if `touse' & `gen'==6 & `desc'==""

qui replace `gen'  = 6     if `touse' & missing(`gen') & inrange(`sicnum', 3650, 3651)
qui replace `desc' = "Toys" if `touse' & `gen'==6 & `desc'==""

qui replace `gen'  = 6     if `touse' & missing(`gen') & `sicnum'==3652
qui replace `desc' = "Toys" if `touse' & `gen'==6 & `desc'==""

qui replace `gen'  = 6     if `touse' & missing(`gen') & `sicnum'==3732
qui replace `desc' = "Toys" if `touse' & `gen'==6 & `desc'==""

qui replace `gen'  = 6     if `touse' & missing(`gen') & inrange(`sicnum', 3930, 3931)
qui replace `desc' = "Toys" if `touse' & `gen'==6 & `desc'==""

qui replace `gen'  = 6     if `touse' & missing(`gen') & inrange(`sicnum', 3940, 3949)
qui replace `desc' = "Toys" if `touse' & `gen'==6 & `desc'==""

qui replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7800, 7829)
qui replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

qui replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7830, 7833)
qui replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

qui replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7840, 7841)
qui replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

qui replace `gen'  = 7     if `touse' & missing(`gen') & `sicnum'==7900
qui replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

qui replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7910, 7911)
qui replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

qui replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7920, 7929)
qui replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

qui replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7930, 7933)
qui replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

qui replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7940, 7949)
qui replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

qui replace `gen'  = 7     if `touse' & missing(`gen') & `sicnum'==7980
qui replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

qui replace `gen'  = 7     if `touse' & missing(`gen') & inrange(`sicnum', 7990, 7999)
qui replace `desc' = "Fun" if `touse' & `gen'==7 & `desc'==""

qui replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2700, 2709)
qui replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

qui replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2710, 2719)
qui replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

qui replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2720, 2729)
qui replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

qui replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2730, 2739)
qui replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

qui replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2740, 2749)
qui replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

qui replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2770, 2771)
qui replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

qui replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2780, 2789)
qui replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

qui replace `gen'  = 8     if `touse' & missing(`gen') & inrange(`sicnum', 2790, 2799)
qui replace `desc' = "Books" if `touse' & `gen'==8 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==2047
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 2391, 2392)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 2510, 2519)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 2590, 2599)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 2840, 2843)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==2844
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3160, 3161)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3170, 3171)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3172
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3190, 3199)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3229
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3260
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3262, 3263)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3269
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3230, 3231)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3630, 3639)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3750, 3751)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3800
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3860, 3861)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3870, 3873)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3910, 3911)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3914
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3915
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & inrange(`sicnum', 3960, 3962)
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3991
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 9     if `touse' & missing(`gen') & `sicnum'==3995
qui replace `desc' = "Hshld" if `touse' & `gen'==9 & `desc'==""

qui replace `gen'  = 10     if `touse' & missing(`gen') & inrange(`sicnum', 2300, 2390)
qui replace `desc' = "Clths" if `touse' & `gen'==10 & `desc'==""

qui replace `gen'  = 10     if `touse' & missing(`gen') & inrange(`sicnum', 3020, 3021)
qui replace `desc' = "Clths" if `touse' & `gen'==10 & `desc'==""

qui replace `gen'  = 10     if `touse' & missing(`gen') & inrange(`sicnum', 3100, 3111)
qui replace `desc' = "Clths" if `touse' & `gen'==10 & `desc'==""

qui replace `gen'  = 10     if `touse' & missing(`gen') & inrange(`sicnum', 3130, 3131)
qui replace `desc' = "Clths" if `touse' & `gen'==10 & `desc'==""

qui replace `gen'  = 10     if `touse' & missing(`gen') & inrange(`sicnum', 3140, 3149)
qui replace `desc' = "Clths" if `touse' & `gen'==10 & `desc'==""

qui replace `gen'  = 10     if `touse' & missing(`gen') & inrange(`sicnum', 3150, 3151)
qui replace `desc' = "Clths" if `touse' & `gen'==10 & `desc'==""

qui replace `gen'  = 10     if `touse' & missing(`gen') & inrange(`sicnum', 3963, 3965)
qui replace `desc' = "Clths" if `touse' & `gen'==10 & `desc'==""

qui replace `gen'  = 11     if `touse' & missing(`gen') & inrange(`sicnum', 8000, 8099)
qui replace `desc' = "Hlth" if `touse' & `gen'==11 & `desc'==""

qui replace `gen'  = 12     if `touse' & missing(`gen') & `sicnum'==3693
qui replace `desc' = "MedEq" if `touse' & `gen'==12 & `desc'==""

qui replace `gen'  = 12     if `touse' & missing(`gen') & inrange(`sicnum', 3840, 3849)
qui replace `desc' = "MedEq" if `touse' & `gen'==12 & `desc'==""

qui replace `gen'  = 12     if `touse' & missing(`gen') & inrange(`sicnum', 3850, 3851)
qui replace `desc' = "MedEq" if `touse' & `gen'==12 & `desc'==""

qui replace `gen'  = 13     if `touse' & missing(`gen') & `sicnum'==2830
qui replace `desc' = "Drugs" if `touse' & `gen'==13 & `desc'==""

qui replace `gen'  = 13     if `touse' & missing(`gen') & `sicnum'==2831
qui replace `desc' = "Drugs" if `touse' & `gen'==13 & `desc'==""

qui replace `gen'  = 13     if `touse' & missing(`gen') & `sicnum'==2833
qui replace `desc' = "Drugs" if `touse' & `gen'==13 & `desc'==""

qui replace `gen'  = 13     if `touse' & missing(`gen') & `sicnum'==2834
qui replace `desc' = "Drugs" if `touse' & `gen'==13 & `desc'==""

qui replace `gen'  = 13     if `touse' & missing(`gen') & `sicnum'==2835
qui replace `desc' = "Drugs" if `touse' & `gen'==13 & `desc'==""

qui replace `gen'  = 13     if `touse' & missing(`gen') & `sicnum'==2836
qui replace `desc' = "Drugs" if `touse' & `gen'==13 & `desc'==""

qui replace `gen'  = 14     if `touse' & missing(`gen') & inrange(`sicnum', 2800, 2809)
qui replace `desc' = "Chems" if `touse' & `gen'==14 & `desc'==""

qui replace `gen'  = 14     if `touse' & missing(`gen') & inrange(`sicnum', 2810, 2819)
qui replace `desc' = "Chems" if `touse' & `gen'==14 & `desc'==""

qui replace `gen'  = 14     if `touse' & missing(`gen') & inrange(`sicnum', 2820, 2829)
qui replace `desc' = "Chems" if `touse' & `gen'==14 & `desc'==""

qui replace `gen'  = 14     if `touse' & missing(`gen') & inrange(`sicnum', 2850, 2859)
qui replace `desc' = "Chems" if `touse' & `gen'==14 & `desc'==""

qui replace `gen'  = 14     if `touse' & missing(`gen') & inrange(`sicnum', 2860, 2869)
qui replace `desc' = "Chems" if `touse' & `gen'==14 & `desc'==""

qui replace `gen'  = 14     if `touse' & missing(`gen') & inrange(`sicnum', 2870, 2879)
qui replace `desc' = "Chems" if `touse' & `gen'==14 & `desc'==""

qui replace `gen'  = 14     if `touse' & missing(`gen') & inrange(`sicnum', 2890, 2899)
qui replace `desc' = "Chems" if `touse' & `gen'==14 & `desc'==""

qui replace `gen'  = 15     if `touse' & missing(`gen') & `sicnum'==3031
qui replace `desc' = "Rubbr" if `touse' & `gen'==15 & `desc'==""

qui replace `gen'  = 15     if `touse' & missing(`gen') & `sicnum'==3041
qui replace `desc' = "Rubbr" if `touse' & `gen'==15 & `desc'==""

qui replace `gen'  = 15     if `touse' & missing(`gen') & inrange(`sicnum', 3050, 3053)
qui replace `desc' = "Rubbr" if `touse' & `gen'==15 & `desc'==""

qui replace `gen'  = 15     if `touse' & missing(`gen') & inrange(`sicnum', 3060, 3069)
qui replace `desc' = "Rubbr" if `touse' & `gen'==15 & `desc'==""

qui replace `gen'  = 15     if `touse' & missing(`gen') & inrange(`sicnum', 3070, 3079)
qui replace `desc' = "Rubbr" if `touse' & `gen'==15 & `desc'==""

qui replace `gen'  = 15     if `touse' & missing(`gen') & inrange(`sicnum', 3080, 3089)
qui replace `desc' = "Rubbr" if `touse' & `gen'==15 & `desc'==""

qui replace `gen'  = 15     if `touse' & missing(`gen') & inrange(`sicnum', 3090, 3099)
qui replace `desc' = "Rubbr" if `touse' & `gen'==15 & `desc'==""

qui replace `gen'  = 16     if `touse' & missing(`gen') & inrange(`sicnum', 2200, 2269)
qui replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

qui replace `gen'  = 16     if `touse' & missing(`gen') & inrange(`sicnum', 2270, 2279)
qui replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

qui replace `gen'  = 16     if `touse' & missing(`gen') & inrange(`sicnum', 2280, 2284)
qui replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

qui replace `gen'  = 16     if `touse' & missing(`gen') & inrange(`sicnum', 2290, 2295)
qui replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

qui replace `gen'  = 16     if `touse' & missing(`gen') & `sicnum'==2297
qui replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

qui replace `gen'  = 16     if `touse' & missing(`gen') & `sicnum'==2298
qui replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

qui replace `gen'  = 16     if `touse' & missing(`gen') & `sicnum'==2299
qui replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

qui replace `gen'  = 16     if `touse' & missing(`gen') & inrange(`sicnum', 2393, 2395)
qui replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

qui replace `gen'  = 16     if `touse' & missing(`gen') & inrange(`sicnum', 2397, 2399)
qui replace `desc' = "Txtls" if `touse' & `gen'==16 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 800, 899)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 2400, 2439)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 2450, 2459)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 2490, 2499)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 2660, 2661)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 2950, 2952)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3200
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3210, 3211)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3240, 3241)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3250, 3259)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3261
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3264
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3270, 3275)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3280, 3281)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3290, 3293)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3295, 3299)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3420, 3429)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3430, 3433)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3440, 3441)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3442
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3446
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3448
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3449
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3450, 3451)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3452
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & inrange(`sicnum', 3490, 3499)
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 17     if `touse' & missing(`gen') & `sicnum'==3996
qui replace `desc' = "BldMt" if `touse' & `gen'==17 & `desc'==""

qui replace `gen'  = 18     if `touse' & missing(`gen') & inrange(`sicnum', 1500, 1511)
qui replace `desc' = "Cnstr" if `touse' & `gen'==18 & `desc'==""

qui replace `gen'  = 18     if `touse' & missing(`gen') & inrange(`sicnum', 1520, 1529)
qui replace `desc' = "Cnstr" if `touse' & `gen'==18 & `desc'==""

qui replace `gen'  = 18     if `touse' & missing(`gen') & inrange(`sicnum', 1530, 1539)
qui replace `desc' = "Cnstr" if `touse' & `gen'==18 & `desc'==""

qui replace `gen'  = 18     if `touse' & missing(`gen') & inrange(`sicnum', 1540, 1549)
qui replace `desc' = "Cnstr" if `touse' & `gen'==18 & `desc'==""

qui replace `gen'  = 18     if `touse' & missing(`gen') & inrange(`sicnum', 1600, 1699)
qui replace `desc' = "Cnstr" if `touse' & `gen'==18 & `desc'==""

qui replace `gen'  = 18     if `touse' & missing(`gen') & inrange(`sicnum', 1700, 1799)
qui replace `desc' = "Cnstr" if `touse' & `gen'==18 & `desc'==""

qui replace `gen'  = 19     if `touse' & missing(`gen') & `sicnum'==3300
qui replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

qui replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3310, 3317)
qui replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

qui replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3320, 3325)
qui replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

qui replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3330, 3339)
qui replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

qui replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3340, 3341)
qui replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

qui replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3350, 3357)
qui replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

qui replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3360, 3369)
qui replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

qui replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3370, 3379)
qui replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

qui replace `gen'  = 19     if `touse' & missing(`gen') & inrange(`sicnum', 3390, 3399)
qui replace `desc' = "Steel" if `touse' & `gen'==19 & `desc'==""

qui replace `gen'  = 20     if `touse' & missing(`gen') & `sicnum'==3400
qui replace `desc' = "FabPr" if `touse' & `gen'==20 & `desc'==""

qui replace `gen'  = 20     if `touse' & missing(`gen') & `sicnum'==3443
qui replace `desc' = "FabPr" if `touse' & `gen'==20 & `desc'==""

qui replace `gen'  = 20     if `touse' & missing(`gen') & `sicnum'==3444
qui replace `desc' = "FabPr" if `touse' & `gen'==20 & `desc'==""

qui replace `gen'  = 20     if `touse' & missing(`gen') & inrange(`sicnum', 3460, 3469)
qui replace `desc' = "FabPr" if `touse' & `gen'==20 & `desc'==""

qui replace `gen'  = 20     if `touse' & missing(`gen') & inrange(`sicnum', 3470, 3479)
qui replace `desc' = "FabPr" if `touse' & `gen'==20 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & inrange(`sicnum', 3510, 3519)
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & inrange(`sicnum', 3520, 3529)
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3530
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3531
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3532
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3533
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3534
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3535
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3536
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3538
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & inrange(`sicnum', 3540, 3549)
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & inrange(`sicnum', 3550, 3559)
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & inrange(`sicnum', 3560, 3569)
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3580
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3581
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3582
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3585
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3586
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & `sicnum'==3589
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 21     if `touse' & missing(`gen') & inrange(`sicnum', 3590, 3599)
qui replace `desc' = "Mach" if `touse' & `gen'==21 & `desc'==""

qui replace `gen'  = 22     if `touse' & missing(`gen') & `sicnum'==3600
qui replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

qui replace `gen'  = 22     if `touse' & missing(`gen') & inrange(`sicnum', 3610, 3613)
qui replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

qui replace `gen'  = 22     if `touse' & missing(`gen') & inrange(`sicnum', 3620, 3621)
qui replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

qui replace `gen'  = 22     if `touse' & missing(`gen') & inrange(`sicnum', 3623, 3629)
qui replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

qui replace `gen'  = 22     if `touse' & missing(`gen') & inrange(`sicnum', 3640, 3644)
qui replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

qui replace `gen'  = 22     if `touse' & missing(`gen') & `sicnum'==3645
qui replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

qui replace `gen'  = 22     if `touse' & missing(`gen') & `sicnum'==3646
qui replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

qui replace `gen'  = 22     if `touse' & missing(`gen') & inrange(`sicnum', 3648, 3649)
qui replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

qui replace `gen'  = 22     if `touse' & missing(`gen') & `sicnum'==3660
qui replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

qui replace `gen'  = 22     if `touse' & missing(`gen') & `sicnum'==3690
qui replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

qui replace `gen'  = 22     if `touse' & missing(`gen') & inrange(`sicnum', 3691, 3692)
qui replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

qui replace `gen'  = 22     if `touse' & missing(`gen') & `sicnum'==3699
qui replace `desc' = "ElcEq" if `touse' & `gen'==22 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==2296
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==2396
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & inrange(`sicnum', 3010, 3011)
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3537
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3647
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3694
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3700
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3710
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3711
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3713
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3714
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3715
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3716
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3792
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & inrange(`sicnum', 3790, 3791)
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 23     if `touse' & missing(`gen') & `sicnum'==3799
qui replace `desc' = "Autos" if `touse' & `gen'==23 & `desc'==""

qui replace `gen'  = 24     if `touse' & missing(`gen') & `sicnum'==3720
qui replace `desc' = "Aero" if `touse' & `gen'==24 & `desc'==""

qui replace `gen'  = 24     if `touse' & missing(`gen') & `sicnum'==3721
qui replace `desc' = "Aero" if `touse' & `gen'==24 & `desc'==""

qui replace `gen'  = 24     if `touse' & missing(`gen') & inrange(`sicnum', 3723, 3724)
qui replace `desc' = "Aero" if `touse' & `gen'==24 & `desc'==""

qui replace `gen'  = 24     if `touse' & missing(`gen') & `sicnum'==3725
qui replace `desc' = "Aero" if `touse' & `gen'==24 & `desc'==""

qui replace `gen'  = 24     if `touse' & missing(`gen') & inrange(`sicnum', 3728, 3729)
qui replace `desc' = "Aero" if `touse' & `gen'==24 & `desc'==""

qui replace `gen'  = 25     if `touse' & missing(`gen') & inrange(`sicnum', 3730, 3731)
qui replace `desc' = "Ships" if `touse' & `gen'==25 & `desc'==""

qui replace `gen'  = 25     if `touse' & missing(`gen') & inrange(`sicnum', 3740, 3743)
qui replace `desc' = "Ships" if `touse' & `gen'==25 & `desc'==""

qui replace `gen'  = 26     if `touse' & missing(`gen') & inrange(`sicnum', 3760, 3769)
qui replace `desc' = "Guns" if `touse' & `gen'==26 & `desc'==""

qui replace `gen'  = 26     if `touse' & missing(`gen') & `sicnum'==3795
qui replace `desc' = "Guns" if `touse' & `gen'==26 & `desc'==""

qui replace `gen'  = 26     if `touse' & missing(`gen') & inrange(`sicnum', 3480, 3489)
qui replace `desc' = "Guns" if `touse' & `gen'==26 & `desc'==""

qui replace `gen'  = 27     if `touse' & missing(`gen') & inrange(`sicnum', 1040, 1049)
qui replace `desc' = "Gold" if `touse' & `gen'==27 & `desc'==""

qui replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1000, 1009)
qui replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

qui replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1010, 1019)
qui replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

qui replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1020, 1029)
qui replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

qui replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1030, 1039)
qui replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

qui replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1050, 1059)
qui replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

qui replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1060, 1069)
qui replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

qui replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1070, 1079)
qui replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

qui replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1080, 1089)
qui replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

qui replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1090, 1099)
qui replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

qui replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1100, 1119)
qui replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

qui replace `gen'  = 28     if `touse' & missing(`gen') & inrange(`sicnum', 1400, 1499)
qui replace `desc' = "Mines" if `touse' & `gen'==28 & `desc'==""

qui replace `gen'  = 29     if `touse' & missing(`gen') & inrange(`sicnum', 1200, 1299)
qui replace `desc' = "Coal" if `touse' & `gen'==29 & `desc'==""

qui replace `gen'  = 30     if `touse' & missing(`gen') & `sicnum'==1300
qui replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

qui replace `gen'  = 30     if `touse' & missing(`gen') & inrange(`sicnum', 1310, 1319)
qui replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

qui replace `gen'  = 30     if `touse' & missing(`gen') & inrange(`sicnum', 1320, 1329)
qui replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

qui replace `gen'  = 30     if `touse' & missing(`gen') & inrange(`sicnum', 1330, 1339)
qui replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

qui replace `gen'  = 30     if `touse' & missing(`gen') & inrange(`sicnum', 1370, 1379)
qui replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

qui replace `gen'  = 30     if `touse' & missing(`gen') & `sicnum'==1380
qui replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

qui replace `gen'  = 30     if `touse' & missing(`gen') & `sicnum'==1381
qui replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

qui replace `gen'  = 30     if `touse' & missing(`gen') & `sicnum'==1382
qui replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

qui replace `gen'  = 30     if `touse' & missing(`gen') & `sicnum'==1389
qui replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

qui replace `gen'  = 30     if `touse' & missing(`gen') & inrange(`sicnum', 2900, 2912)
qui replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

qui replace `gen'  = 30     if `touse' & missing(`gen') & inrange(`sicnum', 2990, 2999)
qui replace `desc' = "Oil" if `touse' & `gen'==30 & `desc'==""

qui replace `gen'  = 31     if `touse' & missing(`gen') & `sicnum'==4900
qui replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

qui replace `gen'  = 31     if `touse' & missing(`gen') & inrange(`sicnum', 4910, 4911)
qui replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

qui replace `gen'  = 31     if `touse' & missing(`gen') & inrange(`sicnum', 4920, 4922)
qui replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

qui replace `gen'  = 31     if `touse' & missing(`gen') & `sicnum'==4923
qui replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

qui replace `gen'  = 31     if `touse' & missing(`gen') & inrange(`sicnum', 4924, 4925)
qui replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

qui replace `gen'  = 31     if `touse' & missing(`gen') & inrange(`sicnum', 4930, 4931)
qui replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

qui replace `gen'  = 31     if `touse' & missing(`gen') & `sicnum'==4932
qui replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

qui replace `gen'  = 31     if `touse' & missing(`gen') & `sicnum'==4939
qui replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

qui replace `gen'  = 31     if `touse' & missing(`gen') & inrange(`sicnum', 4940, 4942)
qui replace `desc' = "Util" if `touse' & `gen'==31 & `desc'==""

qui replace `gen'  = 32     if `touse' & missing(`gen') & `sicnum'==4800
qui replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

qui replace `gen'  = 32     if `touse' & missing(`gen') & inrange(`sicnum', 4810, 4813)
qui replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

qui replace `gen'  = 32     if `touse' & missing(`gen') & inrange(`sicnum', 4820, 4822)
qui replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

qui replace `gen'  = 32     if `touse' & missing(`gen') & inrange(`sicnum', 4830, 4839)
qui replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

qui replace `gen'  = 32     if `touse' & missing(`gen') & inrange(`sicnum', 4840, 4841)
qui replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

qui replace `gen'  = 32     if `touse' & missing(`gen') & inrange(`sicnum', 4880, 4889)
qui replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

qui replace `gen'  = 32     if `touse' & missing(`gen') & `sicnum'==4890
qui replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

qui replace `gen'  = 32     if `touse' & missing(`gen') & `sicnum'==4891
qui replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

qui replace `gen'  = 32     if `touse' & missing(`gen') & `sicnum'==4892
qui replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

qui replace `gen'  = 32     if `touse' & missing(`gen') & `sicnum'==4899
qui replace `desc' = "Telcm" if `touse' & `gen'==32 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7020, 7021)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7030, 7033)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7200
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7210, 7212)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7214
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7215, 7216)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7217
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7219
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7220, 7221)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7230, 7231)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7240, 7241)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7250, 7251)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7260, 7269)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7270, 7290)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7291
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7292, 7299)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7395
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7500
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7520, 7529)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7530, 7539)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7540, 7549)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7600
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7620
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7622
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7623
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & `sicnum'==7629
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7630, 7631)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7640, 7641)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7690, 7699)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 8100, 8199)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 8200, 8299)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 8300, 8399)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 8400, 8499)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 8600, 8699)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 8800, 8899)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 33     if `touse' & missing(`gen') & inrange(`sicnum', 7510, 7515)
qui replace `desc' = "PerSv" if `touse' & `gen'==33 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 2750, 2759)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==3993
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7218
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7300
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7310, 7319)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7320, 7329)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7330, 7339)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7340, 7342)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7349
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7350, 7351)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7352
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7353
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7359
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7360, 7369)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7370, 7372)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7374
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7375
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7376
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7377
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7378
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7379
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7380
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7381, 7382)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7383
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7384
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7385
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 7389, 7390)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7391
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7392
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7393
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7394
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7396
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7397
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7399
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==7519
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==8700
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 8710, 8713)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 8720, 8721)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 8730, 8734)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 8740, 8748)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 8900, 8910)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & `sicnum'==8911
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 8920, 8999)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 34     if `touse' & missing(`gen') & inrange(`sicnum', 4220, 4229)
qui replace `desc' = "BusSv" if `touse' & `gen'==34 & `desc'==""

qui replace `gen'  = 35     if `touse' & missing(`gen') & inrange(`sicnum', 3570, 3579)
qui replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

qui replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3680
qui replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

qui replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3681
qui replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

qui replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3682
qui replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

qui replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3683
qui replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

qui replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3684
qui replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

qui replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3685
qui replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

qui replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3686
qui replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

qui replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3687
qui replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

qui replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3688
qui replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

qui replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3689
qui replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

qui replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==3695
qui replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

qui replace `gen'  = 35     if `touse' & missing(`gen') & `sicnum'==7373
qui replace `desc' = "Comps" if `touse' & `gen'==35 & `desc'==""

qui replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3622
qui replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

qui replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3661
qui replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

qui replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3662
qui replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

qui replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3663
qui replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

qui replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3664
qui replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

qui replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3665
qui replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

qui replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3666
qui replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

qui replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3669
qui replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

qui replace `gen'  = 36     if `touse' & missing(`gen') & inrange(`sicnum', 3670, 3679)
qui replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

qui replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3810
qui replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

qui replace `gen'  = 36     if `touse' & missing(`gen') & `sicnum'==3812
qui replace `desc' = "Chips" if `touse' & `gen'==36 & `desc'==""

qui replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3811
qui replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

qui replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3820
qui replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

qui replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3821
qui replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

qui replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3822
qui replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

qui replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3823
qui replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

qui replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3824
qui replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

qui replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3825
qui replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

qui replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3826
qui replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

qui replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3827
qui replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

qui replace `gen'  = 37     if `touse' & missing(`gen') & `sicnum'==3829
qui replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

qui replace `gen'  = 37     if `touse' & missing(`gen') & inrange(`sicnum', 3830, 3839)
qui replace `desc' = "LabEq" if `touse' & `gen'==37 & `desc'==""

qui replace `gen'  = 38     if `touse' & missing(`gen') & inrange(`sicnum', 2520, 2549)
qui replace `desc' = "Paper" if `touse' & `gen'==38 & `desc'==""

qui replace `gen'  = 38     if `touse' & missing(`gen') & inrange(`sicnum', 2600, 2639)
qui replace `desc' = "Paper" if `touse' & `gen'==38 & `desc'==""

qui replace `gen'  = 38     if `touse' & missing(`gen') & inrange(`sicnum', 2670, 2699)
qui replace `desc' = "Paper" if `touse' & `gen'==38 & `desc'==""

qui replace `gen'  = 38     if `touse' & missing(`gen') & inrange(`sicnum', 2760, 2761)
qui replace `desc' = "Paper" if `touse' & `gen'==38 & `desc'==""

qui replace `gen'  = 38     if `touse' & missing(`gen') & inrange(`sicnum', 3950, 3955)
qui replace `desc' = "Paper" if `touse' & `gen'==38 & `desc'==""

qui replace `gen'  = 39     if `touse' & missing(`gen') & inrange(`sicnum', 2440, 2449)
qui replace `desc' = "Boxes" if `touse' & `gen'==39 & `desc'==""

qui replace `gen'  = 39     if `touse' & missing(`gen') & inrange(`sicnum', 2640, 2659)
qui replace `desc' = "Boxes" if `touse' & `gen'==39 & `desc'==""

qui replace `gen'  = 39     if `touse' & missing(`gen') & inrange(`sicnum', 3220, 3221)
qui replace `desc' = "Boxes" if `touse' & `gen'==39 & `desc'==""

qui replace `gen'  = 39     if `touse' & missing(`gen') & inrange(`sicnum', 3410, 3412)
qui replace `desc' = "Boxes" if `touse' & `gen'==39 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4000, 4013)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4040, 4049)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4100
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4110, 4119)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4120, 4121)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4130, 4131)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4140, 4142)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4150, 4151)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4170, 4173)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4190, 4199)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4200
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4210, 4219)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4230, 4231)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4240, 4249)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4400, 4499)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4500, 4599)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4600, 4699)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4700
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4710, 4712)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4720, 4729)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4730, 4739)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & inrange(`sicnum', 4740, 4749)
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4780
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4782
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4783
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4784
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4785
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 40     if `touse' & missing(`gen') & `sicnum'==4789
qui replace `desc' = "Trans" if `touse' & `gen'==40 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5000
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5010, 5015)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5020, 5023)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5030, 5039)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5040, 5042)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5043
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5044
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5045
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5046
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5047
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5048
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5049
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5050, 5059)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5060
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5063
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5064
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5065
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5070, 5078)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5080
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5081
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5082
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5083
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5084
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5085
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5086, 5087)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5088
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5090
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5091, 5092)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5093
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5094
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5099
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & `sicnum'==5100
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5110, 5113)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5120, 5122)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5130, 5139)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5140, 5149)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5150, 5159)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5160, 5169)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5170, 5172)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5180, 5182)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 41     if `touse' & missing(`gen') & inrange(`sicnum', 5190, 5199)
qui replace `desc' = "Whlsl" if `touse' & `gen'==41 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5200
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5210, 5219)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5220, 5229)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5230, 5231)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5250, 5251)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5260, 5261)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5270, 5271)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5300
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5310, 5311)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5320
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5330, 5331)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5334
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5340, 5349)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5390, 5399)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5400
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5410, 5411)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5412
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5420, 5429)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5430, 5439)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5440, 5449)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5450, 5459)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5460, 5469)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5490, 5499)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5500
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5510, 5529)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5530, 5539)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5540, 5549)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5550, 5559)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5560, 5569)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5570, 5579)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5590, 5599)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5600, 5699)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5700
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5710, 5719)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5720, 5722)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5730, 5733)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5734
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5735
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5736
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5750, 5799)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5900
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5910, 5912)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5920, 5929)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5930, 5932)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5940
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5941
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5942
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5943
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5944
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5945
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5946
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5947
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5948
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5949
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5950, 5959)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5960, 5969)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5970, 5979)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & inrange(`sicnum', 5980, 5989)
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5990
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5992
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5993
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5994
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5995
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 42     if `touse' & missing(`gen') & `sicnum'==5999
qui replace `desc' = "Rtail" if `touse' & `gen'==42 & `desc'==""

qui replace `gen'  = 43     if `touse' & missing(`gen') & inrange(`sicnum', 5800, 5819)
qui replace `desc' = "Meals" if `touse' & `gen'==43 & `desc'==""

qui replace `gen'  = 43     if `touse' & missing(`gen') & inrange(`sicnum', 5820, 5829)
qui replace `desc' = "Meals" if `touse' & `gen'==43 & `desc'==""

qui replace `gen'  = 43     if `touse' & missing(`gen') & inrange(`sicnum', 5890, 5899)
qui replace `desc' = "Meals" if `touse' & `gen'==43 & `desc'==""

qui replace `gen'  = 43     if `touse' & missing(`gen') & `sicnum'==7000
qui replace `desc' = "Meals" if `touse' & `gen'==43 & `desc'==""

qui replace `gen'  = 43     if `touse' & missing(`gen') & inrange(`sicnum', 7010, 7019)
qui replace `desc' = "Meals" if `touse' & `gen'==43 & `desc'==""

qui replace `gen'  = 43     if `touse' & missing(`gen') & inrange(`sicnum', 7040, 7049)
qui replace `desc' = "Meals" if `touse' & `gen'==43 & `desc'==""

qui replace `gen'  = 43     if `touse' & missing(`gen') & `sicnum'==7213
qui replace `desc' = "Meals" if `touse' & `gen'==43 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6000
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6010, 6019)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6020
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6021
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6022
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6023, 6024)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6025
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6026
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6027
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6028, 6029)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6030, 6036)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6040, 6059)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6060, 6062)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6080, 6082)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6090, 6099)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & `sicnum'==6100
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6110, 6111)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6112, 6113)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6120, 6129)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6130, 6139)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6140, 6149)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6150, 6159)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6160, 6169)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6170, 6179)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 44     if `touse' & missing(`gen') & inrange(`sicnum', 6190, 6199)
qui replace `desc' = "Banks" if `touse' & `gen'==44 & `desc'==""

qui replace `gen'  = 45     if `touse' & missing(`gen') & `sicnum'==6300
qui replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

qui replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6310, 6319)
qui replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

qui replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6320, 6329)
qui replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

qui replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6330, 6331)
qui replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

qui replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6350, 6351)
qui replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

qui replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6360, 6361)
qui replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

qui replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6370, 6379)
qui replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

qui replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6390, 6399)
qui replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

qui replace `gen'  = 45     if `touse' & missing(`gen') & inrange(`sicnum', 6400, 6411)
qui replace `desc' = "Insur" if `touse' & `gen'==45 & `desc'==""

qui replace `gen'  = 46     if `touse' & missing(`gen') & `sicnum'==6500
qui replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

qui replace `gen'  = 46     if `touse' & missing(`gen') & `sicnum'==6510
qui replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

qui replace `gen'  = 46     if `touse' & missing(`gen') & `sicnum'==6512
qui replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

qui replace `gen'  = 46     if `touse' & missing(`gen') & `sicnum'==6513
qui replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

qui replace `gen'  = 46     if `touse' & missing(`gen') & `sicnum'==6514
qui replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

qui replace `gen'  = 46     if `touse' & missing(`gen') & `sicnum'==6515
qui replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

qui replace `gen'  = 46     if `touse' & missing(`gen') & inrange(`sicnum', 6517, 6519)
qui replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

qui replace `gen'  = 46     if `touse' & missing(`gen') & inrange(`sicnum', 6520, 6529)
qui replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

qui replace `gen'  = 46     if `touse' & missing(`gen') & inrange(`sicnum', 6530, 6531)
qui replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

qui replace `gen'  = 46     if `touse' & missing(`gen') & `sicnum'==6532
qui replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

qui replace `gen'  = 46     if `touse' & missing(`gen') & inrange(`sicnum', 6540, 6541)
qui replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

qui replace `gen'  = 46     if `touse' & missing(`gen') & inrange(`sicnum', 6550, 6553)
qui replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

qui replace `gen'  = 46     if `touse' & missing(`gen') & inrange(`sicnum', 6590, 6599)
qui replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

qui replace `gen'  = 46     if `touse' & missing(`gen') & inrange(`sicnum', 6610, 6611)
qui replace `desc' = "RlEst" if `touse' & `gen'==46 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & inrange(`sicnum', 6200, 6299)
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6700
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & inrange(`sicnum', 6710, 6719)
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & inrange(`sicnum', 6720, 6722)
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6723
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6724
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6725
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6726
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & inrange(`sicnum', 6730, 6733)
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & inrange(`sicnum', 6740, 6779)
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & inrange(`sicnum', 6790, 6791)
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6792
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6793
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6794
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6795
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6798
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 47     if `touse' & missing(`gen') & `sicnum'==6799
qui replace `desc' = "Fin" if `touse' & `gen'==47 & `desc'==""

qui replace `gen'  = 48     if `touse' & missing(`gen') & inrange(`sicnum', 4950, 4959)
qui replace `desc' = "Other" if `touse' & `gen'==48 & `desc'==""

qui replace `gen'  = 48     if `touse' & missing(`gen') & inrange(`sicnum', 4960, 4961)
qui replace `desc' = "Other" if `touse' & `gen'==48 & `desc'==""

qui replace `gen'  = 48     if `touse' & missing(`gen') & inrange(`sicnum', 4970, 4971)
qui replace `desc' = "Other" if `touse' & `gen'==48 & `desc'==""

qui replace `gen'  = 48     if `touse' & missing(`gen') & inrange(`sicnum', 4990, 4991)
qui replace `desc' = "Other" if `touse' & `gen'==48 & `desc'==""

* ---- END GENERATED RULES ----

    // Variable labels (SAS "label")
    label variable `gen'  "Fama and French 48 Industries"
    label variable `desc' "Fama and French 48 Industries - Description"
end
