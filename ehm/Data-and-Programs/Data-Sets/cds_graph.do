use "monthly cds.dta", clear
format month %tm

twoway line spread5y month, ///
    by(shortname, cols(5) ixaxes note("")) ///
    xlabel(`=ym(2002,1)'(24)`=ym(2014,1)', ///
        format(%tmCCYY) angle(45) labsize(vsmall)) ///
    ytitle("5-year CDS spread") ///
    xtitle("Month")
graph export "monthly_cds.png", replace