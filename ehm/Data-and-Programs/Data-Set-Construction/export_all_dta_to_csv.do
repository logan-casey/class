clear all
set more off

local data_dir "Data-and-Programs/Data-Sets"
local dta_files : dir "`data_dir'" files "*.dta"

foreach f of local dta_files {
    use "`data_dir'/`f'", clear
    local out = subinstr("`f'", ".dta", ".csv", .)
    export delimited using "`data_dir'/`out'", replace
    di as txt "Exported: `f' -> `out'"
}
