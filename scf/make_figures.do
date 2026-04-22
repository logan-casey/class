cd C:\Users\logan\Documents\GitHub\class\scf
use "rscfp2022.dta", clear

** FIGURE 1
preserve
tempfile curve_asset curve_fin curve_networth

* 1) Asset curve: weighted percentiles of asset, weighted mean asset by percentile
xtile p_asset = asset [pw=wgt], nq(100)
collapse (mean) mean_asset=asset [pw=wgt], by(p_asset)
rename p_asset pct
save `curve_asset'

* 2) Financial asset curve: weighted percentiles of fin, weighted mean fin by percentile
restore, preserve
xtile p_fin = fin [pw=wgt], nq(100)
collapse (mean) mean_fin=fin [pw=wgt], by(p_fin)
rename p_fin pct
save `curve_fin'

* 3) Net worth curve: weighted percentiles of networth, weighted mean networth by percentile
restore, preserve
xtile p_networth = networth [pw=wgt], nq(100)
collapse (mean) mean_networth=networth [pw=wgt], by(p_networth)
replace mean_networth = . if mean_networth < 0
rename p_networth pct
save `curve_networth'

* 4) Merge all curves on common percentile index for plotting
use `curve_asset', clear
merge 1:1 pct using `curve_fin', nogen
merge 1:1 pct using `curve_networth', nogen

	
* 5) Plot
twoway ///
    (line mean_asset    pct, lcolor(navy)   lwidth(medthick)) ///
    (line mean_fin      pct, lcolor(maroon) lwidth(medthick)) ///
    (line mean_networth pct, lcolor(forest_green) lwidth(medthick)), ///
	yscale(log) ylabel(100 "100"1000 "1k" 10000 "10k" 100000 "100k" 1000000 "1m" 10000000 "10m", angle(horizontal) labsize(small)) ///
    xtitle("Percentile") ///
    ytitle("Mean") ///
    legend(order(1 "Total assets" 2 "Financial assets" 3 "Net worth") position(6) ring(1) rows(1))
graph export "wealth_dist.png", replace width(2400)

restore

** FIGURE 2: Participation rates by asset class (x-axis = p_asset)
preserve
drop if missing(asset) | asset <= 0

* Treat missing/negative components as zero (asset-only construction)
foreach v in checking saving mma call cds savbnd vehic homeeq nnresre oresre equity bus bond cashli gbmutf obmutf trusts othma {
    gen `v'_pos = cond(missing(`v'), 0, cond(`v' < 0, 0, `v'))
}

gen safe_assets = checking_pos + saving_pos + mma_pos + call_pos + cds_pos + savbnd_pos
gen vehicles = vehic_pos
gen real_estate = homeeq_pos + nnresre_pos + oresre_pos
gen public_equity = equity_pos
gen private_business = bus_pos
gen bonds = bond_pos + cashli_pos + gbmutf_pos + obmutf_pos + trusts_pos + othma_pos

xtile p_asset = asset [pw=wgt], nq(100)

gen hold_safe_assets = safe_assets > 0 if !missing(safe_assets)
gen hold_vehicles = vehicles > 0 if !missing(vehicles)
gen hold_real_estate = real_estate > 0 if !missing(real_estate)
gen hold_public_equity = public_equity > 0 if !missing(public_equity)
gen hold_private_business = hbus > 0 if !missing(hbus)
gen hold_bonds = bonds > 0 if !missing(bonds)

collapse ///
    (mean) sh_hold_safe_assets=hold_safe_assets ///
    (mean) sh_hold_vehicles=hold_vehicles ///
    (mean) sh_hold_real_estate=hold_real_estate ///
    (mean) sh_hold_public_equity=hold_public_equity ///
    (mean) sh_hold_private_business=hold_private_business ///
    (mean) sh_hold_bonds=hold_bonds ///
    [pw=wgt], by(p_asset)

twoway ///
    (lpoly sh_hold_safe_assets p_asset, bwidth(5) lcolor(navy) lwidth(medthick)) ///
    (lpoly sh_hold_vehicles p_asset, bwidth(5) lcolor(maroon) lwidth(medthick)) ///
    (lpoly sh_hold_real_estate p_asset, bwidth(5) lcolor(green) lwidth(medthick)) ///
    (lpoly sh_hold_public_equity p_asset, bwidth(5) lcolor(orange) lwidth(medthick)) ///
    (lpoly sh_hold_private_business p_asset, bwidth(5) lcolor(brown) lwidth(medthick)) ///
    (lpoly sh_hold_bonds p_asset, bwidth(5) lcolor(gs6) lwidth(medthick)), ///
    xlabel(0(10)100, labsize(small)) ///
    ylabel(0(.2)1, format(%3.1f) angle(horizontal) labsize(small)) ///
    xtitle("Percentile of distribution of total assets") ///
    ytitle("Share with holdings") ///
    legend(order(1 "Safe assets" 2 "Vehicles" 3 "Real estate" 4 "Public equity" 5 "Private business" 6 "Bonds") position(6) ring(1) rows(1)) ///
    title("Participation Rates by Asset Class (Kernel Smoothed)")
graph export "participation_by_asset_percentile.png", replace width(2400)

twoway ///
    (line sh_hold_safe_assets p_asset, lcolor(navy) lwidth(medthick)) ///
    (line sh_hold_vehicles p_asset, lcolor(maroon) lwidth(medthick)) ///
    (line sh_hold_real_estate p_asset, lcolor(green) lwidth(medthick)) ///
    (line sh_hold_public_equity p_asset, lcolor(orange) lwidth(medthick)) ///
    (line sh_hold_private_business p_asset, lcolor(brown) lwidth(medthick)) ///
    (line sh_hold_bonds p_asset, lcolor(gs6) lwidth(medthick)), ///
    xlabel(0(10)100, labsize(small)) ///
    ylabel(0(.2)1, format(%3.1f) angle(horizontal) labsize(small)) ///
    xtitle("Percentile of distribution of total assets") ///
    ytitle("Share with holdings") ///
    legend(order(1 "Safe assets" 2 "Vehicles" 3 "Real estate" 4 "Public equity" 5 "Private business" 6 "Bonds") position(6) ring(1) rows(1)) ///
    title("Participation Rates by Asset Class (Raw)")
graph export "participation_by_asset_percentile_raw.png", replace width(2400)


restore

** FIGURE 3: Asset class shares in portfolios (x-axis = p_asset)
preserve
* Treat missing/negative components as zero (asset-only construction)
foreach v in checking saving mma call cds savbnd vehic homeeq nnresre oresre equity bus bond cashli gbmutf obmutf trusts othma {
    gen `v'_pos = cond(missing(`v'), 0, cond(`v' < 0, 0, `v'))
}

gen safe_assets = checking_pos + saving_pos + mma_pos + call_pos + cds_pos + savbnd_pos
gen vehicles = vehic_pos
gen real_estate = homeeq_pos + nnresre_pos + oresre_pos
gen public_equity = equity_pos
gen private_business = bus_pos
gen bonds = bond_pos + cashli_pos + gbmutf_pos + obmutf_pos + trusts_pos + othma_pos

xtile p_asset = asset [pw=wgt], nq(100)

gen totalassets = safe_assets + vehicles + real_estate + public_equity + private_business + bonds
gen ratio_safe_assets = safe_assets / totalassets if totalassets > 0
gen ratio_vehicles = vehicles / totalassets if totalassets > 0
gen ratio_real_estate = real_estate / totalassets if totalassets > 0
gen ratio_public_equity = public_equity / totalassets if totalassets > 0
gen ratio_private_business = private_business / totalassets if totalassets > 0
gen ratio_bonds = bonds / totalassets if totalassets > 0
// gen total_ratio = ratio_safe_assets + ratio_vehicles + ratio_real_estate + ratio_public_equity + ratio_private_business + ratio_bonds

collapse ///
    (mean) sh_safe_assets=ratio_safe_assets ///
    (mean) sh_vehicles=ratio_vehicles ///
    (mean) sh_real_estate=ratio_real_estate ///
    (mean) sh_public_equity=ratio_public_equity ///
    (mean) sh_private_business=ratio_private_business ///
    (mean) sh_bonds=ratio_bonds ///
    [pw=wgt], by(p_asset)

replace sh_safe_assets = 100 * sh_safe_assets
replace sh_vehicles = 100 * sh_vehicles
replace sh_real_estate = 100 * sh_real_estate
replace sh_public_equity = 100 * sh_public_equity
replace sh_private_business = 100 * sh_private_business
replace sh_bonds = 100 * sh_bonds
// replace sh_total = 100 * sh_total

twoway ///
    (lpoly sh_safe_assets p_asset, bwidth(5) lcolor(navy) lwidth(medthick)) ///
    (lpoly sh_vehicles p_asset, bwidth(5) lcolor(maroon) lwidth(medthick)) ///
    (lpoly sh_real_estate p_asset, bwidth(5) lcolor(green) lwidth(medthick)) ///
    (lpoly sh_public_equity p_asset, bwidth(5) lcolor(orange) lwidth(medthick)) ///
    (lpoly sh_private_business p_asset, bwidth(5) lcolor(brown) lwidth(medthick)) ///
    (lpoly sh_bonds p_asset, bwidth(5) lcolor(gs6) lwidth(medthick)), ///
    xlabel(0(10)100, labsize(small)) ///
    ylabel(0(20)100, angle(horizontal) labsize(small)) ///
    xtitle("Percentile of distribution of total assets") ///
    ytitle("Share of total assets (%)") ///
    legend(order(1 "Safe assets" 2 "Vehicles" 3 "Real estate" 4 "Public equity" 5 "Private business" 6 "Bonds") position(6) ring(1) rows(1)) ///
    title("Asset Class Shares in Household Portfolios (Kernel Smoothed)")
graph export "portfolio_shares_by_asset_percentile.png", replace width(2400)

twoway ///
    (line sh_safe_assets p_asset, lcolor(navy) lwidth(medthick)) ///
    (line sh_vehicles p_asset, lcolor(maroon) lwidth(medthick)) ///
    (line sh_real_estate p_asset, lcolor(green) lwidth(medthick)) ///
    (line sh_public_equity p_asset, lcolor(orange) lwidth(medthick)) ///
    (line sh_private_business p_asset, lcolor(brown) lwidth(medthick)) ///
    (line sh_bonds p_asset, lcolor(gs6) lwidth(medthick)), ///
    xlabel(0(10)100, labsize(small)) ///
    ylabel(0(20)100, angle(horizontal) labsize(small)) ///
    xtitle("Percentile of distribution of total assets") ///
    ytitle("Share of total assets (%)") ///
    legend(order(1 "Safe assets" 2 "Vehicles" 3 "Real estate" 4 "Public equity" 5 "Private business" 6 "Bonds") position(6) ring(1) rows(1)) ///
    title("Asset Class Shares in Household Portfolios (Raw)")
graph export "portfolio_shares_by_asset_percentile_raw.png", replace width(2400)


restore
