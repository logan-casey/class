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

label var mean_fin "Financial assets"
label var mean_asset "Total assets"
label var mean_networth "Net worth"

	
* 5) Plot
twoway ///
    (line mean_asset    pct, lcolor(navy)   lwidth(medthick)) ///
    (line mean_fin      pct, lcolor(maroon) lwidth(medthick)) ///
    (line mean_networth pct, lcolor(forest_green) lwidth(medthick)), ///
	yscale(log) ylabel(1000 "1k" 10000 "10k" 100000 "100k" 1000000 "1m" 10000000 "10m", angle(horizontal) labsize(small)) ///
    xtitle("Percentile (weighted; variable-specific ranking)") ///
    ytitle("Weighted mean") ///
    legend(order(1 "Asset" 2 "Fin" 3 "Net worth")) ///
    title("Means across percentiles")
graph export "wealth_dist.png", replace width(2400)

restore

** FIGURE 2: Participation rates by asset class (x-axis = p_asset)
preserve
* Treat missing/negative components as zero (asset-only construction)
foreach v in checking saving mma call cds savbnd vehic homeq nnresre oresre deq bus govtbnd obnd mortbnd cashli gbmutf obmutf trusts othma {
    gen `v'_pos = cond(missing(`v'), 0, cond(`v' < 0, 0, `v'))
}

gen safe_assets = checking_pos + saving_pos + mma_pos + call_pos + cds_pos + savbnd_pos
gen vehicles = vehic_pos
gen real_estate = homeq_pos + nnresre_pos + oresre_pos
gen public_equity = deq_pos
gen private_business = bus_pos
gen bonds = govtbnd_pos + obnd_pos + mortbnd_pos + cashli_pos + gbmutf_pos + obmutf_pos + trusts_pos + othma_pos

xtile p_asset = asset [pw=wgt], nq(100)

gen hold_safe_assets = safe_assets > 0 if !missing(safe_assets)
gen hold_vehicles = vehicles > 0 if !missing(vehicles)
gen hold_real_estate = real_estate > 0 if !missing(real_estate)
gen hold_public_equity = public_equity > 0 if !missing(public_equity)
gen hold_private_business = private_business > 0 if !missing(private_business)
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
    xtitle("Percentile of distribution of total assets (weighted)") ///
    ytitle("Sample-weighted share with holdings") ///
    legend(order(1 "Safe assets" 2 "Vehicles" 3 "Real estate" 4 "Public equity" 5 "Private business" 6 "Bonds")) ///
    title("Participation Rates by Asset Class (Kernel Smoothed)")
graph export "participation_by_asset_percentile.png", replace width(2400)

restore

** FIGURE 3: Asset class shares in portfolios (x-axis = p_asset)
preserve
* Treat missing/negative components as zero (asset-only construction)
foreach v in checking saving mma call cds savbnd vehic homeq nnresre oresre deq bus govtbnd obnd mortbnd cashli gbmutf obmutf trusts othma {
    gen `v'_pos = cond(missing(`v'), 0, cond(`v' < 0, 0, `v'))
}

gen safe_assets = checking_pos + saving_pos + mma_pos + call_pos + cds_pos + savbnd_pos
gen vehicles = vehic_pos
gen real_estate = homeq_pos + nnresre_pos + oresre_pos
gen public_equity = deq_pos
gen private_business = bus_pos
gen bonds = govtbnd_pos + obnd_pos + mortbnd_pos + cashli_pos + gbmutf_pos + obmutf_pos + trusts_pos + othma_pos

xtile p_asset = asset [pw=wgt], nq(100)

gen ratio_safe_assets = safe_assets / asset if asset > 0
gen ratio_vehicles = vehicles / asset if asset > 0
gen ratio_real_estate = real_estate / asset if asset > 0
gen ratio_public_equity = public_equity / asset if asset > 0
gen ratio_private_business = private_business / asset if asset > 0
gen ratio_bonds = bonds / asset if asset > 0

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

twoway ///
    (lpoly sh_safe_assets p_asset, bwidth(5) lcolor(navy) lwidth(medthick)) ///
    (lpoly sh_vehicles p_asset, bwidth(5) lcolor(maroon) lwidth(medthick)) ///
    (lpoly sh_real_estate p_asset, bwidth(5) lcolor(green) lwidth(medthick)) ///
    (lpoly sh_public_equity p_asset, bwidth(5) lcolor(orange) lwidth(medthick)) ///
    (lpoly sh_private_business p_asset, bwidth(5) lcolor(brown) lwidth(medthick)) ///
    (lpoly sh_bonds p_asset, bwidth(5) lcolor(gs6) lwidth(medthick)), ///
    xlabel(0(10)100, labsize(small)) ///
    ylabel(0(20)100, angle(horizontal) labsize(small)) ///
    xtitle("Percentile of distribution of total assets (weighted)") ///
    ytitle("Share of total assets (%)") ///
    legend(order(1 "Safe assets" 2 "Vehicles" 3 "Real estate" 4 "Public equity" 5 "Private business" 6 "Bonds")) ///
    title("Asset Class Shares in Household Portfolios (Kernel Smoothed)")
graph export "portfolio_shares_by_asset_percentile.png", replace width(2400)

restore
