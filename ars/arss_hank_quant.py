"""Quantitative open-economy HANK model from
Exchange Rates and Monetary Policy with Heterogeneous Agents:
Sizing up the Real Income Channel,
by Auclert, Rognlie, Souchier, and Straub

strategy: start with SSJ notebooks/IKC rep kit, translate equations into simple blocks

paper calibration, IRFs (to foreign interest-rate shock in sec 5.3 (AR(1), persistence 0.85,
normalized to a 1% impact depreciation in Q)).

equations used:
- (13), (18), (20), (21), (47)-(53)
- UIP and Fisher relations from core model
- (48) for het labor-income incidence in HH block
- (23) to construct an NFA path from net exports
"""

from __future__ import annotations # recommended for performance

import numpy as np
import sequence_jacobian as sj
from sequence_jacobian import grids
from pathlib import Path

# -----------------------------------------------------------------------------
# DEFINE MODEL BLOCKS
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Household block
# -----------------------------------------------------------------------------

def make_grids_oe(rho_e, sd_e, n_e, min_a, max_a, n_a):
    """Rouwenhorst idiosyncratic income process + asset grid."""
    e_grid, pi_e, Pi = grids.markov_rouwenhorst(rho=rho_e, sigma=sd_e, N=n_e)
    a_grid = grids.asset_grid(min_a, max_a, n_a)
    return e_grid, pi_e, Pi, a_grid


def income_incidence(W, P, N, e_grid, pi_e, zeta, W_ss, P_ss, N_ss):
    """Eq. (48): heterogeneous labor-income incidence to aggregate labor income."""
    agg_ratio = (W / P * N) / (W_ss / P_ss * N_ss)
    expo = 1.0 + zeta * np.log(agg_ratio)
    num = e_grid ** expo
    denom = np.sum(pi_e * num)
    y = (W / P) * N * num / denom
    return y


def build_household_block():
    hh = sj.hetblocks.hh_sim.hh
    hh_oe = hh.add_hetinputs([make_grids_oe, income_incidence])
    return hh_oe.rename("hh_oe")


# -----------------------------------------------------------------------------
# Simple blocks
# -----------------------------------------------------------------------------


@sj.simple
def production(Y):
    """Eq. (13): Y = N."""
    N = Y
    return N


@sj.simple
def exchange_rate_identity(E, Q, P):
    """Eq. (6): E = Q * P (with P* = 1)"""
    exrate_res = E - Q * P
    return exrate_res


@sj.simple
def import_price_pass_through(E):
    """Full pass-through for imports (theta_F = 0): PF = E."""
    PF = E
    return PF


@sj.simple
def price_index(PF, PH, alpha, eta):
    """Eq. (4): price index"""
    P = (alpha * PF ** (1.0 - eta) + (1.0 - alpha) * PH ** (1.0 - eta)) ** (1.0 / (1.0 - eta))
    return P


@sj.simple
def inflation_defs(P, PH, PH_star, W):
    """Inflation definitions for CPI, home prices, export prices, and wages."""
    pi = P / P(-1) - 1.0
    piH = PH / PH(-1) - 1.0
    piH_star = PH_star / PH_star(-1) - 1.0
    piw = W / W(-1) - 1.0
    return pi, piH, piH_star, piw


@sj.simple
def foreign_rate_block(i_star):
    """With P*=1 and no foreign inflation: r* = i*."""
    r_star = i_star
    return r_star


@sj.simple
def monetary_rule_taylor(i, pi, r_ss, eps_m, rho_m, phi_pi):
    """Eq. (20): MP rule."""
    mp_taylor_res = i - (rho_m * i(-1) + (1.0 - rho_m) * (r_ss + phi_pi * pi(+1)) + eps_m)
    return mp_taylor_res


@sj.simple
def monetary_rule_const_r(i, r_ss, pi, eps_m):
    """Base-model policy rule with a constant real-rate anchor."""
    mp_res = i - (r_ss + pi(+1) + eps_m)
    return mp_res


@sj.simple
def fisher_equation(i, r, pi):
    """(unnumbered, linearized) i = r - pi_{t+1}."""
    fisher_res = i - r - pi(+1)
    return fisher_res


@sj.simple
def uip_real(r, r_star, Q):
    """Eq. (7): real UIP."""
    uip_res = (1.0 + r) - (1.0 + r_star) * Q(+1) / Q
    return uip_res


@sj.simple
def piH_nkpc(piH, W, PH, Z, mu, kappa_H, r):
    """Eq. (49): domestic-price Phillips curve."""
    piH_res = piH - (kappa_H * (mu * W / (Z * PH) - 1.0) + piH(+1) / (1.0 + r))
    return piH_res


@sj.simple
def piHstar_nkpc(piH_star, PH, E, PH_star, mu_H_star, kappa_H_star, r):
    """Eq. (51): export-price Phillips curve."""
    piH_star_res = piH_star - (kappa_H_star * (mu_H_star * PH / (E * PH_star) - 1.0) + piH_star(+1) / (1.0 + r))
    return piH_star_res


@sj.simple
def xH_target_law(xH_hat, PH, P, eta, beta, theta_sub, xH_ss, PH_ss, P_ss):
    """Eq. (52): target domestic share in log-deviations from steady state."""
    xH_target_res = (xH_hat / xH_ss).apply(np.log) + (1.0 - beta * theta_sub) * eta * ((PH / P) / (PH_ss / P_ss)).apply(np.log) - beta * theta_sub * (xH_hat(+1) / xH_ss).apply(np.log)
    return xH_target_res


@sj.simple
def xHstar_target_law(xH_star_hat, PH_star, gamma, beta_star, theta_sub, xH_star_ss, PH_star_ss):
    """Eq. (53): target foreign home-good share in log-deviations from steady state."""
    xH_star_target_res = (xH_star_hat / xH_star_ss).apply(np.log) + (1.0 - beta_star * theta_sub) * gamma * (PH_star / PH_star_ss).apply(np.log) - beta_star * theta_sub * (xH_star_hat(+1) / xH_star_ss).apply(np.log)
    return xH_star_target_res


@sj.simple
def delayed_substitution(shareH, shareH_star, xH_hat, xH_star_hat, theta_sub, C, C_star, shareH_ss, shareH_star_ss, xH_ss, xH_star_ss):
    """Eq. (54)-(55): sluggish adjustment of actual shares in log-deviations."""
    shareH_res = (shareH / shareH_ss).apply(np.log) - (theta_sub * (shareH(-1) / shareH_ss).apply(np.log) + (1.0 - theta_sub) * (xH_hat / xH_ss).apply(np.log))
    shareH_star_res = (shareH_star / shareH_star_ss).apply(np.log) - (theta_sub * (shareH_star(-1) / shareH_star_ss).apply(np.log) + (1.0 - theta_sub) * (xH_star_hat / xH_star_ss).apply(np.log))
    CH = shareH * C
    CH_star = shareH_star * C_star
    return shareH_res, shareH_star_res, CH, CH_star


@sj.simple
def goods_market(CH, CH_star, Y):
    """Eq. (21): goods market clearing."""
    goods_mkt = CH + CH_star - Y
    return goods_mkt


@sj.simple
def union_wage_nkpc(piw, N, W, P, C, psi_labor, phi_labor, sigma, mu_w, kappa_w, beta):
    """Eq. (18): wage Phillips curve."""
    uprime = C ** (-sigma)
    vprime = psi_labor * N ** phi_labor
    rhs = kappa_w * (vprime / ((W / (mu_w * P)) * uprime) - 1.0) + beta * piw(+1)
    wnkpc = piw - rhs
    return wnkpc


@sj.simple
def external_accounts(PH, P, Y, C):
    """Eq. (23) trade-balance component: NX = PH/P * Y - C."""
    NX = PH / P * Y - C
    return NX


# -----------------------------------------------------------------------------
# Calibration and model builder
# -----------------------------------------------------------------------------

def base_calibration():
    """Quarterly calibration from Table 2 + section 5.3 shock settings."""
    beta_star = 0.990
    r_ss = 1.0 / beta_star - 1.0

    theta_w = 0.938
    theta_H = 0
    theta_H_star = 0
    alpha = 0.4

    calib = {
        # Preferences / demand
        "sigma": 1.0,
        "eis": 1.0,
        "phi_labor": 2.0,
        "beta": 0.965,
        "beta_star": beta_star,
        "alpha": alpha,
        "eta": 1/(2-alpha), # varies ...
        "gamma": 1/(2-alpha),
        "cbar": 0,  # (not directly needed)
        # Pricing / markups
        "mu": 1.043,
        # Phillips curve coefficients
        "theta_w": theta_w,
        "theta_H": theta_H,
        "theta_H_star": theta_H_star,
        "kappa_w": (1.0 - 0.962 * theta_w) * (1.0 - theta_w) / theta_w,
        "kappa_H": (1.0 - theta_H) * (1.0 - theta_H / (1.0 + r_ss)) / theta_H,
        "kappa_H_star": (1.0 - theta_H_star) * (1.0 - theta_H_star / (1.0 + r_ss)) / theta_H_star,
        # Monetary policy
        "eps_m": 0.0,
        "r_ss": r_ss,
        # Productivity
        "Z": 1.0,
        # Foreign block normalization
        "C_star": 1.0,
        # Idiosyncratic earnings process
        "rho_e": 0.912,
        "sd_e": 0.883,
        "n_e": 7,
        # Asset grid
        "min_a": 0.0,
        "max_a": 1000.0,
        "n_a": 250,
        # Steady-state anchors
        "W_ss": 1.0 / 1.041,
        "P_ss": 1.0,
        "N_ss": 1.0,
        # Shock process (section 5.3)
        "rho_i_star": 0.85,
        "q0_target": 0.01,
    }
    return calib


def quantitative_calibration():
    """Quarterly calibration from Table 2 + section 5.3 shock settings."""
    beta_star = 0.990
    r_ss = 1.0 / beta_star - 1.0

    theta_w = 0.938
    theta_H = 0.66
    theta_H_star = 0.66

    calib = {
        # Preferences / demand
        "sigma": 1.0,
        "eis": 1.0,
        "phi_labor": 2.0,
        "beta": 0.962,
        "beta_star": beta_star,
        "alpha": 0.344,
        "eta": 4.0,
        "gamma": 4.0,
        "theta_sub": 0.976,
        "cbar": 0.085,  # (not directly needed)
        "zeta": -0.196,
        # Pricing / markups
        "mu": 1.041,
        "mu_w": 1.0,
        "mu_H_star": 1.0,
        # Phillips curve coefficients
        "theta_w": theta_w,
        "theta_H": theta_H,
        "theta_H_star": theta_H_star,
        "kappa_w": (1.0 - 0.962 * theta_w) * (1.0 - theta_w) / theta_w,
        "kappa_H": (1.0 - theta_H) * (1.0 - theta_H / (1.0 + r_ss)) / theta_H,
        "kappa_H_star": (1.0 - theta_H_star) * (1.0 - theta_H_star / (1.0 + r_ss)) / theta_H_star,
        # Monetary policy
        "rho_m": 0.8,
        "phi_pi": 1.5,
        "eps_m": 0.0,
        "r_ss": r_ss,
        # Productivity
        "Z": 1.0,
        # Foreign block normalization
        "C_star": 1.0,
        # Idiosyncratic earnings process
        "rho_e": 0.912,
        "sd_e": 0.883,
        "n_e": 7,
        # Asset grid
        "min_a": 0.0,
        "max_a": 1000.0,
        "n_a": 250,
        # Steady-state anchors
        "W_ss": 1.0 / 1.041,
        "P_ss": 1.0,
        "N_ss": 1.0,
        # Shock process (section 5.3)
        "rho_i_star": 0.85,
        "q0_target": 0.01,
    }
    return calib


def build_quant_model(policy_rule="const_r"):
    """Build the quantitative open-economy model.

    policy_rule: "taylor" or "const_r"
    """
    if policy_rule not in {"taylor", "const_r"}:
        raise ValueError(f"Unknown policy_rule='{policy_rule}'. Use 'taylor' or 'const_r'.")

    hh_oe = build_household_block()
    mp_block = monetary_rule_taylor if policy_rule == "taylor" else monetary_rule_const_r
    blocks = [
        hh_oe,
        production,
        exchange_rate_identity,
        import_price_pass_through,
        price_index,
        inflation_defs,
        foreign_rate_block,
        mp_block,
        fisher_equation,
        uip_real,
        piH_nkpc,
        piHstar_nkpc,
        xH_target_law,
        xHstar_target_law,
        delayed_substitution,
        goods_market,
        union_wage_nkpc,
        external_accounts,
    ]
    model = sj.create_model(blocks, name=f"HA_OE_Quant_{policy_rule}")
    return model, hh_oe


def build_steady_state(model, hh_oe, calib):
    """Assemble a consistent steady state under the paper calibration."""
    hh_calib = {
        "rho_e": calib["rho_e"],
        "sd_e": calib["sd_e"],
        "n_e": calib["n_e"],
        "min_a": calib["min_a"],
        "max_a": calib["max_a"],
        "n_a": calib["n_a"],
        "r": calib["r_ss"],
        "beta": calib["beta"],
        "eis": calib["eis"],
        "W": calib["W_ss"],
        "P": calib["P_ss"],
        "N": calib["N_ss"],
        "zeta": calib["zeta"],
        "W_ss": calib["W_ss"],
        "P_ss": calib["P_ss"],
        "N_ss": calib["N_ss"],
    }
    hh_ss = hh_oe.steady_state(hh_calib)
    C_ss = float(hh_ss["C"])

    # Calibrate psi_labor to satisfy wage NKPC at steady state.
    psi_labor = (calib["W_ss"] / (calib["mu_w"] * calib["P_ss"])) * (C_ss ** (-calib["sigma"])) / (calib["N_ss"] ** calib["phi_labor"])

    # Match average import share of 40% and goods-market clearing at Y_ss = 1.
    shareH_ss = 0.60
    CH_ss = shareH_ss * C_ss
    CH_star_ss = 1.0 - CH_ss
    shareH_star_ss = CH_star_ss / calib["C_star"]
    xH_ss = shareH_ss
    xH_star_ss = shareH_star_ss

    ss_guess = {
        **calib,
        "psi_labor": psi_labor,
        "Y": 1.0,
        "W": calib["W_ss"],
        "PH": 1.0,
        "PH_star": 1.0,
        "P": 1.0,
        "Q": 1.0,
        "E": 1.0,
        "i": calib["r_ss"],
        "r": calib["r_ss"],
        "i_star": calib["r_ss"],
        "xH_hat": shareH_ss,
        "xH_star_hat": shareH_star_ss,
        "shareH": shareH_ss,
        "shareH_star": shareH_star_ss,
        "xH_ss": xH_ss,
        "xH_star_ss": xH_star_ss,
        "shareH_ss": shareH_ss,
        "shareH_star_ss": shareH_star_ss,
        "PH_ss": 1.0,
        "PH_star_ss": 1.0,
    }
    ss = model.steady_state(ss_guess)
    return ss


# -----------------------------------------------------------------------------
# IRFs
# -----------------------------------------------------------------------------


def _linear_nfa_from_nx(dnx, r_ss):
    """Linearized Eq. (23) around nfa_ss = 0."""
    T = len(dnx)
    dnfa = np.zeros(T)
    for t in range(T):
        prev = dnfa[t - 1] if t > 0 else 0.0
        dnfa[t] = (1.0 + r_ss) * prev + dnx[t]
    return dnfa


def solve_exchange_rate_irf(T=40, policy_rule="taylor"):
    """Run section-5.3-style foreign-rate shock IRF, normalized to dQ0 = 1%.

    policy_rule: "taylor" or "const_r"
    """
    model, hh_oe = build_quant_model(policy_rule=policy_rule)
    calib = quantitative_calibration()
    ss = build_steady_state(model, hh_oe, calib)

    unknowns = ["Y", "W", "PH", "PH_star", "Q", "E", "i", "r", "xH_hat", "xH_star_hat", "shareH", "shareH_star"]
    mp_target = "mp_taylor_res" if policy_rule == "taylor" else "mp_res"
    targets = [
        "goods_mkt",
        "wnkpc",
        "piH_res",
        "piH_star_res",
        mp_target,
        "fisher_res",
        "uip_res",
        "xH_target_res",
        "xH_star_target_res",
        "shareH_res",
        "shareH_star_res",
        "exrate_res",
    ]

    outputs = [
        "Y",
        "C",
        "W",
        "Q",
        "E",
        "r",
        "i",
        "PH",
        "P",
        "PH_star",
        "CH",
        "CH_star",
        "NX",
        "shareH",
        "shareH_star",
        "xH_hat",
        "xH_star_hat",
        "pi",
        "piH",
        "piH_star",
        "piw",
        "goods_mkt",
        "wnkpc",
    ]

    rho = calib["rho_i_star"]
    base = rho ** np.arange(T)

    # First pass for normalization.
    irf_raw = model.solve_impulse_linear(ss, unknowns, targets, {"i_star": 1e-4 * base}, outputs=outputs)
    q0 = float(irf_raw["Q"][0])
    if np.isclose(q0, 0.0):
        raise RuntimeError("Q impact is numerically zero; cannot normalize shock to target dQ0.")
    scale = calib["q0_target"] / q0

    # Final, normalized IRF.
    irf = model.solve_impulse_linear(ss, unknowns, targets, {"i_star": 1e-4 * scale * base}, outputs=outputs)

    # Ex-post NFA path from linearized current account identity.
    dnfa = _linear_nfa_from_nx(irf["NX"], calib["r_ss"])

    pct = {
        "Y_pct": 100.0 * irf["Y"] / ss["Y"],
        "C_pct": 100.0 * irf["C"] / ss["C"],
        "Q_pct": 100.0 * irf["Q"] / ss["Q"],
        "NX_pctY": 100.0 * irf["NX"] / ss["Y"],
        "NFA_pctY": 100.0 * dnfa / ss["Y"],
        "r_pp": 100.0 * irf["r"],
        "i_pp": 100.0 * irf["i"],
    }

    return {
        "model": model,
        "policy_rule": policy_rule,
        "calibration": calib,
        "ss": ss,
        "irf": irf,
        "derived": pct,
        "shock_scale": scale,
    }


def plot_exchange_rate_irf_figure(result, T_plot=32, savepath="figures/ha_oe_quant_irf.png"):
    """Save a paper-style 8-panel IRF figure."""
    import matplotlib.pyplot as plt

    ss = result["ss"]
    irf = result["irf"]
    d = result["derived"]

    T_plot = min(T_plot, len(d["Y_pct"]))
    t = np.arange(T_plot)

    # Recover levels from linear responses where needed.
    W = ss["W"] + irf["W"]
    P = ss["P"] + irf["P"]
    Y = ss["Y"] + irf["Y"]  # N=Y from Eq. (13)
    E = ss["E"] + irf["E"]
    PH = ss["PH"] + irf["PH"]
    PH_star = ss["PH_star"] + irf["PH_star"]
    CH_star = (ss["CH_star"] + irf["CH_star"])

    # Real wage income and dividends (Eq. 15 object).
    wage_income = W / P * Y
    wage_income_ss = ss["W"] / ss["P"] * ss["Y"]
    wage_income_pctY = 100.0 * (wage_income - wage_income_ss) / ss["Y"]

    dividends = (PH * Y - W * Y) / P + (E * PH_star - PH) / P * CH_star
    dividends_ss = (ss["PH"] * ss["Y"] - ss["W"] * ss["Y"]) / ss["P"] + (ss["E"] * ss["PH_star"] - ss["PH"]) / ss["P"] * ss["CH_star"]
    dividends_pctY = 100.0 * (dividends - dividends_ss) / ss["Y"]

    series = [
        (d["Y_pct"][:T_plot], "Output", "Percent of Yss"),
        (d["C_pct"][:T_plot], "Consumption", "Percent of Css"),
        (d["NX_pctY"][:T_plot], "Net exports", "Percent of Yss"),
        (d["NFA_pctY"][:T_plot], "NFA", "Percent of Yss"),
        (wage_income_pctY[:T_plot], "Real wage income", "Percent of Yss"),
        (dividends_pctY[:T_plot], "Dividends", "Percent of Yss"),
        (d["Q_pct"][:T_plot], "Real exchange rate", "Percent"),
        (d["r_pp"][:T_plot], "Real interest rate", "Percentage points"),
    ]

    fig, axes = plt.subplots(2, 4, figsize=(14, 7), constrained_layout=True)
    axes = axes.ravel()

    for ax, (y, title, ylabel) in zip(axes, series):
        ax.plot(t, y, color="#2F6B4F", lw=2.2)
        ax.axhline(0.0, color="#666666", lw=0.8, alpha=0.6)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Quarters")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25, lw=0.6)

    savepath = Path(savepath)
    savepath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(savepath, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return savepath


if __name__ == "__main__":
    out = solve_exchange_rate_irf(T=40)
    d = out["derived"]
    fig_path = plot_exchange_rate_irf_figure(out, T_plot=32, savepath="figures/ha_oe_quant_irf.png")
    print("Shock scale:", out["shock_scale"])
    print("Q impact (%):", d["Q_pct"][0])
    print("Y first 8 qtrs (%):", np.round(d["Y_pct"][:8], 4))
    print("C first 8 qtrs (%):", np.round(d["C_pct"][:8], 4))
    print("NFA first 8 qtrs (% of Y):", np.round(d["NFA_pctY"][:8], 4))
    print("Saved figure:", fig_path)
