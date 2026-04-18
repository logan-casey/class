"""Base open-economy HA model with paper benchmark calibration.

Implements the baseline system (without quantitative extensions in Eqs. (47)-(53)),
using:
- Eq. (13): production Y = N
- Eq. (18): wage Phillips curve
- Eq. (19): simple policy rule (constant real rate rule)
- Eq. (21): goods market clearing
- Eq. (23): current-account accounting (used ex post for NFA path)
- Real UIP and Fisher relations

The script computes IRFs to an AR(1) foreign interest-rate shock (rho=0.85),
normalized so that the impact depreciation in Q is 1%.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import sequence_jacobian as sj
from sequence_jacobian import grids


# -----------------------------------------------------------------------------
# Household block (Eq. 9 in canonical one-asset form)
# -----------------------------------------------------------------------------


def make_grids_oe(rho_e, sd_e, n_e, min_a, max_a, n_a):
    """Rouwenhorst idiosyncratic process and asset grid."""
    e_grid, pi_e, Pi = grids.markov_rouwenhorst(rho=rho_e, sigma=sd_e, N=n_e)
    a_grid = grids.asset_grid(min_a, max_a, n_a)
    return e_grid, pi_e, Pi, a_grid


def income_baseline(W, P, N, e_grid):
    """Baseline labor income without Eq. (48) incidence extension."""
    y = (W / P) * N * e_grid
    return y


def build_household_block():
    hh = sj.hetblocks.hh_sim.hh
    hh_oe = hh.add_hetinputs([make_grids_oe, income_baseline])
    return hh_oe.rename("hh_oe")


# -----------------------------------------------------------------------------
# Base model blocks
# -----------------------------------------------------------------------------


@sj.simple
def production(Y):
    """Eq. (13): Y = N."""
    N = Y
    return N


@sj.simple
def home_price_flex(W, mu):
    """Eq. (14): flexible-price domestic goods price."""
    PH = mu * W
    return PH


@sj.simple
def import_price_pass_through(E):
    """Baseline PCP import pricing: PF = E."""
    PF = E
    return PF


@sj.simple
def law_of_one_price(PH, E):
    """Eq. (17): P*_H = PH / E."""
    PH_star = PH / E
    return PH_star


@sj.simple
def foreign_demand_home_goods(PH_star, C_star, alpha, gamma):
    """Eq. (12) with P* normalized to 1."""
    CH_star = alpha * PH_star ** (-gamma) * C_star
    return CH_star


@sj.simple
def price_index(PF, PH, alpha, eta):
    """Eq. (4): CPI."""
    P = (alpha * PF ** (1.0 - eta) + (1.0 - alpha) * PH ** (1.0 - eta)) ** (1.0 / (1.0 - eta))
    return P


@sj.simple
def inflation_defs(P, W):
    """CPI and wage inflation definitions."""
    pi = P / P(-1) - 1.0
    piw = W / W(-1) - 1.0
    return pi, piw


@sj.simple
def exchange_rate_identity(E, Q, P):
    """Eq. (6) with P* = 1: E = Q * P."""
    exrate_res = E - Q * P
    return exrate_res


@sj.simple
def foreign_rate_block(i_star):
    """With P* constant: r* = i*."""
    r_star = i_star
    return r_star


@sj.simple
def monetary_rule_const_r(i, r_ss, pi, eps_m):
    """Eq. (19): i_t = r_ss + pi_{t+1} + eps_t."""
    mp_res = i - (r_ss + pi(+1) + eps_m)
    return mp_res


@sj.simple
def fisher_equation(i, r, pi):
    """1 + i = (1 + r)(1 + pi_{t+1})."""
    fisher_res = 1.0 + i - (1.0 + r) * (1.0 + pi(+1))
    return fisher_res


@sj.simple
def uip_real(r, r_star, Q):
    """Eq. (7): real UIP."""
    uip_res = (1.0 + r) - (1.0 + r_star) * Q(+1) / Q
    return uip_res


@sj.simple
def hh_goods_split(C, PH, P, alpha, eta):
    """Eq. (11): domestic demand component from CES intratemporal allocation."""
    CH = (1.0 - alpha) * (PH / P) ** (-eta) * C
    return CH


@sj.simple
def goods_market(CH, CH_star, Y):
    """Eq. (21): CH + CH* = Y."""
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
    """Net exports in consumption units used in Eq. (23)."""
    NX = PH / P * Y - C
    return NX


@sj.simple
def dividends(PH, Y, W, N, P, E, PH_star, CH_star):
    """Eq. (15): real dividends."""
    Div = (PH * Y - W * N) / P + (E * PH_star - PH) / P * CH_star
    return Div


# -----------------------------------------------------------------------------
# Calibration and assembly
# -----------------------------------------------------------------------------


def base_calibration():
    """Benchmark quarterly calibration (Table 2 baseline values)."""
    beta_star = 0.990
    r_ss = 1.0 / beta_star - 1.0
    alpha = 0.4
    eta = 2.0 - alpha

    theta_w = 0.938
    beta = 0.965

    calib = {
        # Preferences
        "sigma": 1.0,
        "eis": 1.0,
        "phi_labor": 2.0,
        "beta": beta,
        "beta_star": beta_star,
        "alpha": alpha,
        "eta": eta,
        "gamma": eta,
        # Prices / markups
        "mu": 1.043,
        "mu_w": 1.0,
        # Wage Phillips slope
        "theta_w": theta_w,
        "kappa_w": (1.0 - beta * theta_w) * (1.0 - theta_w) / theta_w,
        # Policy
        "eps_m": 0.0,
        "r_ss": r_ss,
        # Foreign block normalization
        "C_star": 1.0,
        # Idiosyncratic process
        "rho_e": 0.912,
        "sd_e": 0.883,
        "n_e": 7,
        # Asset grid
        "min_a": 0.0,
        "max_a": 1000.0,
        "n_a": 250,
        # Steady-state anchors
        "W_ss": 1.0 / 1.043,
        "P_ss": 1.0,
        "N_ss": 1.0,
        # Shock process
        "rho_i_star": 0.85,
        "q0_target": 0.01,
    }
    return calib


def build_base_model():
    """Build baseline DAG (without Eqs. 47-53)."""
    hh_oe = build_household_block()
    blocks = [
        hh_oe,
        production,
        home_price_flex,
        import_price_pass_through,
        law_of_one_price,
        foreign_demand_home_goods,
        price_index,
        inflation_defs,
        exchange_rate_identity,
        foreign_rate_block,
        monetary_rule_const_r,
        fisher_equation,
        uip_real,
        hh_goods_split,
        goods_market,
        union_wage_nkpc,
        external_accounts,
        dividends,
    ]
    model = sj.create_model(blocks, name="HA_OE_Base")
    return model, hh_oe


def build_steady_state(model, hh_oe, calib):
    """Construct steady state for baseline model."""
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
    }
    hh_ss = hh_oe.steady_state(hh_calib)
    C_ss = float(hh_ss["C"])

    psi_labor = (calib["W_ss"] / (calib["mu_w"] * calib["P_ss"])) * (C_ss ** (-calib["sigma"])) / (calib["N_ss"] ** calib["phi_labor"])

    ss_guess = {
        **calib,
        "psi_labor": psi_labor,
        "Y": 1.0,
        "W": calib["W_ss"],
        "PH": 1.0,
        "PF": 1.0,
        "PH_star": 1.0,
        "P": 1.0,
        "Q": 1.0,
        "E": 1.0,
        "i": calib["r_ss"],
        "r": calib["r_ss"],
        "i_star": calib["r_ss"],
    }
    ss = model.steady_state(ss_guess)
    return ss


# -----------------------------------------------------------------------------
# IRFs
# -----------------------------------------------------------------------------


def _linear_nfa_from_nx(dnx, r_ss):
    """Linearized Eq. (23) around nfa_ss=0."""
    T = len(dnx)
    dnfa = np.zeros(T)
    for t in range(T):
        prev = dnfa[t - 1] if t > 0 else 0.0
        dnfa[t] = (1.0 + r_ss) * prev + dnx[t]
    return dnfa


def solve_exchange_rate_irf(T=40):
    """Baseline IRFs to foreign rate shock, normalized to dQ0 = 1%."""
    model, hh_oe = build_base_model()
    calib = base_calibration()
    ss = build_steady_state(model, hh_oe, calib)

    unknowns = ["Y", "W", "Q", "E", "i", "r"]
    targets = ["goods_mkt", "wnkpc", "exrate_res", "mp_res", "fisher_res", "uip_res"]

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
        "pi",
        "piw",
        "Div",
        "goods_mkt",
        "wnkpc",
    ]

    rho = calib["rho_i_star"]
    base = rho ** np.arange(T)

    irf_raw = model.solve_impulse_linear(ss, unknowns, targets, {"i_star": 1e-4 * base}, outputs=outputs)
    q0 = float(irf_raw["Q"][0])
    if np.isclose(q0, 0.0):
        raise RuntimeError("Q impact is numerically zero; cannot normalize shock to target dQ0.")
    scale = calib["q0_target"] / q0

    irf = model.solve_impulse_linear(ss, unknowns, targets, {"i_star": 1e-4 * scale * base}, outputs=outputs)
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

    return {"model": model, "calibration": calib, "ss": ss, "irf": irf, "derived": pct, "shock_scale": scale}


def plot_exchange_rate_irf_figure(result, T_plot=32, savepath="figures/ha_oe_base_irf.png"):
    """Save 8-panel baseline IRF figure."""
    import matplotlib.pyplot as plt

    ss = result["ss"]
    irf = result["irf"]
    d = result["derived"]

    T_plot = min(T_plot, len(d["Y_pct"]))
    t = np.arange(T_plot)

    W = ss["W"] + irf["W"]
    P = ss["P"] + irf["P"]
    Y = ss["Y"] + irf["Y"]

    wage_income = W / P * Y
    wage_income_ss = ss["W"] / ss["P"] * ss["Y"]
    wage_income_pctY = 100.0 * (wage_income - wage_income_ss) / ss["Y"]

    dividends_pctY = 100.0 * irf["Div"] / ss["Y"]

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
        ax.plot(t, y, color="#1F5E43", lw=2.2)
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
    fig_path = plot_exchange_rate_irf_figure(out, T_plot=32, savepath="figures/ha_oe_base_irf.png")

    print("Shock scale:", out["shock_scale"])
    print("Q impact (%):", d["Q_pct"][0])
    print("Y first 8 qtrs (%):", np.round(d["Y_pct"][:8], 4))
    print("C first 8 qtrs (%):", np.round(d["C_pct"][:8], 4))
    print("NFA first 8 qtrs (% of Y):", np.round(d["NFA_pctY"][:8], 4))
    print("Saved figure:", fig_path)
