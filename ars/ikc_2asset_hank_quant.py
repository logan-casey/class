"""HA open-economy SSJ scaffold mapped to ha_oe.pdf equations.

This module is a block-level blueprint for equations (1)-(23) and (47)-(53).
It is intentionally lightweight: it defines block names, equation residuals,
and a model assembly routine, while leaving calibration and HH grids to the
caller.
"""

from __future__ import annotations

import numpy as np
import sequence_jacobian as sj


# ---------------------------------------------------------------------------
# Core equations (1)-(23)
# ---------------------------------------------------------------------------


@sj.simple
def consumption_aggregator(cF, cH, alpha, eta):
    """Eq. (3): CES consumption aggregator."""
    c = (
        alpha ** (1 / eta) * cF ** ((eta - 1) / eta)
        + (1 - alpha) ** (1 / eta) * cH ** ((eta - 1) / eta)
    ) ** (eta / (eta - 1))
    return c


@sj.simple
def price_index(PF, PH, alpha, eta):
    """Eq. (4): CPI index."""
    P = (alpha * PF ** (1 - eta) + (1 - alpha) * PH ** (1 - eta)) ** (1 / (1 - eta))
    return P


@sj.simple
def hh_goods_split(c, PF, PH, P, alpha, eta):
    """Eqs. (10)-(11): intratemporal demand split."""
    cF = alpha * (PF / P) ** (-eta) * c
    cH = (1 - alpha) * (PH / P) ** (-eta) * c
    return cF, cH


@sj.simple
def real_exchange_rate(E, P_star, P):
    """Eq. (6): real exchange rate."""
    Q = E * P_star / P
    return Q


@sj.simple
def nominal_asset_pricing(i, i_star, E, P_stock_nom, D_nom):
    """Eq. (5): nominal UIP + nominal stock pricing residuals."""
    uip_nom = (1 + i) - (1 + i_star) * E(+1) / E
    equity_nom = (1 + i) - (P_stock_nom(+1) + D_nom(+1)) / P_stock_nom
    return uip_nom, equity_nom


@sj.simple
def real_asset_pricing(r, r_star, Q, p, d):
    """Eq. (7): real UIP + real stock pricing residuals."""
    uip_real = (1 + r) - (1 + r_star) * Q(+1) / Q
    equity_real = (1 + r) - (p(+1) + d(+1)) / p
    return uip_real, equity_real


@sj.simple
def asset_transition_defs(P_stock, s_prime, BH_prime, BF_prime, P, E, D, s, BH, BF, i_lag, i_star_lag):
    """Eq. (8): definitions of a_prime and a_pre."""
    a_prime = (P_stock * s_prime + BH_prime + E * BF_prime) / P
    a_pre = ((P_stock + D) * s + (1 + i_lag) * BH + (1 + i_star_lag) * E * BF) / P
    return a_prime, a_pre


@sj.simple
def foreign_demand_home_goods(PH_star, P_star, C_star, alpha, gamma):
    """Eq. (12): ROW demand for home goods."""
    CH_star = alpha * (PH_star / P_star) ** (-gamma) * C_star
    return CH_star


@sj.simple
def production(N):
    """Eq. (13): production technology."""
    Y = N
    return Y


@sj.simple
def home_price_flex(W, mu):
    """Eq. (14): flexible-price markup."""
    PH = mu * W
    return PH


@sj.simple
def dividends(PH, Y, W, N, P, E, PH_star, CH_star):
    """Eq. (15): real dividends with unhedged export term."""
    d = (PH * Y - W * N) / P + (E * PH_star - PH) / P * CH_star
    return d


@sj.simple
def foreign_policy(beta_star, B):
    """Eq. (16): foreign interest-rate process."""
    i_star = B / (beta_star * B(+1)) - 1
    r_star = i_star
    return i_star, r_star


@sj.simple
def law_of_one_price(PH, E):
    """Eq. (17): law of one price in baseline PCP case."""
    PH_star = PH / E
    return PH_star


@sj.simple
def wage_inflation(W):
    """Definition around Eq. (18): wage inflation."""
    piw = W / W(-1) - 1
    return piw


@sj.simple
def union_wage_nkpc(piw, N, W, P, UCE, vphi, inv_frisch, mu_w, kappa_w, beta):
    """Eq. (18): wage Phillips curve residual."""
    vprime = vphi * N ** (1 + inv_frisch)
    rhs = kappa_w * (vprime / ((W / (mu_w * P)) * UCE) - 1) + beta * piw(+1)
    wnkpc = piw - rhs
    return wnkpc


@sj.simple
def monetary_rule_const_r(i, r_ss, pi_cpi, eps_m):
    """Eq. (19): constant-real-rate rule residual."""
    mp_res = i - (r_ss + pi_cpi(+1) + eps_m)
    return mp_res


@sj.simple
def monetary_rule_taylor(i, i_lag, r_ss, pi_cpi, eps_m, rho_m, phi):
    """Eq. (20): inertial Taylor rule residual."""
    rhs = rho_m * i_lag + (1 - rho_m) * (r_ss + phi * pi_cpi(+1)) + eps_m
    mp_taylor_res = i - rhs
    return mp_taylor_res


@sj.simple
def goods_market_clearing(CH, CH_star, Y):
    """Eq. (21): goods market clearing residual."""
    goods_mkt = CH + CH_star - Y
    return goods_mkt


@sj.simple
def nfa_accounting(A, p):
    """Eq. (22): net foreign assets accounting."""
    nfa = A - p
    return nfa


@sj.simple
def current_account(nfa, PH, P, Y, C, r):
    """Eq. (23): current account identity residual."""
    nx = (PH / P) * Y - C
    ca_res = nfa - nfa(-1) - (nx + r(-1) * nfa(-1))
    return nx, ca_res


# ---------------------------------------------------------------------------
# Quantitative extensions (47)-(53)
# ---------------------------------------------------------------------------


@sj.simple
def consumption_aggregator_stone_geary(cF, cH, cbar, alpha, eta):
    """Eq. (47): Stone-Geary CES consumption basket."""
    c = (
        alpha ** (1 / eta) * (cF - cbar) ** ((eta - 1) / eta)
        + (1 - alpha) ** (1 / eta) * cH ** ((eta - 1) / eta)
    ) ** (eta / (eta - 1))
    return c


@sj.simple
def income_incidence(e, W, P, N, W_ss, P_ss, N_ss, zeta, incidence_norm):
    """Eq. (48): heterogeneous labor-income incidence."""
    agg_ratio = (W / P * N) / (W_ss / P_ss * N_ss)
    expo = 1 + zeta * np.log(agg_ratio)
    labor_income = (W / P) * N * e ** expo / incidence_norm
    return labor_income


@sj.simple
def piH_nkpc(piH, W, Z, PH, mu, kappaH, r):
    """Eq. (49): domestic-price Phillips curve residual."""
    rhs = kappaH * (mu * W / (Z * PH) - 1) + piH(+1) / (1 + r)
    piH_res = piH - rhs
    return piH_res


@sj.simple
def piF_nkpc(piF, E, PF, muF, kappaF, r):
    """Eq. (50): import-price Phillips curve residual."""
    rhs = kappaF * (muF * E / PF - 1) + piF(+1) / (1 + r)
    piF_res = piF - rhs
    return piF_res


@sj.simple
def piHstar_nkpc(piH_star, PH, E, PH_star, muH_star, kappaH_star, r):
    """Eq. (51): export-price Phillips curve residual."""
    rhs = kappaH_star * (muH_star * PH / (E * PH_star) - 1) + piH_star(+1) / (1 + r)
    piH_star_res = piH_star - rhs
    return piH_star_res


@sj.simple
def xH_target_law(xH_hat, PH, P, eta, beta, theta_sub):
    """Eq. (52): delayed-substitution law for domestic target ratio."""
    xH_res = (
        np.log(xH_hat)
        + (1 - beta * theta_sub) * eta * np.log(PH / P)
        - beta * theta_sub * np.log(xH_hat(+1))
    )
    return xH_res


@sj.simple
def xHstar_target_law(xH_star_hat, PH_star, gamma, beta_star, theta_sub):
    """Eq. (53): delayed-substitution law for foreign target ratio."""
    xH_star_res = (
        np.log(xH_star_hat)
        + (1 - beta_star * theta_sub) * gamma * np.log(PH_star)
        - beta_star * theta_sub * np.log(xH_star_hat(+1))
    )
    return xH_star_res


# ---------------------------------------------------------------------------
# Model assembly helper
# ---------------------------------------------------------------------------


def build_ha_oe_model(
    hh_block,
    use_stone_geary=False,
    use_taylor_inertia=False,
    use_quant_extensions=False,
):
    """Return an SSJ model object with blocks aligned to the equation map.

    Parameters
    ----------
    hh_block:
        Household block, typically `@sj.het` with outputs including C and A.
    use_stone_geary:
        If True, uses Eq. (47) aggregator; else Eq. (3).
    use_taylor_inertia:
        If True, uses Eq. (20); else Eq. (19).
    use_quant_extensions:
        If True, appends Eqs. (49)-(53) extension blocks.
    """
    cons_block = consumption_aggregator_stone_geary if use_stone_geary else consumption_aggregator
    mp_block = monetary_rule_taylor if use_taylor_inertia else monetary_rule_const_r

    blocks = [
        cons_block,
        price_index,
        hh_goods_split,
        real_exchange_rate,
        nominal_asset_pricing,
        real_asset_pricing,
        foreign_policy,
        law_of_one_price,
        foreign_demand_home_goods,
        production,
        home_price_flex,
        dividends,
        wage_inflation,
        union_wage_nkpc,
        mp_block,
        hh_block,
        nfa_accounting,
        current_account,
        goods_market_clearing,
    ]

    if use_quant_extensions:
        blocks += [piH_nkpc, piF_nkpc, piHstar_nkpc, xH_target_law, xHstar_target_law]

    return sj.create_model(blocks, name="HA Open Economy")

