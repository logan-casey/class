# HA-OE Equation List and SSJ Block Map

This document maps equations (1)-(23) and (47)-(53) from `ha_oe.pdf` to a `sequence_jacobian` block design, using naming conventions close to `ssj_ref` and `repkit/main_sec67.ipynb`.

## 1) Necessary Equations

### Core model equations (1)-(23)

| Eq. | Equation role | Suggested code target |
|---|---|---|
| (1) | Full household problem with portfolio choice over equity, home bonds, foreign bonds + borrowing constraint | `hh_portfolio_accounting` (accounting block; mostly reduced away by (8)-(9)) |
| (2) | Period utility: CRRA over CES consumption basket, disutility of labor | `hh_utility_params` (parameters used by HH and union blocks) |
| (3) | CES aggregator over foreign and home goods | `consumption_aggregator` (`@sj.simple`) |
| (4) | CPI index from import/home prices | `price_index` (`@sj.simple`) |
| (5) | Nominal no-arbitrage: UIP and stock pricing in nominal terms | `nominal_asset_pricing` (`@sj.simple`, residuals) |
| (6) | Real exchange rate definition \(Q_t = E_t P_t^* / P_t\) | `real_exchange_rate` (`@sj.simple`) |
| (7) | Real UIP and real stock pricing | `real_asset_pricing` (`@sj.simple`, residuals) |
| (8) | Definitions of beginning/end-of-period real assets \(a_t^p, a_t'\) | `asset_transition_defs` (`@sj.simple`) |
| (9) | Canonical HH dynamic problem in \((a_t^p, e)\): budget + borrowing constraint | `hh_oe` (`@sj.het`) |
| (10) | Import demand from intratemporal CES FOC | `hh_goods_split` (`@sj.simple`) |
| (11) | Home-good demand from intratemporal CES FOC | `hh_goods_split` (`@sj.simple`) |
| (12) | Foreign demand for home goods | `foreign_demand_home_goods` (`@sj.simple`) |
| (13) | Production technology \(Y_t = N_t\) | `production` (`@sj.simple`) |
| (14) | Flexible-price markup pricing \(P_{H,t} = \mu W_t\) | `home_price_flex` (`@sj.simple`) |
| (15) | Real dividends including unhedged FX export term | `dividends` (`@sj.simple`) |
| (16) | Foreign interest rate process from discount shifter \(B_t\) | `foreign_policy` (`@sj.simple`) |
| (17) | Law of one price for exports \(P^*_{H,t} = P_{H,t}/E_t\) | `law_of_one_price` (`@sj.simple`) |
| (18) | Wage Phillips curve (Calvo wages) | `union_wage_nkpc` (`@sj.simple`, residual) |
| (19) | Constant-real-rate monetary rule | `monetary_rule_const_r` (`@sj.simple`) |
| (20) | Alternative Taylor rule with inertia | `monetary_rule_taylor` (`@sj.simple`) |
| (21) | Home goods market clearing \(C_{H,t}+C^*_{H,t}=Y_t\) | `goods_market_clearing` (`@sj.simple`, residual) |
| (22) | Net foreign assets definition \(nfa_t=A_t-p_t\) | `nfa_accounting` (`@sj.simple`) |
| (23) | Current account identity / NFA law of motion | `current_account` (`@sj.simple`, residual) |

### Quantitative extensions (47)-(53)

| Eq. | Equation role | Suggested code target |
|---|---|---|
| (47) | Stone-Geary CES basket with import subsistence term | `consumption_aggregator_stone_geary` (`@sj.simple`) |
| (48) | Heterogeneous labor-income incidence elasticity to aggregate wage bill | `income_incidence` (`hetinput` used by `hh_oe`) |
| (49) | Domestic-price Phillips curve (Calvo \(P_H\)) | `piH_nkpc` (`@sj.simple`, residual) |
| (50) | Import-price Phillips curve (Calvo \(P_F\), incomplete pass-through) | `piF_nkpc` (`@sj.simple`, residual) |
| (51) | Export-price Phillips curve (Calvo \(P_H^*\)) | `piHstar_nkpc` (`@sj.simple`, residual) |
| (52) | Delayed substitution dynamics for domestic target ratio \(\hat x_{H,t}\) | `xH_target_law` (`@sj.simple` or `@sj.solved`) |
| (53) | Delayed substitution dynamics for foreign target ratio \(\hat x^*_{H,t}\) | `xHstar_target_law` (`@sj.simple` or `@sj.solved`) |

## 2) Block-Level Mapping (Sequence-Jacobian)

### A. Household side

| Block name | Type | Main equations |
|---|---|---|
| `hh_oe` | `@sj.het` | (2), (9) |
| `income_incidence` | `hetinput` | (48) |
| `hh_goods_split` | `@sj.simple` | (10), (11) |
| `consumption_aggregator` or `consumption_aggregator_stone_geary` | `@sj.simple` | (3) or (47) |
| `asset_transition_defs` | `@sj.simple` | (8), with optional accounting to (1) |

### B. Prices, returns, and policy

| Block name | Type | Main equations |
|---|---|---|
| `price_index` | `@sj.simple` | (4) |
| `real_exchange_rate` | `@sj.simple` | (6) |
| `nominal_asset_pricing` | `@sj.simple` (residuals) | (5) |
| `real_asset_pricing` | `@sj.simple` (residuals) | (7) |
| `foreign_policy` | `@sj.simple` | (16) |
| `law_of_one_price` | `@sj.simple` | (17) |
| `union_wage_nkpc` | `@sj.simple` (residual) | (18) |
| `monetary_rule_const_r` or `monetary_rule_taylor` | `@sj.simple` | (19) or (20) |
| `piH_nkpc`, `piF_nkpc`, `piHstar_nkpc` | `@sj.simple` (residuals) | (49)-(51) |
| `xH_target_law`, `xHstar_target_law` | `@sj.simple`/`@sj.solved` | (52)-(53) |

### C. Production, external demand, and equilibrium

| Block name | Type | Main equations |
|---|---|---|
| `production` | `@sj.simple` | (13) |
| `home_price_flex` (baseline) | `@sj.simple` | (14) |
| `dividends` | `@sj.simple` | (15) |
| `foreign_demand_home_goods` | `@sj.simple` | (12) |
| `goods_market_clearing` | `@sj.simple` (residual) | (21) |
| `nfa_accounting` | `@sj.simple` | (22) |
| `current_account` | `@sj.simple` (residual) | (23) |

## 3) Suggested DAG Assembly

```python
import sequence_jacobian as sj

# HH core (with optional extension hooks)
hh_ext = hh_oe.add_hetinputs([income_incidence, make_grids_oe])

blocks = [
    # preferences/demand
    consumption_aggregator,      # or consumption_aggregator_stone_geary
    hh_goods_split,

    # prices and returns
    price_index, real_exchange_rate,
    nominal_asset_pricing, real_asset_pricing,
    foreign_policy, law_of_one_price,

    # production and dividends
    production, home_price_flex, dividends,
    foreign_demand_home_goods,

    # policy and wage setting
    union_wage_nkpc,
    monetary_rule_const_r,       # or monetary_rule_taylor

    # extension blocks (optional quantitative section)
    piH_nkpc, piF_nkpc, piHstar_nkpc,
    xH_target_law, xHstar_target_law,

    # household and equilibrium accounting
    hh_ext, asset_transition_defs,
    nfa_accounting, current_account, goods_market_clearing,
]

ha_oe_model = sj.create_model(blocks, name="HA Open Economy")
```

## 4) Practical Notes for Coding

1. In implementation, equation (9) is the key HH dynamic block; (1) and most of portfolio detail become accounting identities unless you explicitly model endogenous portfolio composition.
2. Treat forward-looking equations as residual blocks and solve with model-level unknowns/targets (as in `main_sec67`).
3. Use equation switches:
   - baseline demand: (3), baseline pricing: (14), policy: (19)
   - quantitative extensions: replace with (47), (49)-(53), and optionally (20)
