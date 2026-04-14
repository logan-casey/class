#!/usr/bin/env python3
"""Translate Sections (1) and (3) of demand_estimates_final.do into Python.

Dependencies:
    pip install pandas numpy pyreadstat statsmodels linearmodels pyblp
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import pyreadstat
import pyblp
import statsmodels.formula.api as smf
from linearmodels.iv import IV2SLS


INTER_CERTNOS = [
    3510,
    628,
    33869,
    3511,
    6548,
    7213,
    6384,
    57957,
    867,
    18409,
    6672,
    12368,
    57890,
    9846,
    588,
    4297,
    6557,
    17534,
    29950,
]

BRAND_CERTNOS = [
    628,
    867,
    3510,
    3511,
    6384,
    6548,
    6672,
    7213,
    9846,
    12368,
    17534,
    18409,
    29950,
    33869,
    57890,
    57957,
]


@dataclass
class Section3Result:
    coeff_unins: pd.Series
    coeff_ins: pd.Series
    table_unins: list[IV2SLS]
    table_ins: list[IV2SLS]
    pyblp_unins: object | None
    pyblp_ins: object | None


def read_dta(path: Path) -> pd.DataFrame:
    df, _ = pyreadstat.read_dta(path)
    return df


def stata_month_to_datetime(month_value: pd.Series) -> pd.Series:
    """Convert Stata monthly date (months since 1960-01-01) to datetime."""
    base = pd.Timestamp("1960-01-01")
    out = pd.to_datetime(
        [base + pd.DateOffset(months=int(m)) if pd.notna(m) else pd.NaT for m in month_value]
    )
    return pd.Series(out, index=month_value.index)


def stata_daily_to_datetime(day_value: pd.Series) -> pd.Series:
    """Convert Stata daily date (days since 1960-01-01) to datetime."""
    return pd.to_datetime(day_value, unit="D", origin="1960-01-01", errors="coerce")


def quarter_index_from_daily(day_value: pd.Series) -> pd.Series:
    dt = stata_daily_to_datetime(day_value)
    return ((dt.dt.year - 1960) * 4 + (dt.dt.quarter - 1)).astype("Int64")


def _wls_predict(df: pd.DataFrame, y: str, interactions: Iterable[str], weight_col: str, out_col: str) -> None:
    rhs = " + ".join(list(interactions) + ["C(q)", "C(certno)", "offdom"])
    formula = f"{y} ~ {rhs}"
    model = smf.wls(formula=formula, data=df, weights=np.sqrt(df[weight_col])).fit(cov_type="HC1")
    df[out_col] = model.predict(df)


def _iv_fit(
    df: pd.DataFrame,
    dep: str,
    endog: list[str],
    instruments: list[str],
    weight_col: str,
) -> IV2SLS:
    exog = "offdom + C(certno) + C(q)"
    endog_block = " + ".join(endog)
    instr_block = " + ".join(instruments)
    formula = f"{dep} ~ 0 + {exog} + [{endog_block} ~ {instr_block}]"
    fit = IV2SLS.from_formula(formula, data=df, weights=np.sqrt(df[weight_col])).fit(cov_type="robust")
    return fit


def _run_pyblp_linear(
    df: pd.DataFrame,
    share_col: str,
    rate_col: str,
    inst_cols: list[str],
    market_col: str,
) -> object | None:
    """Optional pyblp linear logit run for the same demand equation."""
    work = df[[market_col, "certno", share_col, rate_col, "hazard", "offdom", *inst_cols]].copy()
    work = work.rename(columns={market_col: "market_ids", "certno": "firm_ids", share_col: "shares", rate_col: "prices"})
    for i, col in enumerate(inst_cols):
        work[f"demand_instruments{i}"] = work[col]

    formulation = pyblp.Formulation("0 + prices + hazard + offdom", absorb="C(firm_ids) + C(market_ids)")
    try:
        problem = pyblp.Problem(formulation, work)
        return problem.solve()
    except Exception:
        return None


def build_section1_data(
    deposit_data: Path,
    monthly_cds_data: Path,
    monthly_swaps_data: Path,
    cd_rate_data: Path,
    cds_to_hazard_data: Path,
) -> pd.DataFrame:
    df = read_dta(deposit_data)
    df["eq_ratio"] = df["eq"] / df["dep"]

    monthly_cds = read_dta(monthly_cds_data)
    df = df.merge(monthly_cds, on=["shortname", "month"], how="left", validate="m:1")

    swaps = read_dta(monthly_swaps_data)
    df["date2"] = df["date"]
    df["date"] = stata_month_to_datetime(df["month"])
    swaps["date"] = stata_daily_to_datetime(swaps["date"])
    df = df.merge(swaps, on="date", how="left", validate="m:1")
    df = df.drop(columns=["date"]).rename(columns={"date2": "date"})

    cd = read_dta(cd_rate_data)
    df["charter_nbr"] = df["certno"]
    df = df.merge(cd, on=["charter_nbr", "month"], how="left", validate="1:1")

    df["cds"] = df["spread5y"].round(4)
    hazard = read_dta(cds_to_hazard_data)
    df = df.merge(hazard, on="cds", how="inner", validate="m:1")

    df["d_share_unins"] = np.log(df["share_unins"]) - np.log(df["o_unins_share"])
    df["d_share_ins"] = np.log(df["share_ins"]) - np.log(df["o_ins_share"])
    df["loan_exp"] = df["lnlsnet"] / df["asset"]
    df["iv_ln1"] = df["ntre"] / df["asset"]
    df["iv_ln2"] = df["ntlnls"] / df["asset"] - df["ntre"] / df["asset"]
    df["iv_cmo"] = df["sccol"] / df["sc"]
    df["iv_cmo2"] = df["idsccmo"] / df["sc"]

    for cert in INTER_CERTNOS:
        df[f"inter_{cert}"] = df["cmt1y"] * (df["certno"] == cert).astype(float)

    df["q"] = quarter_index_from_daily(df["date"])
    df["ins_rt_spd"] = df["ins_rate"] - df["cmt1y"]
    df["unins_rt_spd"] = df["unins_rate"] - df["cmt1y"]
    return df


def run_section3(df: pd.DataFrame) -> Section3Result:
    interactions = [f"inter_{c}" for c in INTER_CERTNOS]

    # 3.a Uninsured deposits
    df["rt"] = df["unins_rate"] - df["cmt1y"]
    _wls_predict(df, "unins_rate", interactions, "t_depunins", "p_iv_unins_fe")

    unins1 = _iv_fit(df, "d_share_unins", ["rt", "hazard"], ["p_iv_unins_fe", "iv_cmo"], "t_depunins")
    unins2 = _iv_fit(df, "d_share_unins", ["rt", "hazard"], ["p_iv_unins_fe", "iv_ln1", "iv_ln2"], "t_depunins")
    unins3 = _iv_fit(df, "d_share_unins", ["rt", "hazard"], ["p_iv_unins_fe", "iv_ln1", "iv_ln2", "iv_cmo"], "t_depunins")
    df["unins_error"] = unins3.resids
    coeff = unins3.params.copy()

    # 3.b Insured deposits
    df["rt"] = df["ins_rate"] - df["cmt1y"]
    _wls_predict(df, "ins_rate", interactions, "t_depins", "p_iv_ins_fe")

    ins1 = _iv_fit(df, "d_share_ins", ["rt"], ["p_iv_ins_fe"], "t_depins")
    df["ins_error"] = ins1.resids
    coeff2 = ins1.params.copy()

    ins2 = _iv_fit(df, "d_share_ins", ["rt", "hazard"], ["p_iv_ins_fe", "iv_cmo"], "t_depins")
    ins3 = _iv_fit(df, "d_share_ins", ["rt", "hazard"], ["p_iv_ins_fe", "iv_ln1", "iv_ln2"], "t_depins")
    ins4 = _iv_fit(df, "d_share_ins", ["rt", "hazard"], ["p_iv_ins_fe", "iv_ln1", "iv_ln2", "iv_cmo"], "t_depins")

    # 3.c Store demand parameters
    q_terms_unins = [k for k in coeff.index if k.startswith("C(q)")]
    q_terms_ins = [k for k in coeff2.index if k.startswith("C(q)")]

    df["date_u"] = coeff[q_terms_unins[-1]] if q_terms_unins else np.nan
    df["delta_u"] = 0.0
    for cert in BRAND_CERTNOS:
        key = f"C(certno)[{cert}]"
        key_alt = f"C(certno)[{float(cert)}]"
        val = coeff[key] if key in coeff else coeff.get(key_alt, 0.0)
        df.loc[df["certno"] == cert, "delta_u"] = float(val)
    df["office_u"] = float(coeff.get("offdom", np.nan))

    df["date_i"] = coeff2[q_terms_ins[-1]] if q_terms_ins else np.nan
    df["delta_i"] = 0.0
    for cert in BRAND_CERTNOS:
        key = f"C(certno)[{cert}]"
        key_alt = f"C(certno)[{float(cert)}]"
        val = coeff2[key] if key in coeff2 else coeff2.get(key_alt, 0.0)
        df.loc[df["certno"] == cert, "delta_i"] = float(val)
    df["office_i"] = float(coeff2.get("offdom", np.nan))

    df["alpha_unins"] = float(coeff.get("rt", np.nan))
    df["alpha_ins"] = float(coeff2.get("rt", np.nan))
    df["gamma"] = float(coeff.get("hazard", np.nan))

    pyblp_unins = _run_pyblp_linear(
        df=df,
        share_col="share_unins",
        rate_col="unins_rt_spd",
        inst_cols=["p_iv_unins_fe", "iv_ln1", "iv_ln2", "iv_cmo"],
        market_col="month",
    )
    pyblp_ins = _run_pyblp_linear(
        df=df,
        share_col="share_ins",
        rate_col="ins_rt_spd",
        inst_cols=["p_iv_ins_fe", "iv_ln1", "iv_ln2", "iv_cmo"],
        market_col="month",
    )

    return Section3Result(
        coeff_unins=coeff,
        coeff_ins=coeff2,
        table_unins=[unins1, unins2, unins3],
        table_ins=[ins1, ins2, ins3, ins4],
        pyblp_unins=pyblp_unins,
        pyblp_ins=pyblp_ins,
    )


def _print_table_summaries(results: Section3Result) -> None:
    print("\nUninsured 2SLS models (Section 3.a):")
    for i, fit in enumerate(results.table_unins, start=1):
        print(f"  Model {i}: nobs={fit.nobs}, r2={fit.rsquared:.4f}")

    print("\nInsured 2SLS models (Section 3.b):")
    for i, fit in enumerate(results.table_ins, start=1):
        print(f"  Model {i}: nobs={fit.nobs}, r2={fit.rsquared:.4f}")

    if results.pyblp_unins is not None:
        print("\npyblp uninsured model solved.")
    else:
        print("\npyblp uninsured model did not solve (check market/share structure).")

    if results.pyblp_ins is not None:
        print("pyblp insured model solved.")
    else:
        print("pyblp insured model did not solve (check market/share structure).")


def parse_args() -> argparse.Namespace:
    root_default = Path("Data-and-Programs")
    data_default = root_default / "Data-Sets"

    parser = argparse.ArgumentParser(description="Sections 1 and 3 translation for demand_estimates_final.do")
    parser.add_argument("--deposit-data", type=Path, default=data_default / "deposit_data_final.dta")
    parser.add_argument("--monthly-cds-data", type=Path, default=data_default / "monthly cds.dta")
    parser.add_argument("--monthly-swaps-data", type=Path, default=data_default / "monthly_swaps_and_treasuries_rates.dta")
    parser.add_argument("--cd-rate-data", type=Path, default=data_default / "depost_rate_data_1yr_cd_wide.dta")
    parser.add_argument("--cds-to-hazard-data", type=Path, default=data_default / "cds-to-hazard.dta")
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("Data-and-Programs/Demand-Estimates/demand_estimates_sections1_3_output.csv"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    df = build_section1_data(
        deposit_data=args.deposit_data,
        monthly_cds_data=args.monthly_cds_data,
        monthly_swaps_data=args.monthly_swaps_data,
        cd_rate_data=args.cd_rate_data,
        cds_to_hazard_data=args.cds_to_hazard_data,
    )
    results = run_section3(df)
    _print_table_summaries(results)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)
    print(f"\nSaved translated Section 1 + 3 dataset to: {args.output_csv}")


if __name__ == "__main__":
    main()
