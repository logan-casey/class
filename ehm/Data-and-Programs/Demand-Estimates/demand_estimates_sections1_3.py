#!/usr/bin/env python3
"""Translate Sections (1) and (3) of demand_estimates_final.do into Python.

This version is CSV-native (uses exported .csv files instead of .dta files).
Dependencies:
    pip install pandas numpy statsmodels linearmodels pyblp
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
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


def read_csv_data(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def parse_stata_month(series: pd.Series) -> pd.Series:
    """Parse Stata monthly dates from either int months or strings like 2002m3."""
    if pd.api.types.is_numeric_dtype(series):
        return pd.to_numeric(series, errors="coerce").round().astype("Int64")

    s = series.astype(str).str.strip()
    m = s.str.extract(r"^(?P<year>\d{4})m(?P<mon>\d{1,2})$", expand=True)
    out = pd.Series(pd.NA, index=series.index, dtype="Int64")

    matched = m["year"].notna()
    if matched.any():
        y = pd.to_numeric(m.loc[matched, "year"], errors="coerce")
        mon = pd.to_numeric(m.loc[matched, "mon"], errors="coerce")
        out.loc[matched] = ((y - 1960) * 12 + (mon - 1)).astype("Int64")

    if (~matched).any():
        dt = pd.to_datetime(s.loc[~matched], errors="coerce")
        out.loc[~matched] = ((dt.dt.year - 1960) * 12 + (dt.dt.month - 1)).astype("Int64")

    return out


def month_to_datetime(stata_month: pd.Series) -> pd.Series:
    base = pd.Timestamp("1960-01-01")
    vals = pd.to_numeric(stata_month, errors="coerce")
    out = [
        (base + pd.DateOffset(months=int(v))) if pd.notna(v) else pd.NaT
        for v in vals
    ]
    return pd.Series(pd.to_datetime(out), index=stata_month.index)


def parse_stata_daily(series: pd.Series) -> pd.Series:
    """Parse Stata daily dates from int days or strings like 31mar2002."""
    if isinstance(series, pd.DataFrame):
        series = series.iloc[:, 0]
    if pd.api.types.is_datetime64_any_dtype(series):
        return pd.to_datetime(series, errors="coerce")
    if pd.api.types.is_numeric_dtype(series):
        return pd.to_datetime(series, unit="D", origin="1960-01-01", errors="coerce")

    s = series.astype(str).str.strip()
    d = pd.to_datetime(s, format="%d%b%Y", errors="coerce")
    iso = pd.to_datetime(s, format="%Y-%m-%d", errors="coerce")
    return d.fillna(iso)


def quarter_index_from_daily(day_value: pd.Series) -> pd.Series:
    dt = parse_stata_daily(day_value)
    return ((dt.dt.year - 1960) * 4 + (dt.dt.quarter - 1)).astype("Int64")


def _wls_predict(df: pd.DataFrame, y: str, interactions: Iterable[str], weight_col: str, out_col: str) -> None:
    rhs = " + ".join(list(interactions) + ["C(q)", "C(certno)", "offdom"])
    model = smf.wls(formula=f"{y} ~ {rhs}", data=df, weights=np.sqrt(df[weight_col])).fit(cov_type="HC1")
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
    return IV2SLS.from_formula(formula, data=df, weights=np.sqrt(df[weight_col])).fit(cov_type="robust")


def _run_pyblp_linear(
    df: pd.DataFrame,
    share_col: str,
    rate_col: str,
    inst_cols: list[str],
    market_col: str,
) -> object | None:
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
    cd_rate_data: Path | None,
    cds_to_hazard_data: Path,
) -> pd.DataFrame:
    df = read_csv_data(deposit_data)
    df["month_key"] = parse_stata_month(df["month"])
    df["eq_ratio"] = df["eq"] / df["dep"]

    monthly_cds = read_csv_data(monthly_cds_data)
    monthly_cds["month_key"] = parse_stata_month(monthly_cds["month"])
    monthly_cds = monthly_cds[["shortname", "month_key", "spread1y", "spread5y"]]
    df = df.merge(monthly_cds, on=["shortname", "month_key"], how="left", validate="m:1")

    swaps = read_csv_data(monthly_swaps_data)
    swaps["swap_date"] = parse_stata_daily(swaps["date"]).dt.normalize()
    swaps = swaps.drop(columns=["date"])
    df["date_merge"] = month_to_datetime(df["month_key"]).dt.normalize()
    df = df.merge(swaps, left_on="date_merge", right_on="swap_date", how="left", validate="m:1")
    df = df.drop(columns=["swap_date", "date_merge"])

    if cd_rate_data is not None and cd_rate_data.exists():
        cd = read_csv_data(cd_rate_data)
        if "month" in cd.columns:
            cd["month_key"] = parse_stata_month(cd["month"])
        if "charter_nbr" not in cd.columns and "certno" in cd.columns:
            cd["charter_nbr"] = pd.to_numeric(cd["certno"], errors="coerce")
        df["charter_nbr"] = pd.to_numeric(df["certno"], errors="coerce")
        df = df.merge(cd, on=["charter_nbr", "month_key"], how="left", validate="1:1")

    # Accept legacy wide names from CD-rate file.
    if "ins_rate" not in df.columns and "rate1" in df.columns:
        df["ins_rate"] = df["rate1"]
    if "unins_rate" not in df.columns and "rate2" in df.columns:
        df["unins_rate"] = df["rate2"]

    required_now = ["ins_rate", "unins_rate", "cmt1y"]
    missing_now = [c for c in required_now if c not in df.columns]
    if missing_now:
        rate_like = [c for c in df.columns if "rate" in c.lower()]
        raise ValueError(
            "Missing required columns after Section 1 merges: "
            f"{missing_now}. Current rate-like columns: {rate_like[:25]}. "
            "If needed, inspect your CD-rate file and pass it with --cd-rate-data. "
            "Expected columns include ins_rate/unins_rate or rate1/rate2, plus keys charter_nbr and month."
        )

    df["cds"] = df["spread5y"].round(4)
    hazard = read_csv_data(cds_to_hazard_data)
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
    df["month"] = df["month_key"]
    return df


def run_section3(df: pd.DataFrame) -> Section3Result:
    interactions = [f"inter_{c}" for c in INTER_CERTNOS]

    df["rt"] = df["unins_rate"] - df["cmt1y"]
    _wls_predict(df, "unins_rate", interactions, "t_depunins", "p_iv_unins_fe")
    unins1 = _iv_fit(df, "d_share_unins", ["rt", "hazard"], ["p_iv_unins_fe", "iv_cmo"], "t_depunins")
    unins2 = _iv_fit(df, "d_share_unins", ["rt", "hazard"], ["p_iv_unins_fe", "iv_ln1", "iv_ln2"], "t_depunins")
    unins3 = _iv_fit(df, "d_share_unins", ["rt", "hazard"], ["p_iv_unins_fe", "iv_ln1", "iv_ln2", "iv_cmo"], "t_depunins")
    df["unins_error"] = unins3.resids
    coeff = unins3.params.copy()

    df["rt"] = df["ins_rate"] - df["cmt1y"]
    _wls_predict(df, "ins_rate", interactions, "t_depins", "p_iv_ins_fe")
    ins1 = _iv_fit(df, "d_share_ins", ["rt"], ["p_iv_ins_fe"], "t_depins")
    df["ins_error"] = ins1.resids
    coeff2 = ins1.params.copy()
    ins2 = _iv_fit(df, "d_share_ins", ["rt", "hazard"], ["p_iv_ins_fe", "iv_cmo"], "t_depins")
    ins3 = _iv_fit(df, "d_share_ins", ["rt", "hazard"], ["p_iv_ins_fe", "iv_ln1", "iv_ln2"], "t_depins")
    ins4 = _iv_fit(df, "d_share_ins", ["rt", "hazard"], ["p_iv_ins_fe", "iv_ln1", "iv_ln2", "iv_cmo"], "t_depins")

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
    print("\npyblp uninsured model solved." if results.pyblp_unins is not None else "\npyblp uninsured model did not solve.")
    print("pyblp insured model solved." if results.pyblp_ins is not None else "pyblp insured model did not solve.")


def parse_args() -> argparse.Namespace:
    data_default = Path("Data-and-Programs/Data-Sets")
    parser = argparse.ArgumentParser(description="Sections 1 and 3 translation for demand_estimates_final.do (CSV inputs)")
    parser.add_argument("--deposit-data", type=Path, default=data_default / "deposit_data_final.csv")
    parser.add_argument("--monthly-cds-data", type=Path, default=data_default / "monthly cds.csv")
    parser.add_argument("--monthly-swaps-data", type=Path, default=data_default / "monthly_swaps_and_treasuries_rates.csv")
    parser.add_argument("--cd-rate-data", type=Path, default=Path("Data-and-Programs/Data-Sets/depost_rate_data_1yr_cd_wide.csv"))
    parser.add_argument("--cds-to-hazard-data", type=Path, default=data_default / "cds-to-hazard.csv")
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
