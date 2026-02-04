#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path
from typing import Dict, Optional, Tuple


# ---------------------------
# Helpers
# ---------------------------
def is_missing(x: str) -> bool:
    return (not x) or x in ("NA", "NaN", "nan", ".", "NULL", "null")


def to_int(x: str) -> Optional[int]:
    if is_missing(x):
        return None
    try:
        # allow "123.0"
        if "." in x:
            v = float(x)
            if v != v:  # NaN
                return None
            return int(v)
        return int(x)
    except Exception:
        return None


def trait_from_filename(path: str) -> str:
    base = os.path.basename(path)
    if base.endswith("_meta_out.tsv.gz"):
        return base[: -len("_meta_out.tsv.gz")]
    if base.endswith(".tsv.gz"):
        return base[: -len(".tsv.gz")]
    if base.endswith(".gz"):
        return base[: -len(".gz")]
    return os.path.splitext(base)[0]


def trait_code(trait: str) -> str:
    # finngen_R12_AB1_AMOEBIASIS -> AB1_AMOEBIASIS
    if trait.startswith("finngen_R12_"):
        return trait[len("finngen_R12_"):]
    return trait


# ---------------------------
# Read afreq (rsid -> ALT_FREQS)
# ---------------------------
def load_afreq_map(afreq_path: Path) -> Dict[str, str]:
    """
    Expected columns:
    #CHROM  ID  REF ALT PROVISIONAL_REF? ALT_FREQS OBS_CT
    """
    m: Dict[str, str] = {}
    with afreq_path.open("rt", encoding="utf-8", newline="") as f:
        header = f.readline()
        if not header:
            raise ValueError(f"Empty afreq file: {afreq_path}")

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            rsid = parts[1]
            if is_missing(rsid):
                continue
            freq = parts[5]
            if not is_missing(freq):
                m[rsid] = freq

    if not m:
        raise ValueError(f"No rsid->freq entries found in: {afreq_path}")
    return m


# ---------------------------
# Read mapping file for trait totals
# ---------------------------
def find_col(cols: list[str], name: str) -> int:
    # exact or case-insensitive
    if name in cols:
        return cols.index(name)
    lower = {c.lower(): i for i, c in enumerate(cols)}
    if name.lower() in lower:
        return lower[name.lower()]
    raise ValueError(f"Column '{name}' not found. Columns: {cols}")


def load_trait_totals(mapping_tsv: Path, tcode: str) -> Tuple[int, int]:
    """
    Reads mapping file and returns:
      fg_total = fg_n_cases + fg_n_controls
      ukbb_total = ukbb_n_cases + ukbb_n_controls

    We assume mapping file contains columns named exactly:
      fg_n_cases, fg_n_controls, ukbb_n_cases, ukbb_n_controls
    Trait column is auto-detected among common names; if your file uses a specific
    trait column name, set --trait_col.
    """
    with mapping_tsv.open("rt", encoding="utf-8", newline="") as f:
        header = f.readline()
        if not header:
            raise ValueError(f"Empty mapping file: {mapping_tsv}")
        cols = header.rstrip("\n").split("\t")

        # required N component columns
        i_fgc = find_col(cols, "fg_n_cases")
        i_fgu = find_col(cols, "fg_n_controls")
        i_ukc = find_col(cols, "ukbb_n_cases")
        i_uku = find_col(cols, "ukbb_n_controls")

        # trait column: try common ones
        trait_candidates = ["trait", "endpoint", "phenotype", "pheno", "phenocode", "FINNGEN_ENDPOINT", "finngen_endpoint"]
        i_trait = None
        for cand in trait_candidates:
            try:
                i_trait = find_col(cols, cand)
                break
            except Exception:
                pass
        if i_trait is None:
            # fallback: first column containing 'trait' or 'endpoint' or 'pheno'
            for i, c in enumerate(cols):
                cl = c.lower()
                if ("trait" in cl) or ("endpoint" in cl) or ("pheno" in cl):
                    i_trait = i
                    break
        if i_trait is None:
            raise ValueError(f"Could not detect trait column in mapping file. Columns: {cols}")

        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_trait, i_fgc, i_fgu, i_ukc, i_uku):
                continue

            t = parts[i_trait]
            if is_missing(t):
                continue

            # match: exact code or contains code
            if not (t == tcode or t.endswith(tcode) or tcode in t):
                continue

            fgc = to_int(parts[i_fgc])
            fgu = to_int(parts[i_fgu])
            ukc = to_int(parts[i_ukc])
            uku = to_int(parts[i_uku])

            if fgc is None or fgu is None:
                raise ValueError(f"Mapping file has missing FinnGen n for trait {tcode}: fg_n_cases/controls")

            fg_total = fgc + fgu

            # UKBB might be 0 for FinnGen-only traits; treat missing as 0
            uk_total = 0
            if ukc is not None and uku is not None:
                uk_total = ukc + uku

            return fg_total, uk_total

    raise ValueError(f"Trait {tcode} not found in mapping file: {mapping_tsv}")


# ---------------------------
# Sumstats conversion
# ---------------------------
SUMSTATS_REQUIRED = [
    "REF",
    "ALT",
    "rsid",
    "all_inv_var_meta_beta",
    "all_inv_var_meta_sebeta",
    "all_inv_var_meta_p",
    "all_meta_N",  # tag only (1 or 2)
]


def parse_sumstats_header(header_line: str) -> Dict[str, int]:
    cols = header_line.rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(cols)}
    missing = [c for c in SUMSTATS_REQUIRED if c not in idx]
    if missing:
        raise ValueError(f"Missing required columns in sumstats: {missing}\nHeader={cols}")
    return idx


def convert(sumstats_gz: Path, afreq: Path, mapping_tsv: Path, out_dir: Path) -> Path:
    freq_map = load_afreq_map(afreq)

    trait = trait_from_filename(str(sumstats_gz))
    tcode = trait_code(trait)

    fg_total, uk_total = load_trait_totals(mapping_tsv, tcode)
    N1 = fg_total
    N2 = fg_total + uk_total if uk_total > 0 else fg_total  # still non-missing

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{trait}.ma"

    with gzip.open(sumstats_gz, "rt", encoding="utf-8", newline="") as fin, \
         out_path.open("wt", encoding="utf-8", newline="") as fout:

        header = fin.readline()
        if not header:
            raise ValueError(f"Empty sumstats file: {sumstats_gz}")
        idx = parse_sumstats_header(header)

        i_ref = idx["REF"]
        i_alt = idx["ALT"]
        i_rsid = idx["rsid"]
        i_b = idx["all_inv_var_meta_beta"]
        i_se = idx["all_inv_var_meta_sebeta"]
        i_p = idx["all_inv_var_meta_p"]
        i_tag = idx["all_meta_N"]

        max_need = max(i_ref, i_alt, i_rsid, i_b, i_se, i_p, i_tag)

        fout.write("SNP\tA1\tA2\tfreq\tb\tse\tp\tn\n")

        for line in fin:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max_need:
                continue

            rsid = parts[i_rsid]
            if is_missing(rsid):
                continue

            freq = freq_map.get(rsid)
            if freq is None:
                continue

            b = parts[i_b]
            se = parts[i_se]
            p = parts[i_p]
            if is_missing(b) or is_missing(se) or is_missing(p):
                continue

            ref = parts[i_ref]
            alt = parts[i_alt]
            if is_missing(ref) or is_missing(alt):
                continue

            tag = to_int(parts[i_tag])
            if tag == 1:
                n = N1
            elif tag == 2:
                n = N2
            else:
                # Ensure never missing: default to the largest available sample
                n = N2

            fout.write(f"{rsid}\t{alt}\t{ref}\t{freq}\t{b}\t{se}\t{p}\t{n}\n")

    return out_path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sumstats_gz", required=True, type=Path)
    ap.add_argument("--afreq", required=True, type=Path)
    ap.add_argument("--mapping_tsv", required=True, type=Path)
    ap.add_argument("--out_dir", default=Path("."), type=Path)
    args = ap.parse_args()

    out = convert(args.sumstats_gz, args.afreq, args.mapping_tsv, args.out_dir)
    print(str(out))


if __name__ == "__main__":
    main()
