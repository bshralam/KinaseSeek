#!/usr/bin/env python
"""Download EGFR ligands from ChEMBL via REST API.

Example:
    python scripts/01_download_chembl.py --out data/egfr_chembl.csv --min_pchembl 6.0
"""
import argparse, sys, time
import pandas as pd
import requests
from tqdm import tqdm

EGFR_TARGET = "CHEMBL203"  # Human EGFR

def fetch_activities(min_pchembl: float = 0.0, max_per_mol: int = 1, limit: int = 5000) -> pd.DataFrame:
    base = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
    params = {
        "target_chembl_id": EGFR_TARGET,
        "standard_type": "IC50",
        "limit": 1000,
    }
    rows = []
    fetched = 0
    seen_mols = {}

    with tqdm(total=limit, desc="Downloading", unit="rec") as pbar:
        url = base
        while url and fetched < limit:
            r = requests.get(url, params=params if url == base else None, timeout=60)
            if r.status_code != 200:
                print(f"[WARN] HTTP {r.status_code}: {r.text[:200]}", file=sys.stderr)
                break
            data = r.json()
            for rec in data.get("activities", []):
                smi = rec.get("canonical_smiles") or rec.get("smiles")
                if not smi:
                    continue
                pchembl = rec.get("pchembl_value")
                try:
                    pchembl = float(pchembl) if pchembl is not None else None
                except:
                    pchembl = None
                if pchembl is not None and pchembl < min_pchembl:
                    continue

                mol_id = rec.get("molecule_chembl_id")
                if mol_id:
                    seen_mols.setdefault(mol_id, 0)
                    if seen_mols[mol_id] >= max_per_mol:
                        continue
                    seen_mols[mol_id] += 1

                rows.append({
                    "molecule_chembl_id": mol_id,
                    "canonical_smiles": smi,
                    "pchembl_value": pchembl,
                    "standard_value": rec.get("standard_value"),
                    "standard_units": rec.get("standard_units"),
                    "assay_chembl_id": rec.get("assay_chembl_id"),
                    "record_id": rec.get("record_id"),
                })
                fetched += 1
                pbar.update(1)
                if fetched >= limit:
                    break

            url = data.get("page_meta", {}).get("next", None)
            if url is None:
                break
            time.sleep(0.1)

    df = pd.DataFrame(rows).drop_duplicates()
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument("--min_pchembl", type=float, default=0.0, help="Keep rows with pChEMBL >= this")
    ap.add_argument("--max_per_mol", type=int, default=1, help="Max records per unique ChEMBL molecule")
    ap.add_argument("--limit", type=int, default=5000, help="Max records to fetch")
    args = ap.parse_args()

    df = fetch_activities(min_pchembl=args.min_pchembl, max_per_mol=args.max_per_mol, limit=args.limit)
    if df.empty:
        print("[WARN] No records fetched. Try relaxing filters.", file=sys.stderr)
    df.to_csv(args.out, index=False)
    print(f"[OK] Saved {len(df)} records to {args.out}")

if __name__ == "__main__":
    main()
