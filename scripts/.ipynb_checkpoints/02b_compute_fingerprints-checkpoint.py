#!/usr/bin/env python
"""Compute Morgan fingerprints from a CSV of SMILES.

Input CSV must contain a SMILES column (default: canonical_smiles).
Outputs a parquet with columns: [smiles, chembl_id, fp_text, radius, nBits]
"""
import argparse
import pandas as pd
from rdkit import DataStructs
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from utils.chem import mol_from_smiles, morgan_fp

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input CSV with SMILES")
    ap.add_argument("--out", required=True, help="Output parquet file")
    ap.add_argument("--smiles_col", default="canonical_smiles")
    ap.add_argument("--id_col", default="molecule_chembl_id")
    ap.add_argument("--radius", type=int, default=2)
    ap.add_argument("--nBits", type=int, default=2048)
    args = ap.parse_args()

    df = pd.read_csv(args.inp)
    if args.smiles_col not in df.columns:
        raise SystemExit(f"SMILES column '{args.smiles_col}' not found in {args.inp}")

    out_rows = []
    for _, row in df.iterrows():
        smi = str(row[args.smiles_col])
        mid = row[args.id_col] if args.id_col in df.columns else None
        mol = mol_from_smiles(smi)
        if mol is None:
            continue
        fp = morgan_fp(mol, radius=args.radius, nBits=args.nBits)
        out_rows.append({
            "smiles": smi,
            "chembl_id": mid,
            "fp_text": DataStructs.BitVectToText(fp),
            "radius": args.radius,
            "nBits": args.nBits
        })
    out_df = pd.DataFrame(out_rows)
    out_df.to_parquet(args.out, index=False)
    print(f"[OK] Wrote {len(out_df)} fingerprints -> {args.out}")

if __name__ == "__main__":
    main()
