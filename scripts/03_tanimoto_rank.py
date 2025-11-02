#!/usr/bin/env python
"""Rank compounds by Tanimoto similarity to a query (or set of queries)."""
import argparse
import numpy as np
import pandas as pd
from rdkit import DataStructs
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from utils.chem import mol_from_smiles, morgan_fp

def bitvect_from_text(text: str):
    return DataStructs.CreateFromBitString(text)

def fingerprints_from_parquet(parquet_path: str):
    df = pd.read_parquet(parquet_path)
    fps = [bitvect_from_text(t) for t in df["fp_text"].tolist()]
    return df, fps

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fps", required=True, help="Parquet of fingerprints (from 02)")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--query", help="Single query SMILES")
    group.add_argument("--query_csv", help="CSV with query SMILES")
    ap.add_argument("--query_smiles_col", default="smiles", help="Column name for queries if --query_csv is used")
    ap.add_argument("--topk", type=int, default=50)
    ap.add_argument("--out", default="results/top_hits.csv")
    ap.add_argument("--radius", type=int, default=2)
    ap.add_argument("--nBits", type=int, default=2048)
    args = ap.parse_args()

    lib_df, lib_fps = fingerprints_from_parquet(args.fps)

    # build query fps
    query_fps = []
    if args.query:
        mol = mol_from_smiles(args.query)
        if mol is None:
            raise SystemExit("Invalid query SMILES.")
        query_fps.append(morgan_fp(mol, radius=args.radius, nBits=args.nBits))
    else:
        qdf = pd.read_csv(args.query_csv)
        if args.query_smiles_col not in qdf.columns:
            raise SystemExit(f"Query SMILES column '{args.query_smiles_col}' missing in {args.query_csv}")
        for s in qdf[args.query_smiles_col].astype(str).tolist():
            mol = mol_from_smiles(s)
            if mol is None:
                continue
            query_fps.append(morgan_fp(mol, radius=args.radius, nBits=args.nBits))
        if not query_fps:
            raise SystemExit("No valid queries after parsing CSV.")

    # score = max Tanimoto over all queries
    scores = []
    for fp in lib_fps:
        sims = [DataStructs.TanimotoSimilarity(fp, q) for q in query_fps]
        scores.append(max(sims) if sims else 0.0)

    out = lib_df.copy()
    out["tanimoto_score"] = scores
    out = out.sort_values("tanimoto_score", ascending=False).head(args.topk)
    out.to_csv(args.out, index=False)
    print(f"[OK] Saved top {len(out)} hits -> {args.out}")

if __name__ == "__main__":
    main()
