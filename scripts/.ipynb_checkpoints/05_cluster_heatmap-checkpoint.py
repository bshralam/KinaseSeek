#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list

def smiles_to_fps(smiles_list, radius=2, nBits=2048):
    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
            fps.append(fp)
    return fps

def tanimoto_matrix(fps):
    n = len(fps)
    sim = np.zeros((n, n))
    for i in range(n):
        sim[i, i] = 1.0
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i+1:])
        sim[i, i+1:] = sims
        sim[i+1:, i] = sims
    return sim

def cluster_and_plot(sim, names, out_png, title):
    dist = 1 - sim
    Z = linkage(dist, method="ward")
    order = leaves_list(Z)
    sim_sorted = sim[order][:, order]
    plt.figure(figsize=(8,8))
    plt.imshow(sim_sorted, cmap="viridis", interpolation="nearest")
    plt.colorbar(label="Tanimoto similarity")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"[OK] Saved heatmap -> {out_png}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_csv", required=True)
    ap.add_argument("--smiles_col", default="smiles")
    ap.add_argument("--max_n", type=int, default=100)
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--title", default="Clustered Tanimoto Heatmap")
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv).drop_duplicates(subset=[args.smiles_col]).dropna(subset=[args.smiles_col])
    df = df.head(args.max_n)
    fps = smiles_to_fps(df[args.smiles_col].tolist())
    sim = tanimoto_matrix(fps)
    cluster_and_plot(sim, df[args.smiles_col].tolist(), args.out_png, args.title)

if __name__ == "__main__":
    main()

