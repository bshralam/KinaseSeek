#!/usr/bin/env python
import pandas as pd, matplotlib.pyplot as plt

zinc = pd.read_csv("data/top_ranked_hits_prod.csv")
osim = pd.read_csv("results/top_hits_osimertinib.csv")

plt.figure(figsize=(6,4))
plt.hist(zinc["avg_tanimoto"], bins=20, alpha=0.6, label="EGFR vs ZINC (avg)")
plt.hist(osim["tanimoto_score"], bins=20, alpha=0.6, label="EGFR vs Osimertinib (max)")
plt.xlabel("Tanimoto similarity")
plt.ylabel("Count")
plt.legend()
plt.title("Comparison of Global vs Local Similarity Distributions")
plt.tight_layout()
plt.savefig("results/similarity_comparison.png", dpi=200)
plt.show()
