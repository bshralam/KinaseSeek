# KinaseSeek: A Case Study in Virtual Screening for Kinase Inhibitors

1) Objective: This project is a self-directed case study demonstrating a ligand-based virtual screening workflow for kinase inhibitors, focused on EGFR, a key oncology target.

2) Key Features:

- Fingerprint generation: Morgan (ECFP, radius = 2) fingerprints computed with RDKit for all EGFR ligands.
- Two complementary similarity analyses:
  i) Global baseline: Average Tanimoto similarity of EGFR ligands to a 50 K ZINC fragment set, used as a negative control to show overall structural novelty.
  ii) Reference-guided search: Maximum Tanimoto similarity to the known inhibitor Osimertinib, used to find analog clusters and local structure–activity trends.
- Clustering and visualization: Hierarchical clustering and heatmaps to show analog families, and histograms to compare global vs. local similarity patterns.

3) Tools Used: Python, RDKit, pandas, Matplotlib, PyMol.

4) Figure Summary:
   
a) Histogram comparing the Tanimoto similarity distribution between EGFR inhibitors and a ZINC fragment subset (global average similarity) versus Osimertinib (reference-guided similarity).
The blue curve shows chemical diversity, while the orange curve highlights closer analogs related to the known EGFR inhibitor.

<img width="1200" height="800" alt="image" src="https://github.com/user-attachments/assets/b83b31ef-3ca8-4202-a1ac-8fd5626a309a" />

b) Heatmap of pairwise Tanimoto similarity among 300 ZINC molecules.
The mostly dark blue pattern shows very low similarity, serving as a baseline for chemical diversity.

<img width="1600" height="1600" alt="image" src="https://github.com/user-attachments/assets/3a8b8373-1f76-487d-962d-f72bfac4f7d9" />

c) Heatmap of pairwise similarity among the top EGFR ligands similar to Osimertinib.
Bright clusters show groups of related analogs sharing common scaffolds.

<img width="1600" height="1600" alt="image" src="https://github.com/user-attachments/assets/96946581-2b9f-498a-861e-190d72d4d22b" />


