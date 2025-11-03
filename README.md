# KinaseSeek: A Case Study in Virtual Screening for Kinase Inhibitors

1) Objective: This project is a self-directed case study demonstrating a complete structure-based virtual screening workflow to identify potential inhibitors for a specific target.
The goal is to showcase the full process: from target and library selection to docking, post-processing, and hit analysis. The entire workflow is documented in a single Jupyter Notebook.

2) Key Features:
- Target selection and preparation.
- Ligand library preparation.
- Virtual screening.
- Hit analysis and visualisation.

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


