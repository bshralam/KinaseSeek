A cheminformatics pipeline for ligand-based virtual screening using publicly available chemical data (ChEMBL and ZINC). 

Objectives:

1) Parse and clean real molecular datasets from ChEMBL(kinase inhibitors) and ZINC(virtual screening library).
2) Load molecular structures from SMILES using RDKit.
3) Compute Morgan fingerprints.
4) Calculate average  Tanimoto similarity for each ZINC molecule vs. kinase references.
5) Output a ranked list of potential kinase-like hits.

