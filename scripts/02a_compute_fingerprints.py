from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import pandas as pd

# Load known kinase inhibitors
df = pd.read_csv("../data/test.csv")
df = df.dropna(subset=['canonical_smiles'])  # Remove rows with missing SMILES
df['mol'] = df['canonical_smiles'].apply(lambda s: Chem.MolFromSmiles(str(s)))
ref_fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) for m in df['mol'] if m]

# Load ZINC SMILES
with open("../data/zinc_subset_production.smi", 'r') as f:
    smiles_list = [line.strip() for line in f]

zinc_mols = [Chem.MolFromSmiles(s) for s in smiles_list]
zinc_fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) for m in zinc_mols if m]

# Compute average Tanimoto similarity to reference inhibitors
def avg_similarity(fp, refs):
    sims = [DataStructs.TanimotoSimilarity(fp, ref) for ref in refs]
    return sum(sims) / len(sims)

ranked = [(smiles_list[i], avg_similarity(fp, ref_fps)) for i, fp in enumerate(zinc_fps)]
ranked.sort(key=lambda x: x[1], reverse=True)

# Save results
with open("../data/top_ranked_hits_test.csv", "w") as f:
    f.write("smiles,avg_tanimoto\n")
    for smi, score in ranked[:20]:
        f.write(f"{smi},{score:.4f}\n")

