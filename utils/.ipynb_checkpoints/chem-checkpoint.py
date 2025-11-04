from typing import Iterable, Optional, List
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

def mol_from_smiles(smiles: str) -> Optional[Chem.Mol]:
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.SanitizeMol(mol)
    return mol

def morgan_fp(mol: Chem.Mol, radius: int = 2, nBits: int = 2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)

def tanimoto(a, b) -> float:
    return DataStructs.TanimotoSimilarity(a, b)

def pairwise_tanimoto(fps: List) -> np.ndarray:
    n = len(fps)
    mat = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        mat[i, i] = 1.0
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i+1:])
        mat[i, i+1:] = sims
        mat[i+1:, i] = sims
    return mat
