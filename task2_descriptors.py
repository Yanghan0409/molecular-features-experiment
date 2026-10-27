# -*- coding: utf-8 -*-
"""
任务二：计算三个分子的描述符
MolWt (分子量, g/mol), MolLogP (脂水分配系数), 氢键供体数 (HBD)
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

mols_data = [
    ("乙醇", "CCO"),
    ("苯", "c1ccccc1"),
    ("苯酚", "Oc1ccccc1"),
]

print(f"{'分子':<6} {'MolWt (g/mol)':<16} {'MolLogP':<12} {'氢键供体数':<10}")
print("-" * 50)

for name, smi in mols_data:
    mol = Chem.MolFromSmiles(smi)
    molwt = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    print(f"{name:<6} {molwt:<16.4f} {logp:<12.4f} {hbd:<10d}")
