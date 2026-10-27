# -*- coding: utf-8 -*-
"""
任务三：生成 Morgan 分子指纹，并导出 descriptors.csv 和 fingerprints.csv
参数: radius=2, fpSize=2048, includeChirality=False
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import csv
import os

mols_data = [
    ("乙醇", "CCO"),
    ("苯", "c1ccccc1"),
    ("苯酚", "Oc1ccccc1"),
]

base_dir = r"C:\Users\杨涵\Doubao\chats\2026-10-27\new-chat-2"

# ========== 1. 计算描述符 ==========
descriptors_rows = []
for name, smi in mols_data:
    mol = Chem.MolFromSmiles(smi)
    molwt = round(Descriptors.MolWt(mol), 4)
    logp = round(Descriptors.MolLogP(mol), 4)
    hbd = Descriptors.NumHDonors(mol)
    descriptors_rows.append([name, molwt, logp, hbd])

# ========== 2. 生成 Morgan 指纹 ==========
fps = {}
for name, smi in mols_data:
    mol = Chem.MolFromSmiles(smi)
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=2048, useChirality=False
    )
    fps[name] = fp
    print(f"{name}: 指纹长度 = {fp.GetNumBits()} 位, 置位数量 = {fp.GetNumOnBits()}")

# 乙醇指纹前 20 位
ethanol_fp = fps["乙醇"]
first_20 = [ethanol_fp[i] for i in range(20)]
print(f"\n乙醇指纹前 20 位: {''.join(map(str, first_20))}")
print(f"逐位: {first_20}")

# ========== 3. 导出 descriptors.csv ==========
desc_path = os.path.join(base_dir, "descriptors.csv")
with open(desc_path, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.writer(f)
    writer.writerow(["分子", "MolWt", "MolLogP", "氢键供体数"])
    writer.writerows(descriptors_rows)
print(f"\ndescriptors.csv 已保存: {desc_path}")

# ========== 4. 导出 fingerprints.csv ==========
fp_path = os.path.join(base_dir, "fingerprints.csv")
with open(fp_path, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.writer(f)
    header = ["分子"] + [f"bit_{i}" for i in range(2048)]
    writer.writerow(header)
    for name, _ in mols_data:
        row = [name] + [int(fps[name][i]) for i in range(2048)]
        writer.writerow(row)
print(f"fingerprints.csv 已保存: {fp_path}")
