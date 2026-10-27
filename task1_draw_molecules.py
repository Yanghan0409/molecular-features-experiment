# -*- coding: utf-8 -*-
"""
任务一：从 SMILES 生成结构图
用 RDKit 读取乙醇、苯、苯酚的 SMILES，生成包含三个分子的二维结构图并标注名称。
"""

from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

# 三个分子：名称 + SMILES
mols_data = [
    ("乙醇 (Ethanol)", "CCO"),
    ("苯 (Benzene)", "c1ccccc1"),
    ("苯酚 (Phenol)", "Oc1ccccc1"),
]

# 解析 SMILES 为 RDKit 分子对象
mols = []
legends = []
for name, smi in mols_data:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"[ERROR] 无法解析 SMILES: {smi} ({name})")
        continue
    mols.append(mol)
    legends.append(name)
    print(f"[OK] {name}: {smi} -> 原子数={mol.GetNumAtoms()}, 键数={mol.GetNumBonds()}")

# 生成 1 行 3 列的结构图
img = Draw.MolsToGridImage(
    mols,
    molsPerRow=3,
    subImgSize=(400, 300),
    legends=legends,
    useSVG=False,
)

# 保存为 PNG 图片
output_path = r"C:\Users\杨涵\Doubao\chats\2026-10-27\new-chat-2\molecules_task1.png"
img.save(output_path)
print(f"\n结构图已保存至: {output_path}")
