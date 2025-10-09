# imports
import os
import csv
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from ersilia_pack_utils.core import write_out, read_smiles

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# my model
def my_model(smiles_list):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, radius=3, nBits=2048) for mol in mols]
    array_fps = [np.array(fp) for fp in fps]
    return array_fps

cols, smiles_list = read_smiles(input_file)
# run model
outputs = my_model(smiles_list)

print(outputs)
print(type(outputs[0]))

#check input and output have the same lenght
input_len = len(smiles_list)
output_len = len(outputs)
assert input_len == output_len

headers= headers = ["dim_{0}".format(str(i).zfill(4)) for i in range(len(outputs[0]))]

write_out(outputs,headers,output_file,dtype='int32')
