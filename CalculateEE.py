import pandas as pd
import psi4
import pubchempy as pcp
import numpy as np
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
# Set the memory limit to 8 GB
#psi4.set_memory(int(8e9))
'''

# Read the data from the excel file
df = pd.read_excel("MoleculeList.xlsx", sheet_name="Sheet1")

# Get the column containing the strings
column = df["GetSDF"]

# Convert the column into a list
string_list = column.tolist()

# Remove the newline characters from the strings and convert them into a proper list
proper_list = [string.strip() for string in string_list]
cid_numbers = [int(line.split("(")[1].split(")")[0]) for line in proper_list]
'''
# Loop through the list of CID numbers
results = []
excitation_energies = []

cid_numbers = [9161, 1983, 175, 176, 180, 7410, 8040, 2244, 5897, 187, 6060, 6326, 6581, 6579, 6581, 7855, 135398513, 190, 60961, 196, 5816, 3385, 3672, 5950, 10037, 24409, 2124, 142249, 135413522, 2761507, 6115, 8418, 5354495, 16002, 6322, 2537, 5944, 83127, 8606, 767, 10219816, 5280489, 439570, 5281515, 29011, 289, 135403650, 45356234, 24646, 1732, 6372, 6327, 24638, 5280795, 5997, 24461, 23976, 444539, 73805498, 5754, 222786, 14769, 323, 586, 342, 442634, 6629, 969516, 9999, 9250, 8078, 121304016, 5780, 92934, 3035, 3036, 16318, 9997, 5462310, 189691, 11, 668, 12661, 1068, 22201, 1493, 31222, 6950, 5816, 83127, 5757, 5870, 5870, 6324, 3283, 174, 6354, 24458, 3386, 284, 69507, 7362, 370, 5780, 33032, 33032, 750, 135398634, 14496915, 8370, 8058, 91820602, 9321, 260, 768, 14917, 24841, 784, 533, 402, 123332, 961, 126, 7420, 135, 11671, 61739, 14828, 14945, 3776, 329438, 10407437, 7153, 6134, 84571, 5363390, 3893, 6106, 22311, 727, 5280450, 5761, 92934, 4980, 867, 9937450, 1615, 7153, 7153, 4055, 667490, 297, 7294, 878, 6329, 1486, 7416, 1486, 62329, 156391, 931, 10154219, 943, 943, 946, 948, 379, 5870, 5870, 445639, 71081, 971, 6354, 19862327, 85989101, 1486, 12406, 6518, 8003, 14792, 995, 6869, 6140, 6371, 24404, 9162, 6954, 6953, 23670670, 5994, 145742, 4928, 6334, 6581, 6581, 1032, 6335, 62857, 1049, 24261, 11726, 6998, 21637538, 5951, 66241, 4445035, 159367, 15407, 1099, 5362487, 2733526, 875, 6013, 27588, 15625, 39929, 8028, 8030, 2723790, 6989, 14985, 137843, 6351, 8081, 8184, 1174, 31295, 5641, 3016, 62923, 54670067, 171548, 962, 7929]
#cid_numbers = [241, 297]
for cid in cid_numbers:
#for cid in cid_numbers[35:71]:

    compound = pcp.Compound.from_cid(cid)
    smiles = compound.isomeric_smiles
    
    # Generate a 2D molecule using RDKit
    ob_mol = Chem.MolFromSmiles(smiles)
    ob_mol = Chem.AddHs(ob_mol)
    AllChem.EmbedMolecule(ob_mol)
    
    # Get the molecular geometry in xyz format
    xyz = Chem.MolToXYZBlock(ob_mol)
    
    # Calculate the excitation energy using psi4
    #change this and alter excel
    #put code on github, export into chemcompute 
    #Nanohub = backup
    #def2 basis set
    #B3LYP reference
    """
    psi4.set_options({
        "basis": "sto-3g",
        "reference": "uhf",
        "scf_type": "pk",
        "mp2_type": "conv",
        "e_convergence": 1e-8
    })


    print(energy)
    excitation_energies.append((cid, energy))
    n_orbitals = wfn.nmo()
    print(n_orbitals)
    eigenvalues = wfn.epsilon_a()
    print(eigenvalues)
    """
    psi4.set_options({
        "basis": "def2-TZVP",
        "reference": "B3LYP*"
    })
    xyz_string = Chem.MolToXYZBlock(ob_mol)
    mol = psi4.core.Molecule.from_string(xyz_string)
    #energy, wfn = psi4.energy("mp2/sto-3g", molecule=mol, return_wfn=True)
    energy, scf_wfn = psi4.energy("B3LYP/def2-TZVP", molecule = mol, return_wfn=True)

    #HOMO = scf_wfn.epsilon_a_subset("AO", "ALL")[scf_wfn.nalpha()]
    HOMO = scf_wfn.epsilon_a_subset("AO", "ALL").np[scf_wfn.nalpha()-1]

    LUMO = scf_wfn.epsilon_a_subset("AO", "ALL").np[scf_wfn.nalpha()]
    print(HOMO)
    print(LUMO)
    orbitalEnergy = LUMO - HOMO
    print(orbitalEnergy)
    results.append((cid, orbitalEnergy))
    print(results)
   # Exit the loop after all CID numbers have been processed
    if cid == cid_numbers[-1]:
        break
    psi4.set_output_file('output.dat', False)

# Print the results
#for result in excitation_energies:
 #   print(f"CID {result[0]}: Excitation energies = {result[1]}")