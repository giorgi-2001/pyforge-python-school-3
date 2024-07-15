from rdkit import Chem


def substructure_search(mols, mol):
    return [test_mol for test_mol in mols if Chem.MolFromSmiles(test_mol).HasSubstructMatch(Chem.MolFromSmiles(mol))]

        
print(substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1"))