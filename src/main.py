from rdkit import Chem


def substructure_search(mols, mol): 
    result = []
    try:
        for test_mol in mols:
            if Chem.MolFromSmiles(test_mol).HasSubstructMatch(Chem.MolFromSmiles(mol)):
                result.append(test_mol)
        return result
    except:
        print("Something went wrong")

        
print(substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1"))