from rdkit import Chem
from fastapi import FastAPI, HTTPException, UploadFile
from models import MolRecord
import json


molecule_db = [
    {
        "id": 0, 
        "smiles": "N[C@@H](CCCNC(=N)N)C(=O)O",
        "name": "L-Arginine"
    },

    {
        "id": 1,
        "smiles": "CN(CC(=O)O)C(=N)N",
        "name": "Creatine"
    },
    {
        "id": 2,
        "name": "Epinephrine",
        "smiles": "CC(C1=CC(=C(C=C1)O)O)NCC(C2=CC=CC=C2)O"
    }
]


def substructure_search(mol_records, mol):
    result = []
    query_mol = Chem.MolFromSmiles(mol)

    for record in mol_records:
        mol = Chem.MolFromSmiles(record.get("smiles"))
        if mol.HasSubstructMatch(query_mol):
            result.append(record)
    return result


def auto_increment(db):
    id = 0
    if db:
        id = db[-1].get("id") + 1
    return id


def valid_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return bool(mol)



app = FastAPI()



@app.post("/api/v1/molecules", tags=["Molecule Routes"], status_code=201)
def add_molecule(molecule: MolRecord):
    '''
    Adds a new record of a molecule to the database:
    - Requires a **molecule object**.
    - The **molecule object** must contain the **name** attribute, which is the trivial name of the molecule.
    - The **molecule object** must contain the **smiles** attribute, which is the SMILES expression of the molecule.
    - The **id** for the molecule is generated automatically.
    - The request returns the newly created **molecule dictionary** as the response.
    '''
    
    if not valid_smiles(molecule.smiles): 
        raise HTTPException(status_code=400, detail="Invalid SMILES format")
    
    mol_record = molecule.model_dump()

    # Automatic ID generation
    # If database is empty starting ID will be 0, else ID will increase by 1
    mol_record["id"] = auto_increment(molecule_db)
    molecule_db.append(mol_record)
    return mol_record





# Try with Guanidino group - C(=N)N
@app.get('/api/v1/molecules/search', tags=["Molecule Routes"])
def get_molecules_by_substructure(substructure: str): 
    '''
    Retrieves all records of molecules that contain the provided substructure:
    - Requires passing a query **substructure**.
    - The substructure must be in a valid **SMILES** format.
    - Every molecule in the database that contains the **substructure** will be returned in the response.
    - If no such molecules exist in the database, an **empty list** will be returned.
    - If the provided substructure does not have a valid **SMILES** format, an HTTPException will be raised.
    '''
    if not valid_smiles(substructure): 
        raise HTTPException(status_code=400, detail="Invalid SMILES format")
    
    return substructure_search(molecule_db, substructure)

        
    





@app.get("/api/v1/molecules/{id}", tags=["Molecule Routes"])
def get_molecule_by_id(id: int):
    '''
    Retrieves a single record of a molecule based on the ID:
    - This route requires passing of molecule **id**.
    - The **id** is a unique identifier for the molecule.
    - The request will attempt to find the molecule based on the provided **id**.
    - If the molecule exists, it will be returned in the response.
    - If the molecule does not exist, an HTTPException will be raised.
    '''
    
    molecule_found = False

    for mol in molecule_db:
        if mol["id"] == id:
            molecule_found = True
            return mol
    else:
        if not molecule_found:
            raise HTTPException(status_code=404, detail="Molecule does not exist")
        




@app.put("/api/v1/molecules/{id}", tags=["Molecule Routes"])
def update_molecule(id: int, molecule: MolRecord):
    '''
    Updates record of a molecule:
    - Requires valid **molecule object** and molecule **id**.
    - The request will attempt to find the molecule based on the provided **id**.
    - If the molecule exists, it will be **replaced** by the provided **molecule dictionary**
    - Response will return the updated record.
    - If the molecule does not exist, an HTTPException will be raised.
    '''

    if not valid_smiles(molecule.smiles): 
        raise HTTPException(status_code=400, detail="Invalid SMILES format")

    molecule_found = False

    for index, mol in enumerate(molecule_db):
        if mol["id"] == id:
            molecule_found = True
            mol_record = molecule.model_dump()
            mol_record["id"] = id
            molecule_db[index] = mol_record
            return mol_record
    else:
        if not molecule_found:
            raise HTTPException(status_code=404, detail="Molecule does not exist")





@app.delete("/api/v1/molecules/{id}", tags=["Molecule Routes"])
def delete_molecule(id: int):
    '''
    Deletes a single record of a molecule based on the ID:
    - This route requires passing the molecule **id**.
    - The request will attempt to find the molecule based on the provided **id**.
    - If the molecule exists, it will be **deleted** and the response will confirm the deletion.
    - If the molecule does not exist, an HTTPException will be raised.
    '''

    molecule_found = False

    for index, mol in enumerate(molecule_db):
        if mol["id"] == id:
            molecule_found = True
            molecule_db.pop(index)
            return mol
    else:
        if not molecule_found:
            raise HTTPException(status_code=404, detail="Molecule does not exist")
        





@app.get("/api/v1/molecules", tags=["Molecule Routes"])
def get_molecules(skip: int=0, limit: int | None = None):
    '''
    Retrieves whole list of Molecules
    - Takes two optional query parameters: **skip** and **limit**
    - Returns all molecules if no parameter was provided
    - Returns **slice** of the molecule database if any of the parameters is provided
    '''

    return molecule_db[skip : skip + (limit or len(molecule_db))]




# For testing this route there is a db.json file in src directory
@app.post("/api/v1/molecules/file", tags=["Molecule Routes"])
async def upload_file(file: UploadFile):

    '''
    Takes a JSON file as an input and adds its content to the **molecule database**
    - File upload only work with JSON string format\n
    Example of input:\n
        [{
            "id": 0, 
            "smiles": "N[C@@H](CCCNC(=N)N)C(=O)O",
            "name": "L-Arginine"
        },

        {
            "id": 1,
            "smiles": "CN(CC(=O)O)C(=N)N",
            "name": "Creatine"
        }]
    '''

    content = await file.read()
    await file.close()
    parsed_data = json.loads(content)

    for mol_record in parsed_data:
        if not valid_smiles(mol_record.get("smiles")): 
            raise HTTPException(status_code=400, detail="Invalid SMILES format")
        
        id = auto_increment(molecule_db) # Automatic ID generation
        mol_record["id"] = id
        molecule_db.append(mol_record)

    return molecule_db