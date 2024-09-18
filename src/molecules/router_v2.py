from fastapi import APIRouter, HTTPException, status, UploadFile, Depends
from rdkit import Chem
from .dao import MoleculeDAO
from .schemas import MoleculeAdd, MoleculeResponse, valid_smiles
from typing import List, Annotated
import json


def substructure_search(mol_list: List[MoleculeResponse], smiles: str):
    result = []

    if not valid_smiles(smiles):
        raise ValueError

    query_mol = Chem.MolFromSmiles(smiles)

    for record in mol_list:
        mol = Chem.MolFromSmiles(record.smiles)
        if mol.HasSubstructMatch(query_mol):
            result.append(record)
    return result


NOT_FOUND_EXEPTION = HTTPException(
    status.HTTP_404_NOT_FOUND,
    detail="Molecule was not found",
)


dao_dependencie = Annotated[MoleculeDAO, Depends(MoleculeDAO)]


router = APIRouter()


@router.get("/")
async def get_molecules(dao: dao_dependencie) -> List[MoleculeResponse]:
    return await dao.find_all_molecules()


@router.get("/search")
async def get_molecules_by_substructure(
    smiles: str, dao: dao_dependencie
) -> List[MoleculeResponse]:
    try:
        mol_list = await dao.find_all_molecules()
        return substructure_search(mol_list, smiles)
    except ValueError:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail="Invalid SMILES structure"
        )


@router.get("/{id}")
async def get_molecule_by_id(
    id: int, dao: dao_dependencie
) -> MoleculeResponse | None:
    molecule = await dao.find_full_data(id)
    if molecule is None:
        raise NOT_FOUND_EXEPTION
    return molecule


@router.post("/", status_code=status.HTTP_201_CREATED)
async def add_molecule(
    molecule_to_add: MoleculeAdd,
    dao: dao_dependencie
) -> dict:
    molecule_id = await dao.add_molecule(
        **molecule_to_add.model_dump()
    )
    return {"message": f"Molecule {molecule_id} was added to database"}


@router.post("/file")
async def upload_file(file: UploadFile, dao: dao_dependencie):
    try:
        content = await file.read()
        await file.close()
        parsed_data = json.loads(content)
        valid_data = [MoleculeAdd(**mol) for mol in parsed_data]
        instances = [mol.model_dump() for mol in valid_data]
        return await dao.add_many_molecules(instances)
    except json.JSONDecodeError:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail="Invalid JSON"
        )
    except Exception as err:
        print(err)
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail="Could not parse the file data"
        )


@router.put("/{id}")
async def update_molecule(
    id: int,
    molecule_data: MoleculeAdd,
    dao: dao_dependencie
) -> dict:
    molecule_id = await dao.update_molecule_by_id(
        id, **molecule_data.model_dump()
    )
    if molecule_id is None:
        raise NOT_FOUND_EXEPTION
    return {"message": f"Molecule {molecule_id} was updated"}


@router.delete("/{id}")
async def delete_molecule(id: int, dao: dao_dependencie) -> dict:
    molecule_id = await dao.delete_molecule_by_id(id)
    if molecule_id is None:
        raise NOT_FOUND_EXEPTION
    return {"message": f"Molecule {molecule_id} was deleted"}
