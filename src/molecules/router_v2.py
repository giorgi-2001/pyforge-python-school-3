from fastapi import APIRouter, HTTPException, status, UploadFile, Depends, Query
from rdkit import Chem
from ..logger import logger
from .dao import MoleculeDAO
from .schemas import (
    MoleculeAdd, MoleculeResponse,
    MolResWithPagination, PaginationData,
    valid_smiles
)
from typing import List, Annotated, AsyncGenerator
import json
import math


async def substructure_search(mol_list: AsyncGenerator, smiles: str):
    if not valid_smiles(smiles):
        raise ValueError

    query_mol = Chem.MolFromSmiles(smiles)

    async for record in mol_list:
        mol = Chem.MolFromSmiles(record.smiles)
        if mol.HasSubstructMatch(query_mol):
            yield record
            


async def paginate_sub_search(
    mol_generator: AsyncGenerator,
    page: int = 1,
    limit: int = 2
):
    result = []
    start_index = (page - 1) * limit
    end_index = page * limit
    current_index = 0
    while current_index < end_index:
        try:
            mol = await anext(mol_generator)
        except StopAsyncIteration:
            break

        if current_index >= start_index:
            result.append(mol)

        current_index += 1
        
    return result


NOT_FOUND_EXEPTION = HTTPException(
    status.HTTP_404_NOT_FOUND,
    detail="Molecule was not found",
)


dao_dependencie = Annotated[MoleculeDAO, Depends(MoleculeDAO)]


router = APIRouter()


@router.get("/")
async def list_all_molecules(
    dao: dao_dependencie,
    page: int = Query(ge=1, default=1),
    limit: int = Query(ge=1, le=100, default=2)
) -> MolResWithPagination:
    count = await dao.molecule_count()
    skip_value = (page - 1) * limit
    all_pages = math.ceil(count / limit)
    has_next_page = page < all_pages 

    mol_list_generator = dao.find_all_molecules(
        offset=skip_value, limit=limit
    )

    return {
        "pagination_data": PaginationData(
            current_page=page, page_count=all_pages,
            molecule_count=count, page_size=limit,
            has_next_page=has_next_page
        ),

        "molecules": [molecule async for molecule in mol_list_generator]
    }



@router.get("/search")
async def get_molecules_by_substructure(
    smiles: str, dao: dao_dependencie,
    page: int = Query(ge=1, default=1),
    limit: int = Query(ge=1, le=100, default=2)
) -> List[MoleculeResponse]:
    try:
        mol_list_generator = dao.find_all_molecules()
        mol_result_generator = substructure_search(mol_list_generator, smiles)
        result = await paginate_sub_search(
            mol_generator=mol_result_generator,
            page=page, limit=limit
        )
        return result
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
