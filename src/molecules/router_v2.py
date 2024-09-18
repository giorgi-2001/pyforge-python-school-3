from fastapi import (
    APIRouter, HTTPException, status,
    UploadFile, Depends, Query
)

from celery.result import AsyncResult
from celery import Task
from .dao import MoleculeDAO

from .schemas import (
    MoleculeAdd,
    MolResWithPagination,
    PaginationData,
    valid_smiles
)

from .redis_cache import get_cached_result, set_cache
from src.molecules.tasks import add_task
from ..celery_worker import celery_app
from typing import Annotated
import json
import math


NOT_FOUND_EXEPTION = HTTPException(
    status.HTTP_404_NOT_FOUND,
    detail="Molecule was not found",
)


dao_dependencie = Annotated[MoleculeDAO, Depends(MoleculeDAO)]
task_dependencie = Annotated[Task, Depends(lambda: add_task)]


router = APIRouter()


@router.get("/")
async def list_all_molecules(
    dao: dao_dependencie,
    page: int = Query(ge=1, default=1),
    limit: int = Query(ge=1, le=100, default=20)
) -> MolResWithPagination:
    cache_key = f"page={str(page)}&limit={str(limit)}"

    result: None | dict = get_cached_result(cache_key)

    if result:
        result.update({"source": "cache"})
        return result

    count = await dao.molecule_count()
    skip_value = (page - 1) * limit
    all_pages = math.ceil(count / limit)
    has_next_page = page < all_pages

    mol_list_generator = dao.find_all_molecules(
        offset=skip_value, limit=limit
    )

    pagination_data = PaginationData(
            current_page=page, page_count=all_pages,
            molecule_count=count, page_size=limit,
            has_next_page=has_next_page
        )

    response = {
        "source": "database",
        "pagination_data": pagination_data.model_dump(),
        "molecules": [
            molecule.model_dump() async for molecule in mol_list_generator
        ]
    }

    set_cache(cache_key, response, 180)
    return response


@router.post("/search")
async def add_search_task(
    add_task: task_dependencie,
    smiles: str,
    page: int = Query(ge=1, default=1),
    limit: int = Query(ge=1, le=100, default=20)
) -> dict:
    if not valid_smiles(smiles):
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail="Invalid SMILES structure"
        )

    task = add_task.delay(smiles, page, limit)
    return {"task_id": task.id, "status": task.status}


@router.get("/search/{task_id}")
async def get_search_results(task_id: str) -> dict:
    task_result = AsyncResult(task_id, app=celery_app)
    if task_result.state == 'PENDING':
        return {"task_id": task_id, "status": "Task is still processing"}
    elif task_result.state == 'SUCCESS':
        return {
            "task_id": task_id,
            "status": "Task completed",
            "result": task_result.result
        }
    else:
        return {"task_id": task_id, "status": task_result.state}


@router.get("/{id}")
async def get_molecule_by_id(
    id: int, dao: dao_dependencie
) -> dict:
    cache_key = f"id={id}"
    molecule = get_cached_result(cache_key)

    if molecule:
        return {"source": "cache", "result": molecule}

    molecule = await dao.find_full_data(id)

    if molecule is None:
        raise NOT_FOUND_EXEPTION

    set_cache(cache_key, molecule.model_dump())
    return {"source": "database", "result": molecule.model_dump()}


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
        await dao.add_many_molecules(instances)
        return {"message": "File was processed. Molecules were uploaded."}
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
