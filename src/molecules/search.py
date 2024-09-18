from rdkit import Chem
from src.molecules.schemas import valid_smiles
from src.molecules.dao import MoleculeDAO
from typing import AsyncGenerator
import asyncio
import functools
import threading


async def substructure_search(mol_list: AsyncGenerator, smiles: str):
    if not valid_smiles(smiles):
        raise ValueError("Invalid SMILES")

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
            result.append(mol.model_dump())

        current_index += 1

    return result


def _run_coroutine_in_thread(func, args, kwargs):
    result_container = {}
    exception_container = {}

    def run():
        try:
            new_loop = asyncio.new_event_loop()
            asyncio.set_event_loop(new_loop)
            result = new_loop.run_until_complete(func(*args, **kwargs))
            result_container["result"] = result
        except Exception as e:
            exception_container["exception"] = e
        finally:
            new_loop.close()

    new_thread = threading.Thread(target=run)
    new_thread.start()
    new_thread.join()

    exception = exception_container.get("exception")

    if exception:
        raise exception

    return result_container.get("result")


def async_to_sync(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            loop = asyncio.get_event_loop()
        except RuntimeError:
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)

        if loop.is_running():
            return _run_coroutine_in_thread(func, args, kwargs)

        return loop.run_until_complete(func(*args, **kwargs))

    return wrapper


@async_to_sync
async def search_algorithm(
    smiles: str,
    page: int,
    limit: int,
    Dao=MoleculeDAO
) -> dict:
    dao = Dao()
    all_db_mols = dao.find_all_molecules()

    try:
        contains_substructure = substructure_search(all_db_mols, smiles)
    except ValueError:
        raise

    result = await paginate_sub_search(
        contains_substructure,
        page=page,
        limit=limit
    )

    return result
