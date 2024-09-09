from src.celery_worker import celery_app
from src.molecules.search import search_algorithm


@celery_app.task
def add_task(smiles: str, page: int, limit: int):
    return search_algorithm(smiles, page, limit)
