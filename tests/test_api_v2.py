from sqlalchemy.ext.asyncio import create_async_engine, async_sessionmaker
from fastapi.testclient import TestClient
from src.molecules.models import Base
from src.molecules.router_v2 import MoleculeDAO
from src.main import app
from io import BytesIO
import pytest
import json


DB_URL = "sqlite+aiosqlite:///:memory:"


engine = create_async_engine(DB_URL)
mock_async_session_maker = async_sessionmaker(engine, expire_on_commit=False)


async def setup_function():
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)


async def teardown_function():
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.drop_all)


class MockDAO(MoleculeDAO):
    session_maker = mock_async_session_maker


app.dependency_overrides[MoleculeDAO] = MockDAO


BASE_URL = "/api/v2/molecules"


@pytest.fixture(scope="module")
def client():
    return TestClient(app)


# Testing API

@pytest.mark.asyncio
async def test_list_all_molecules(client: TestClient):
    await setup_function()
    response = client.get(BASE_URL)
    await teardown_function()
    assert response.status_code == 200
    data = response.json()
    assert data["molecules"] == []


@pytest.mark.asyncio
async def test_add_molecule(client: TestClient):
    await setup_function()
    mol = {"name": "methane", "smiles": "C"}
    response = client.post(BASE_URL, json=mol)
    await teardown_function()
    assert response.status_code == 201
    assert response.json() == {"message": "Molecule 1 was added to database"}


@pytest.mark.asyncio
async def test_get_molecules_by_substructure(client: TestClient):
    await setup_function()
    mol = {"name": "ethanol", "smiles": "CCO"}
    client.post(BASE_URL, json=mol)
    response = client.get(BASE_URL + "/search?smiles=O")
    await teardown_function()
    assert response.status_code == 200
    list_with_mol = response.json()
    mol_from_list = list_with_mol[0]
    assert mol_from_list.get("id") == 1
    assert mol_from_list.get("name") == "ethanol"
    assert mol_from_list.get("smiles") == "CCO"


@pytest.mark.asyncio
async def test_get_molecule_by_id(client: TestClient):
    await setup_function()
    mol = {"name": "methane", "smiles": "C"}
    client.post(BASE_URL, json=mol)
    response = client.get(BASE_URL + "/1")
    await teardown_function()
    assert response.status_code == 200
    methane = response.json()
    assert methane.get("id") == 1
    assert methane.get("name") == "methane"
    assert methane.get("smiles") == "C"


@pytest.mark.asyncio
async def test_update_molecule(client: TestClient):
    await setup_function()
    mol = {"name": "methane", "smiles": "C"}
    updated_mol = {"name": "ethanol", "smiles": "CCO"}
    client.post(BASE_URL, json=mol)
    response = client.put(BASE_URL + "/1", json=updated_mol)
    await teardown_function()
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule 1 was updated"}


@pytest.mark.asyncio
async def test_delete_molecule(client: TestClient):
    await setup_function()
    mol = {"name": "methane", "smiles": "C"}
    client.post(BASE_URL, json=mol)
    response = client.delete(BASE_URL + "/1")
    await teardown_function()
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule 1 was deleted"}


@pytest.mark.asyncio
async def test_upload_file(client: TestClient):
    await setup_function()
    with open("./tests/mols.json", "rb") as file:
        response = client.post(BASE_URL + "/file", files={
            "file": ("molecules", file, "multipart/formdata")
        })
    await teardown_function()
    assert response.status_code == 200


# Testing 404 errors

@pytest.mark.asyncio
async def test_get_molecule_by_id_with_404(client: TestClient):
    await setup_function()
    response = client.get(BASE_URL + "/1")
    await teardown_function()
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule was not found"}


@pytest.mark.asyncio
async def test_delete_molecule_with_404(client: TestClient):
    await setup_function()
    response = client.delete(BASE_URL + "/1")
    await teardown_function()
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule was not found"}


@pytest.mark.asyncio
async def test_update_molecule_with_404(client: TestClient):
    await setup_function()
    mol = {"name": "methane", "smiles": "C"}
    response = client.put(BASE_URL + "/1", json=mol)
    await teardown_function()
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule was not found"}


# Testing pydantic validation errors

@pytest.mark.asyncio
async def test_add_molecule_pydantic(client: TestClient):
    await setup_function()
    response = client.post(BASE_URL, json="")
    await teardown_function()
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_update_molecule_pydantic(client: TestClient):
    await setup_function()
    mol = {"name": "methane", "smiles": "C"}
    updated_mol = {"name": "ethanol", "smiles": "&&"}
    client.post(BASE_URL, json=mol)
    response = client.put(BASE_URL + "/1", json=updated_mol)
    await teardown_function()
    assert response.status_code == 422


# Testing file upload errors
@pytest.mark.asyncio
async def test_upload_file_with_invalid_smiles(client: TestClient):
    await setup_function()
    mol = {"name": "methane", "smiles": "C&&^"}
    mock_file_content = json.dumps([mol]).encode()
    mock_file = BytesIO(mock_file_content)
    response = client.post(BASE_URL + "/file", files={
        "file": ("molecules", mock_file, "application/json")
    })
    await teardown_function()
    assert response.status_code == 400
    assert response.json() == {"detail": "Could not parse the file data"}


@pytest.mark.asyncio
async def test_upload_file_with_invalid_file_format(client: TestClient):
    await setup_function()
    mock_file_content = b"molecule"
    mock_file = BytesIO(mock_file_content)
    response = client.post("/api/v1/molecules/file", files={
        "file": ("molecules", mock_file, "text/plain")
    })
    await teardown_function()
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid JSON"}


# Testing substructure search with Invalid smiles

@pytest.mark.asyncio
async def test_get_molecules_by_sub_invalid_smiles(client: TestClient):
    await setup_function()
    response = client.get(BASE_URL + "/search?smiles=&&")
    await teardown_function()
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES structure"}
