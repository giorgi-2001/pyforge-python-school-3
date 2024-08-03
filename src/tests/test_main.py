import pytest
from fastapi.testclient import TestClient
from app.main import app, molecule_db, auto_increment, valid_smiles, substructure_search
from app.models import MolRecord
import json
from io import BytesIO


empty_db = []
mock_db = molecule_db[:]
contains_guanidino_group = molecule_db[:2]


@pytest.mark.parametrize("input, expected", [(empty_db, 0), (mock_db, 3)])
def test_auto_increment(input, expected):
    assert auto_increment(input) == expected


@pytest.mark.parametrize(
        "input, expected", 
        [("C(=N)C", True), ("&@", False), ("", False)]
    )
def test_valid_smiles(input, expected):
    assert valid_smiles(input) == expected


# Substructure search tests

def test_substructure_search_guanidino():
    assert substructure_search(mock_db, "C(=N)N") == contains_guanidino_group


def test_substructure_search_empty_str():
    assert substructure_search(mock_db, "") == []


def test_substructure_search_invalid_smiles():
    with pytest.raises(Exception):
        assert substructure_search(mock_db, "&@")


# testing model

@pytest.fixture()
def aspirin():
    return MolRecord(
        name="acetylsalicylic acid",
        smiles="CC(=O)OC1=CC=CC=C1C(=O)O"
    )


def test_mol_record(aspirin):
    assert aspirin.name == "acetylsalicylic acid"
    assert aspirin.smiles == "CC(=O)OC1=CC=CC=C1C(=O)O"



# Testing API Routes

client = TestClient(app)

nicotine = {
    "original": {
        "name": "Nicotine",
        "smiles": "CN1CCC[C@H]1C2=CN=CC=C2"
    },

    "with_id": {
        "id": 0,
        "name": "Nicotine",
        "smiles": "CN1CCC[C@H]1C2=CN=CC=C2"
    },

    "updated": {
        "id": 0,
        "name": "Nicotine updated",
        "smiles": "CN1CCC[C@H]1C2=CN=CC=C2"
    }
}


def setup_function():
    molecule_db.clear()


def test_get_server():
    response = client.get("/")
    assert response.status_code == 200
    assert response.json() == {"server_id": "1"}


def test_get_molecules():
    response = client.get("/api/v1/molecules")
    assert response.status_code == 200
    assert response.json() == []


def test_get_molecules_by_substructure():
    client.post("/api/v1/molecules", json=nicotine["original"])
    response = client.get("/api/v1/molecules/search?substructure=C2=CN=CC=C2")
    assert response.status_code == 200
    assert response.json() == [nicotine["with_id"]]


def test_get_molecule_by_id():
    client.post("/api/v1/molecules", json=nicotine["original"])
    response = client.get("/api/v1/molecules/0")
    assert response.status_code == 200
    assert response.json() == nicotine["with_id"]


def test_add_molecule():
    response = client.post("/api/v1/molecules", json=nicotine["original"])
    assert response.status_code == 201
    assert response.json() == nicotine["with_id"]


def test_update_molecule():
    client.post("/api/v1/molecules", json=nicotine["original"])
    response = client.put("/api/v1/molecules/0", json=nicotine["updated"])
    assert response.status_code == 200
    assert response.json() == nicotine["updated"]


def test_delete_molecule():
    client.post("/api/v1/molecules", json=nicotine["original"])
    response = client.delete("/api/v1/molecules/0")
    assert response.status_code == 200
    assert response.json() == nicotine["with_id"]


def test_upload_file():
    with open("../db.json", "rb") as file:
        response = client.post("/api/v1/molecules/file", files={
            "file": ("molecules", file, "application/json")
        })
    assert response.status_code == 200
    assert response.json() == molecule_db


# Testing 404 Exeptions 

def test_get_molecule_by_id_with_404_exeption():
    response = client.get("/api/v1/molecules/1")
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule does not exist"}


def test_update_molecule_with_404_exeption():
    response = client.put("/api/v1/molecules/1", json=nicotine["updated"])
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule does not exist"}


def test_delete_molecule_with_404_exeption():
    response = client.delete("/api/v1/molecules/1")
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule does not exist"}


# Testing pydntic exeptions

def test_get_molecules_by_substructure_with_pydantic_exeption():
    response = client.get("/api/v1/molecules/search")
    assert response.status_code == 422

def test_add_molecule_with_pydantic_exeption():
    response = client.post("/api/v1/molecules", json={})
    assert response.status_code == 422

def test_update_molecule_with_pydantic_exeption():
    response = client.put("/api/v1/molecules/0", json={})
    assert response.status_code == 422


# Testing Invalid SMILES expetions

invalid_smiles_molecule = {
    "name": "invalid",
    "smiles": "&A"
}


def test_get_molecules_by_substructure_with_invalid_smiles():
    response = client.get("/api/v1/molecules/search?substructure=&A")
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES format"}


def test_add_molecule_with_invalid_smiles():
    response = client.post("/api/v1/molecules", json=invalid_smiles_molecule)
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES format"}


def test_update_molecule_with_invalid_smiles():
    response = client.put("/api/v1/molecules/0", json=invalid_smiles_molecule)
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES format"}


def test_upload_file_with_invalid_smiles():
    mock_file_content = json.dumps([invalid_smiles_molecule]).encode()
    mock_file = BytesIO(mock_file_content)
    response = client.post("/api/v1/molecules/file", files={
        "file": ("molecules", mock_file, "application/json")
    })
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES format"}


def test_upload_file_with_invalid_file_format():
    mock_file_content = b"molecule"
    mock_file = BytesIO(mock_file_content)
    response = client.post("/api/v1/molecules/file", files={
        "file": ("molecules", mock_file, "text/plain")
    })
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid JSON"}