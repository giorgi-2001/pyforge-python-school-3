from pydantic import BaseModel, Field, field_validator
from datetime import datetime
from rdkit import Chem


def valid_smiles(smiles: str):
    if not smiles:  # cheking for empty string
        return False
    mol = Chem.MolFromSmiles(smiles)
    return bool(mol)


class MoleculeAdd(BaseModel):
    name: str = Field(
        ..., min_length=1, max_length=100, description="Drug name"
    )
    smiles: str = Field(
        ..., min_length=1, max_length=100,
        description="structure of chemical molecules"
    )

    @field_validator("smiles")
    def validate_smile(cls, smiles):
        if not valid_smiles(smiles):
            raise ValueError("Invalid SMILES structure")
        return smiles


class MoleculeResponse(BaseModel):
    name: str
    smiles: str
    id: int
    created_at: datetime
    updated_at: datetime