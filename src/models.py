from pydantic import BaseModel


class MolRecord(BaseModel):
    name: str
    smiles: str
