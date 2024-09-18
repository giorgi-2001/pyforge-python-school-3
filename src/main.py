from fastapi import FastAPI
from os import getenv
from .molecules.router_v1 import router as molecule_router_v1
from .molecules.router_v2 import router as molecule_router_v2


app = FastAPI()


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


app.include_router(
  molecule_router_v1,
  prefix="/api/v1/molecules",
  tags=["Molecules v1"]
)


app.include_router(
  molecule_router_v2,
  prefix="/api/v2/molecules",
  tags=["Molecules v2"]
)
