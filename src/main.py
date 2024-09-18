from fastapi import FastAPI, Request, Response
from fastapi.responses import StreamingResponse
from .molecules.router_v1 import router as molecule_router_v1
from .molecules.router_v2 import router as molecule_router_v2
from .logger import logger, LogFactory
from os import getenv
import time
import json


app = FastAPI()


@app.middleware("http")
async def log_events(request: Request, call_next):
    start_time = time.time()
    response: StreamingResponse = await call_next(request)

    response_body = b""

    async for chunk in response.body_iterator:
        response_body += chunk

    process_time = time.time() - start_time

    log = LogFactory(
      status_code=response.status_code,
      ip_adress=request.client.host,
      method=request.method,
      path=request.url.path,
      process_time=process_time
    )

    if response.status_code >= 400:
        json_str = response_body.decode()
        msg_obj: dict = json.loads(json_str)
        log.set_message(msg_obj.get("detail"))
        logger.error(log.message_log)
    elif response.status_code < 400 and request.method != "GET":
        if response_body:
            json_str = response_body.decode()
            msg_obj: dict = json.loads(json_str)
            log.set_message(msg_obj.get("message"))
        logger.info(log.message_log)
    else:
        logger.info(log.base_log)

    return Response(
      content=response_body, status_code=response.status_code,
      headers=dict(response.headers), media_type=response.media_type
    )


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
