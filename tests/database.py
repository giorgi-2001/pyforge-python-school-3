from sqlalchemy.ext.asyncio import create_async_engine, async_sessionmaker
from src.molecules.models import Base
from src.molecules.dao import MoleculeDAO
from src.molecules.redis_cache import redis_client


DB_URL = "sqlite+aiosqlite:///:memory:"


engine = create_async_engine(DB_URL)
mock_async_session_maker = async_sessionmaker(engine, expire_on_commit=False)


class MockDAO(MoleculeDAO):
    session_maker = mock_async_session_maker


async def setup_function():
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)


async def teardown_function():
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.drop_all)
    redis_client.flushall()
