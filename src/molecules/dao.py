from sqlalchemy import delete, update
from sqlalchemy.future import select
from .models import Molecule
from ..database import async_session_maker


class MoleculeDAO:
    model = Molecule
    session_maker = async_session_maker

    @classmethod
    async def find_all_molecules(cls):
        async with cls.session_maker() as session:
            query = select(cls.model)
            molecules = await session.execute(query)
            return molecules.scalars().all()

    @classmethod
    async def find_full_data(cls, molecule_id):
        async with cls.session_maker() as session:
            # Query to get molecule info
            query = select(cls.model).filter_by(id=molecule_id)
            result = await session.execute(query)
            molecule_info = result.scalar_one_or_none()

            # If molecule is not found, return None
            if not molecule_info:
                return None

            return molecule_info

    @classmethod
    async def add_molecule(cls, **molecule_data: dict):
        async with cls.session_maker() as session:
            async with session.begin():
                new_molecule = cls.model(**molecule_data)
                session.add(new_molecule)
                await session.flush()
                new_molecule_id = new_molecule.id
                await session.commit()
                return new_molecule_id

    @classmethod
    async def add_many_molecules(cls, instances: list[dict]):
        async with cls.session_maker() as session:
            async with session.begin():
                new_instances = [cls.model(**values) for values in instances]
                session.add_all(new_instances)
                await session.flush()
                await session.commit()
                return new_instances

    @classmethod
    async def delete_molecule_by_id(cls, molecule_id: int):
        async with cls.session_maker() as session:
            async with session.begin():
                query = select(cls.model).filter_by(id=molecule_id)
                result = await session.execute(query)
                molecule_to_delete = result.scalar_one_or_none()

                if not molecule_to_delete:
                    return None

                # Delete the molecule
                await session.execute(delete(cls.model)
                                      .filter_by(id=molecule_id))

                await session.commit()
                return molecule_id

    @classmethod
    async def update_molecule_by_id(
        cls, molecule_id: int,
        **molecule_data: dict
    ):
        async with cls.session_maker() as session:
            async with session.begin():
                query = select(cls.model).filter_by(id=molecule_id)
                result = await session.execute(query)
                molecule_to_update = result.scalar_one_or_none()

                if not molecule_to_update:
                    return None

                update_query = (
                    update(cls.model)
                    .where(cls.model.id == molecule_id)
                    .values(**molecule_data)
                    .execution_options(synchronize_session="fetch")
                )

                await session.execute(update_query)

                await session.commit()
                return molecule_id
