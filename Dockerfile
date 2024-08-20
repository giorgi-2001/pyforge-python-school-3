FROM continuumio/miniconda3

WORKDIR /app

COPY requirements.txt .

RUN pip install -r requirements.txt

COPY ./src ./src

COPY ./tests ./tests

EXPOSE 8000

CMD ["uvicorn", "src.main:app", "--host", "0.0.0.0", "--port", "8000", "--reload"]