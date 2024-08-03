FROM continuumio/miniconda3

WORKDIR /src

COPY requirments.txt .

RUN pip install -r requirments.txt

COPY . .

EXPOSE 8000

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--reload"]