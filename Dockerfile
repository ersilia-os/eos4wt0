FROM bentoml/model-server:0.11.0-py310
MAINTAINER ersilia

RUN pip install rdkit==2023.9.2

WORKDIR /repo
COPY . /repo
