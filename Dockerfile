FROM bentoml/model-server:0.11.0-py310
MAINTAINER ersilia

RUN pip install rdkit==2023.9.2
RUN pip install ersilia-pack-utils==0.1.5

WORKDIR /repo
COPY . /repo
