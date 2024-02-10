FROM ghcr.io/scientificcomputing/fenics-gmsh:2023-11-15

WORKDIR /repo


COPY . /repo
RUN python3 -m pip install --upgrade pip && python3 -m pip install -r requirements.txt
