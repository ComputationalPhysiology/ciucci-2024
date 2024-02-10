# Material from Ciucci 2024

This repository contains the source code to reproduce the modeling results from the paper
> Paper here



## Reproducible results
You can find a report of the results by going to the repository homepage at <https://computationalphysiology.github.io/ciucci-2024/myst/report-1>.

Note that the results displayed here are re-generated every time a push is make to the repository (see https://github.com/ComputationalPhysiology/ciucci-2024/actions/workflows/deploy_docs.yml). Note that all figures are also uploaded as artifacts during the build.


## Install dependencies
To run the code in this repo you need FEniCS and gmsh with OpenCascade installed. The easiest way to do this, is to use the [following docker image](https://github.com/scientificcomputing/packages/pkgs/container/fenics-gmsh).
Next, you need to install the requirements,
```
python3 -m pip install -r requirements.txt
```

## LV
TBW

### Preprocessing
```
python3 main.py preprocess-lv -o meshes/lv --psize-ref 3.0
```

### Running simulations
```
python3 main.py run-lv -i meshes/lv -o results/lv
```

### Postprocessing
```
python3 main.py postprocess-lv -i meshes/lv -r results/lv/ -o figures/lv
```
Get stats using
```
python3 main.py postprocess-lv -i meshes/lv -r results/lv/ -o figures/lv --print-stats
```
Create paraview files with stresses imposed on deformed geometry
```
python3 main.py postprocess-lv -i meshes/lv -r results/lv/ -o figures/lv --create-paraview
```

## Cylinder

### Preprocessing
```
python3 main.py preprocess-cylinder -o meshes/cylinder
```

### Simulation
```
python3 main.py run-cylinder -i meshes/cylinder -o results/cylinder
```

### Simulation (twitch)
```
python3 main.py run-cylinder-twitch -i meshes/cylinder -o results/cylinder-twitch
```

### Postprocessing
```
python3 main.py postprocess-cylinder -i meshes/cylinder -r results/cylinder -o figures/cylinder
```

### Postprocessing (twitch)
```
python3 main.py postprocess-cylinder-twitch -i meshes/cylinder -r results/cylinder-twitch -o figures/cylinder-twitch
```



## Citation
If you use the code in this repository, please cite:
TBW


## Author
Henrik Finsberg and Sam Wall


## License
MIT
