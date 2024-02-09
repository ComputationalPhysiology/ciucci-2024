# Material from Ciucci 2024

This repository contains the source code to reproduce the modeling results from the paper
> Paper here



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
TBW


## Cylinder

### Preprocessing
```
```

### Simulation
```
python3 main.py run-cylinder -i meshes/cylinder -o results/cylinder
```


## Citation
If you use the code in this repository, please cite:
TBW


## Author
Henrik Finsberg and Sam Wall


## License
MIT
