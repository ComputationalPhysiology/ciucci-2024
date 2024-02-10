# Material from Ciucci 2024

This repository contains the source code to reproduce the modeling results from the paper
> Paper here



## Reproducible results
You can find a report of the results by going to the repository homepage at <https://computationalphysiology.github.io/ciucci-2024/myst/report-1>.

Note that the results displayed here are re-generated every time a push is made to the repository (see https://github.com/ComputationalPhysiology/ciucci-2024/actions/workflows/deploy_docs.yml). Note that all figures are also uploaded as artifacts during the build.


## Install dependencies
To run the code in this repo you need FEniCS and gmsh with OpenCascade installed. The easiest way to do this is to either use the provided docker image, i.e
```
docker pull ghcr.io/computationalphysiology/ciucci-2024:v1.0.0
```
You can create a new container using the command
```
docker run --name ciucci -w /home/shared -v $PWD:/home/shared -it ghcr.io/computationalphysiology/ciucci-2024:v1.0.0
```
which will also share your current directory with the container. If you are interested you can also check out the [Dockerfile](Dockerfile) if you want to know how the image was created.

## Running scripts

There are three types of simulation scripts.

1. A contracting cylinder, with only one time point
2. A contraction cylinder with a full twitch
3. An unloaded (zero pressure) and loaded (pressure = 15 kPA) contracting left ventricle.

We refer to the [documentation](https://computationalphysiology.github.io/ciucci-2024/myst/report-1) for more information about the models and numerical experiments.

The main script for running all the commands are located inside the `code` directory and you should first make sure to navigate to this directory, i.e
```
cd code
```

You can list all the available options by typing
```
python3 main.py --help
```
Here is an example of how it looks inside the container
```
root@b2637926d9de:/home/shared/code# python3 main.py --help
usage: main.py [-h] [--dry-run] {preprocess-lv,preprocess-cylinder,run-lv,run-cylinder,run-cylinder-twitch,postprocess-lv,postprocess-cylinder,postprocess-cylinder-twitch} ...

positional arguments:
  {preprocess-lv,preprocess-cylinder,run-lv,run-cylinder,run-cylinder-twitch,postprocess-lv,postprocess-cylinder,postprocess-cylinder-twitch}
    preprocess-lv       Create left ventricle mesh
    preprocess-cylinder
                        Create cylinder mesh
    run-lv              Run simulations with left ventricle model
    run-cylinder        Run simulations with cylinder model
    run-cylinder-twitch
                        Run simulations with cylinder model
    postprocess-lv      Postprocess LV results
    postprocess-cylinder
                        Postprocess cylinder results
    postprocess-cylinder-twitch
                        Postprocess cylinder results

options:
  -h, --help            show this help message and exit
  --dry-run             Just print the command and do not run it
```


Below we show the basic command. Note that you can also use the flag `--help` to see all the different options for each command, e.g
```
root@b2637926d9de:/home/shared/code# python3 main.py preprocess-cylinder --help
usage: main.py preprocess-cylinder [-h] [-o MESH_FOLDER] [-c CHAR_LENGTH] [-l LENGTH] [-r RADIUS]

options:
  -h, --help            show this help message and exit
  -o MESH_FOLDER, --mesh-folder MESH_FOLDER
  -c CHAR_LENGTH, --char_length CHAR_LENGTH
  -l LENGTH, --length LENGTH
  -r RADIUS, --radius RADIUS
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


## LV


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



## Citation
If you use the code in this repository, please cite:
TBW


## Author
Henrik Finsberg and Sam Wall


## License
MIT
