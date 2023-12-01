# Folder structure

The code is structured as follow
```
AURORA/
├── data/
│   └── 20220926/...
│   └── ...
│
├── docs/...
│ 
├── e_cascading_data/
│   └── N2/...
│   └── O2/...
│   └── O/...
│
├── scripts/
│   └── Control_script.jl
│   └── ...
│ 
└── src/
    └── main.jl
    └── cascading.jl
    └── setup.jl
    └── ...
```
The folder `data/` contains the subfolders where simulation results are saved.

The folder `docs/` contains all the necessary scripts to power this documentation.

The folder `e_cascading_data/` is where the cascading data produced by the simulations are saved for future use by the program itself. The cascading data are saved by species in the subfolders `N2/`, `O2/` and `O/`. You should not need to venture into this folder.

The folder `scripts/` contains the scripts for the user to start the simulations and plot some of the results.

The folder `src/` contains the source code of the model.
