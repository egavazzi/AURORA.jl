var documenterSearchIndex = {"docs":
[{"location":"manual_folder-structure/#Folder-structure","page":"Folder structure","title":"Folder structure","text":"","category":"section"},{"location":"manual_folder-structure/","page":"Folder structure","title":"Folder structure","text":"The code is structured as follow","category":"page"},{"location":"manual_folder-structure/","page":"Folder structure","title":"Folder structure","text":"AURORA/\n├── data/\n│   └── 20220926/...\n│   └── ...\n│\n├── docs/...\n│ \n├── e_cascading_data/\n│   └── N2/...\n│   └── O2/...\n│   └── O/...\n│\n├── scripts/\n│   └── Control_script.jl\n│   └── ...\n│ \n└── src/\n    └── main.jl\n    └── cascading.jl\n    └── setup.jl\n    └── ...","category":"page"},{"location":"manual_folder-structure/","page":"Folder structure","title":"Folder structure","text":"The folder data/ contains the subfolders where simulation results are saved.","category":"page"},{"location":"manual_folder-structure/","page":"Folder structure","title":"Folder structure","text":"The folder docs/ contains all the necessary scripts to power this documentation.","category":"page"},{"location":"manual_folder-structure/","page":"Folder structure","title":"Folder structure","text":"The folder e_cascading_data/ is where the cascading data produced by the simulations are saved for future use by the program itself. The cascading data are saved by species in the subfolders N2/, O2/ and O/. You should not need to venture into this folder.","category":"page"},{"location":"manual_folder-structure/","page":"Folder structure","title":"Folder structure","text":"The folder scripts/ contains the scripts for the user to start the simulations and plot some of the results.","category":"page"},{"location":"manual_folder-structure/","page":"Folder structure","title":"Folder structure","text":"The folder src/ contains the source code of the model.","category":"page"},{"location":"manual_get-started/#Get-started","page":"Get started","title":"Get started","text":"","category":"section"},{"location":"manual_get-started/","page":"Get started","title":"Get started","text":"Simulations are started using the function calculate_e_transport(...). The function takes in many parameters, so it can be easier to use the script named Control_script.jl situated in the scripts/ folder. The script is pre-filled and you just need to modify the values of the parameters. You can also use the script as a template to make your own control scripts.","category":"page"},{"location":"manual_get-started/","page":"Get started","title":"Get started","text":"After you have modified the Control_script.jl and saved it, you just have to  execute it. This can be done from the Julia REPL or from the command line. ","category":"page"},{"location":"manual_get-started/#Starting-simulation-from-the-Julia-REPL","page":"Get started","title":"Starting simulation from the Julia REPL","text":"","category":"section"},{"location":"manual_get-started/","page":"Get started","title":"Get started","text":"warning: Activating the AURORA environment\nTo be able to use AURORA.jl, the repository environment needs to be activated. This can be done for example by starting Julia from the AURORA.jl/ folder using the command$> julia --project=.Or by entering the Pkg REPL using the ] key and typingpkg> activate .","category":"page"},{"location":"manual_get-started/","page":"Get started","title":"Get started","text":"info: Using VS Code\nIf you are using VS Code with the Julia extension, the local environment should be automatically activated when you open the AURORA.jl/ folder.","category":"page"},{"location":"manual_get-started/","page":"Get started","title":"Get started","text":"Once the AURORA.jl environment activated, you can start simulations from the Julia REPL with the command","category":"page"},{"location":"manual_get-started/","page":"Get started","title":"Get started","text":"julia> include(\"scripts/Control_script.jl\")","category":"page"},{"location":"manual_get-started/","page":"Get started","title":"Get started","text":"If you are using VS Code, you can also use the \"Execute active File in REPL\" button.","category":"page"},{"location":"manual_get-started/","page":"Get started","title":"Get started","text":"The results will be saved in a folder under data/ along with the parameters used to run the simulation.","category":"page"},{"location":"manual_get-started/#Starting-simulation-from-the-command-line","page":"Get started","title":"Starting simulation from the command line","text":"","category":"section"},{"location":"manual_get-started/","page":"Get started","title":"Get started","text":"Move to the AURORA.jl/ folder. Then, execute the Control_script.jl using the command ","category":"page"},{"location":"manual_get-started/","page":"Get started","title":"Get started","text":"$> julia --project=@. scripts/Control_script.jl ","category":"page"},{"location":"manual_get-started/","page":"Get started","title":"Get started","text":"The results will be saved in a folder under data/ along with the parameters used to run the simulation.","category":"page"},{"location":"manual/#Installation","page":"-","title":"Installation","text":"","category":"section"},{"location":"manual/","page":"-","title":"-","text":"Clone the repository (e.g. with Git) or download and extract the .zip (available under the green code button on the GitHub page)\nOpen a terminal, move into the AURORA.jl, and start Julia with the command","category":"page"},{"location":"manual/","page":"-","title":"-","text":"$> julia","category":"page"},{"location":"manual/","page":"-","title":"-","text":"Then, activate the repository and install the packages required by AURORA.jl using the commands","category":"page"},{"location":"manual/","page":"-","title":"-","text":"julia> using Pkg\njulia> Pkg.activate(\".\")\njulia> Pkg.instantiate() # this might take a while...","category":"page"},{"location":"manual/","page":"-","title":"-","text":"AURORA.jl is now ready to use!","category":"page"},{"location":"manual/#Folder-structure","page":"-","title":"Folder structure","text":"","category":"section"},{"location":"manual/","page":"-","title":"-","text":"The code is structured as follow","category":"page"},{"location":"manual/","page":"-","title":"-","text":"AURORA/\n├── data/\n│   └── 20220926/...\n│   └── ...\n│\n├── docs/...\n│ \n├── e_cascading_data/\n│   └── N2/...\n│   └── O2/...\n│   └── O/...\n│\n├── scripts/\n│   └── Control_script.jl\n│   └── ...\n│ \n└── src/\n    └── main.jl\n    └── cascading.jl\n    └── setup.jl\n    └── ...","category":"page"},{"location":"manual/","page":"-","title":"-","text":"The folder data/ contains the subfolders where simulation results are saved.","category":"page"},{"location":"manual/","page":"-","title":"-","text":"The folder docs/ contains all the necessary scripts to power this documentation.","category":"page"},{"location":"manual/","page":"-","title":"-","text":"The folder e_cascading_data/ is where the cascading data produced by the simulations are saved for future use by the program itself. The cascading data are saved by species in the subfolders N2/, O2/ and O/. You should not need to venture into this folder.","category":"page"},{"location":"manual/","page":"-","title":"-","text":"The folder scripts/ contains the scripts for the user to start the simulations and plot some of the results.","category":"page"},{"location":"manual/","page":"-","title":"-","text":"The folder src/ contains the source code of the model.","category":"page"},{"location":"manual/#Get-started","page":"-","title":"Get started","text":"","category":"section"},{"location":"manual/","page":"-","title":"-","text":"Simulations are started using the function calculate_e_transport(...). The function takes in many parameters, so it can be easier to use the script named Control_script.jl situated in the scripts/ folder. The script is pre-filled and you just need to modify the values of the parameters. You can also use the script as a template to make your own control scripts.","category":"page"},{"location":"manual/","page":"-","title":"-","text":"After you have modified the Control_script.jl and saved it, you just have to  execute it. This can be done from the Julia REPL or from the command line. ","category":"page"},{"location":"manual/#Starting-simulation-from-the-Julia-REPL","page":"-","title":"Starting simulation from the Julia REPL","text":"","category":"section"},{"location":"manual/","page":"-","title":"-","text":"warning: Activating the AURORA environment\nTo be able to use AURORA.jl, the repository environment needs to be activated. This can be done for example by starting Julia from the AURORA.jl/ folder using the command$> julia --project=.Or by entering the Pkg REPL using the ] key and typingpkg> activate .","category":"page"},{"location":"manual/","page":"-","title":"-","text":"info: Using VS Code\nIf you are using VS Code with the Julia extension, the local environment should be automatically activated when you open the AURORA.jl/ folder.","category":"page"},{"location":"manual/","page":"-","title":"-","text":"Once the AURORA.jl environment activated, you can start simulations from the Julia REPL with the command","category":"page"},{"location":"manual/","page":"-","title":"-","text":"julia> include(\"scripts/Control_script.jl\")","category":"page"},{"location":"manual/","page":"-","title":"-","text":"If you are using VS Code, you can also use the \"Execute active File in REPL\" button.","category":"page"},{"location":"manual/","page":"-","title":"-","text":"The results will be saved in a folder under data/ along with the parameters used to run the simulation.","category":"page"},{"location":"manual/#Starting-simulation-from-the-command-line","page":"-","title":"Starting simulation from the command line","text":"","category":"section"},{"location":"manual/","page":"-","title":"-","text":"Move to the AURORA.jl/ folder. Then, execute the Control_script.jl using the command ","category":"page"},{"location":"manual/","page":"-","title":"-","text":"$> julia --project=@. scripts/Control_script.jl ","category":"page"},{"location":"manual/","page":"-","title":"-","text":"The results will be saved in a folder under data/ along with the parameters used to run the simulation.","category":"page"},{"location":"manual_installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"manual_installation/","page":"Installation","title":"Installation","text":"Clone the repository (e.g. with Git) or download and extract the .zip (available under the green code button on the GitHub page)\nOpen a terminal, move into the AURORA.jl, and start Julia with the command","category":"page"},{"location":"manual_installation/","page":"Installation","title":"Installation","text":"$> julia","category":"page"},{"location":"manual_installation/","page":"Installation","title":"Installation","text":"Then, activate the repository and install the packages required by AURORA.jl using the commands","category":"page"},{"location":"manual_installation/","page":"Installation","title":"Installation","text":"julia> using Pkg\njulia> Pkg.activate(\".\")\njulia> Pkg.instantiate() # this might take a while...","category":"page"},{"location":"manual_installation/","page":"Installation","title":"Installation","text":"AURORA.jl is now ready to use!","category":"page"},{"location":"#This-is-the-AURORA.jl-documentation","page":"Home","title":"This is the AURORA.jl documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Work in progress\nThis project is still in progress. Be aware that the documentation might not be up to date.","category":"page"}]
}