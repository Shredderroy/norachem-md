## Dockerized run (NOTE: NVIDIA-DOCKER!!! NOT just Docker)
to build docker container
```
$ ./build.bat 
```

To use dockerized verions in normalMD_example folder run 
```
$ nvidia-docker run --rm -v $(pwd)/:/eval/ -w="/eval/"  norachem-md:v0
```

## How to Run on Conda (prior to Docker)
```
[/]$ conda env create -f environment.yml
[/]$ conda activate SIM
(SIM)[/]$ cd constMD_example
(SIM)[constMD_example/]$ python ../src/SIM.py

```
## Output Result Files
```

normalMD_example/
    ├── COM
    ├── SUB
    └── RESULT
        ├── com_dry.crd                 # MD trajectory (use Chimera* to visualize)                                           
        ├── com_dry.prmtop              # MD parameter file (use Chimera* to visualize)
        ├── com_dry_ensemble.pdb        # MD snapshots all in pdb
        ├── Esummary_Rmsd.csv           # MD analysis (time_in_nanosecond, substrate_binding_energy, Backbone_RMSD, AllAtoms_RMSD)
        ├── Emin_com_dry.pdb            # Energy minimized sructure in pdb
        └── rms_byres.out               # Per Residue RMSD
        * https://sites.engineering.ucsb.edu/~shell/che210d/Visualization.pdf
```
