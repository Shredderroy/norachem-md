## SIM
*  *AmberTools* for Parameter Preparation + Trajectory Post-Processing
*  *OpenMM* for Molecular Dynamics Simulation 
*  *MMPBSA* for Substrate Binding Energy Calculation

## Input.json
```
{
    "constMD": {
        "run_this": 1,          #######  either constraint MD
        "random_seed": 42,
        "md_duration_in_ns": 0.1, 
        "time_step_in_fs": 2.0, 
        "enzyme_pdb": "enzyme.pdb",             
        "docked_substrate_pdb": "res.pdb",   
        "substrate_net_charge": 0,
        "constraint_enzyme_atom": "ASP170@OD2",
        "constraint_substrate_atom": "C3"
    },
    "normalMD": {
        "run_this": 0,          #######  or normal MD w/o constraint
        "random_seed": 42,
        "md_duration_in_ns": 0.1, 
        "time_step_in_fs": 2.0, 
        "enzyme_pdb": "enzyme.pdb",           
        "docked_substrate_pdb": "res.pdb",
        "substrate_net_charge": 0
    }
}
```

## Input Examples Files
```
constMD_example
    ├── enzyme.pdb
    ├── input.json
    └── res.pdb    
normalMD_example
    ├── enzyme.pdb
    ├── input.json
    └── res.pdb
```

## Dockerized run (NOTE: NVIDIA-DOCKER!!! NOT just Docker)
```
$ aws ecr get-login-password --region us-west-2 | docker login --username AWS --password-stdin 015042944066.dkr.ecr.us-west-2.amazonaws.com
$ nvidia-docker run --rm -v $(pwd)/:/eval/ -w="/eval/"  015042944066.dkr.ecr.us-west-2.amazonaws.com/openmmsim/openmmsim:v1.0
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
constMD_example_output/
    ├── COM/                # generated complex = (enzyme + substrate) parameters            
    ├── SUB/                # generated substrate parameters
    └── RESULT/             # final outputs
        ├── com_dry.crd                 # MD trajectory (use Chimera* to visualize)                                           
        ├── com_dry.prmtop              # MD parameter file (use Chimera* to visualize)
        ├── com_dry_ensemble.pdb        # MD snapshots all in pdb
        ├── Esummary_Rmsd.csv           # MD analysis (time_in_nanosecond, substrate_binding_energy, Backbone_RMSD, AllAtoms_RMSD)
        ├── Emin_com_dry.pdb            # Energy minimized sructure in pdb
        └── rms_byres.out               # Per Residue RMSD
        * https://sites.engineering.ucsb.edu/~shell/che210d/Visualization.pdf

normalMD_example_output/
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
To use dockerized verions in normalMD_example folder run 
```
nvidia-docker run --rm -v $(pwd)/:/eval/ -w="/eval/"  openmmsim:v0.1
```
