import json
import os
from os.path import exists
import subprocess
import re
import random

def getConstAtomPairIndexList(constEnzAtom, constSubAtom):

    constAtomPairIndexList = list()
    with open("COM/com_dry.pdb") as f:
        for line in f.readlines():
            atomNum = line[7:11].replace(' ', '')
            atomType = line[13:16].replace(' ', '')
            resName = line[17:20].replace(' ', '')
            resNum = line[23:26].replace(' ', '')

            if resName == 'SUB':
                if atomType == constSubAtom:
                    constAtomPairIndexList.append(int(atomNum) - 1)
            else:
                resNameNumAtom_enz = f"{resName}{resNum}@{atomType}"
                if resNameNumAtom_enz == constEnzAtom:     
                    constAtomPairIndexList.append(int(atomNum) - 1)
    return constAtomPairIndexList

def get_reslen(rec_pdb):
    cmd  = f"pdb_wc -r {rec_pdb} | grep residues | cut -f2"
    tot_res_pdb = int(os.popen(cmd).read().strip())
    return tot_res_pdb

"""
0         1         2
0123456789012345678901234567890123456789
ATOM      1  N   PRO     1       0.943  -2.086  -0.237  1.00  0.00
ATOM      2  CD  PRO     1       0.158  -3.126  -0.850  1.00  0.00
ATOM      3  HD2 PRO     1      -0.550  -2.704  -1.563  1.00  0.00
ATOM      4  HD3 PRO     1      -0.387  -3.660  -0.071  1.00  0.00
ATOM      5  CG  PRO     1       1.198  -3.982  -1.510  1.00  0.00
ATOM      6  HG2 PRO     1       0.814  -4.354  -2.460  1.00  0.00
ATOM      7  HG3 PRO     1       1.440  -4.825  -0.862  1.00  0.00
ATOM      8  CB  PRO     1       2.402  -3.094  -1.711  1.00  0.00
ATOM      9  HB2 PRO     1       2.464  -2.797  -2.758  1.00  0.00
ATOM     10  HB3 PRO     1       3.305  -3.637  -1.434  1.00  0.00
....
ATOM     61  HB3 ASP     5      -1.658  -4.338   0.737  1.00  0.00
ATOM     62  CG  ASP     5      -0.092  -5.269   1.882  1.00  0.00
ATOM     63  OD1 ASP     5       0.190  -5.570   3.072  1.00  0.00
ATOM     64  OD2 ASP     5       0.399  -5.806   0.853  1.00  0.00
ATOM     65  C   ASP     5       0.604  -2.885   0.449  1.00  0.00
ATOM     66  O   ASP     5       1.658  -3.512   0.581  1.00  0.00
TER   
ATOM     67 MG   MG      6       0.000   0.000   0.000  1.00  0.00
TER   
END   
"""

def run_openmm(inputD):

    properties_line = "properties = {'CudaPrecision': 'mixed'}"

    ######################################################
    # nano E-9; pico E-12; femto E-15
    # ----------------------------------------------------
    # DEFAULT TIME-SCALE
    #   md_duration_in_ns = 0.1 ns
    #   time_step_in_fs   = 2 fs
    #   md_steps          = 1E(-10+15)*(md_duration_in_ns / time_step_in_fs)  
    #                     = 1E+5     *(md_duration_in_ns / time_step_in_fs)
    #                     = 1E+5     *(1.0/2.0) 
    #                     = 5E+4 = 50,000
    #   snapshot_interval = 2ps 
    #                     = 2E-12
    #                     = 2E(+3-15) 
    #                     = (1E+3 steps) * 2E-15
    #                     = snapshot_interval_steps_fixed_2ps * time_step_in_fs
    #   snapshot_interval_steps_fixed_2ps = 2ps / time_step_in_fs
    #                                     = 2E-9 / time_step_in_fs
    #                                     = 2E(-9+15) / time_step_in_fs
    #                                     = 2E+3 / time_step_in_fs
    #                                     = 2E+3 / 2(ps)
    #                                     = 1,000
    ######################################################
    # md_steps = 1E+6 * (md_duration_in_ns / time_step_in_fs) 
    # snapshot_interval_steps_fixed_2ps = 2E+3 / (time_step_in_fs)

    md_steps = int( 1E+6 * float(inputD['md_duration_in_ns'])/float(inputD['time_step_in_fs']) )
    snapshot_interval_steps_fixed_2ps =  int( 2E+3 / float(inputD['time_step_in_fs']))
    print()
    print("md_duration_in_ns", inputD['md_duration_in_ns'])
    print("time_step_in_fs", inputD['time_step_in_fs'])
    print("md_steps = 1E+6 (md_duration_in_ns/time_step_in_fs) = ", md_steps)
    print(f"snapshot_interval_steps_fixed_2ps = 2E+3 * time_step_in_fs = ", snapshot_interval_steps_fixed_2ps)
    print()
    
    time_step_in_fs_line = f" {inputD['time_step_in_fs']}*femtoseconds)\n"
    dcd_line = f"simulation.reporters.append(DCDReporter('com_wat.dcd', {snapshot_interval_steps_fixed_2ps}))\n"
    state_line = f"simulation.reporters.append(StateDataReporter('openMM.log', {snapshot_interval_steps_fixed_2ps}, step=True, potentialEnergy=True, temperature=True, separator=','))\n"
    #simulation.reporters.append(DCDReporter('com_wat.dcd', 1000))
    #simulation.reporters.append(StateDataReporter('openMM.log', 1000, step=True, potentialEnergy=True, temperature=True, separator=','))

    if inputD['option'] == "constMD":
        const_enz_atom = inputD["const_enz_atom"]
        const_sub_atom = inputD["const_sub_atom"]
        constAtomPairIndexList = getConstAtomPairIndexList(const_enz_atom, const_sub_atom)
        constraints_line = f"system.addConstraint({constAtomPairIndexList[0]}, {constAtomPairIndexList[1]}, 2.1*angstroms)"
    else:
        constraints_line = "# no constraint"


    input_openmm = f"""from __future__ import print_function
import openmm as mm
from openmm import *
from openmm import app
from openmm.app import *
from sys import stdout
from simtk.unit import *

# load in Amber input files
prmtop = AmberPrmtopFile('COM/com_wat.prmtop')
inpcrd = AmberInpcrdFile('COM/com_wat.inpcrd')

# prepare system and integrator
system = prmtop.createSystem(nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True,
    ewaldErrorTolerance=0.0005)
integrator = mm.LangevinIntegrator(300*kelvin, 1.0/picoseconds,
{time_step_in_fs_line}
#   2.0*femtoseconds) # default
integrator.setConstraintTolerance(0.00001)
integrator.setRandomNumberSeed({inputD['random_seed']})

# prepare simulation
platform = mm.Platform.getPlatformByName('CUDA')
{properties_line}

#######################################################
# constMD or normalMD
{constraints_line}
#######################################################
simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
simulation.context.setPositions(inpcrd.positions)

# minimize
print('Minimizing...')
simulation.minimizeEnergy()

print('Saving Minimzed PDB...')
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('tmp.pdb', 'w'))
###########################

# equilibrate for 100 steps
simulation.context.setVelocitiesToTemperature(300*kelvin)
print('Equilibrating...')
simulation.step(100)


# append reporters
#simulation.reporters.append(DCDReporter('com_wat.dcd', 1000))
#simulation.reporters.append(StateDataReporter('openMM.log', 1000, step=True, potentialEnergy=True, temperature=True, separator=','))
{dcd_line}
{state_line}

# run production simulation
print('Running Production...')
simulation.step({md_steps})

print('Done!')
""" 
    with open("openMM.py", 'w') as f:
        f.write(input_openmm)

    # run
    proc = subprocess.run(f"python openMM.py > openMM.log", shell=True)    
    if proc.returncode == 0:
        print(f"openMM Simulation is done")
        print()

def save_dry_traj_ensemble_pdbs():
    if not os.path.exists('RESULT'):
        os.makedirs('RESULT')

    # save_dry_emin_pdb
    subprocess.run(f'grep -e ATOM tmp.pdb -e SUB > RESULT/Emin_com_dry.pdb', shell=True)

    # save dry trajectory & ensemble pdbs
    reslen = get_reslen("COM/com_dry.pdb")
    cpptraj_input = f"""
cpptraj -p  COM/com_wat.prmtop << EOF

trajin com_wat.dcd

# Centeing snapshots
autoimage
center :1-{reslen} mass origin
image origin center familiar com :1-{reslen}

# Strip water and ions
strip :WAT
strip :Na+
strip :Cl-

# stripped water/ion trajectory
trajout RESULT/com_dry.crd nobox
trajout RESULT/com_dry_ensemble.pdb pdb

go
quit
EOF
"""
    subprocess.run(cpptraj_input, shell=True)


def trajectory_analysis():
    if not os.path.exists('RESULT/com_dry.prmtop'):
        os.system('cp COM/com_dry.prmtop RESULT/.')

    reslen = get_reslen("COM/com_dry.pdb")
    cpptraj_input = f"""
cpptraj -p  RESULT/com_dry.prmtop << EOF
trajin RESULT/com_dry.crd
reference COM/com_dry.inpcrd

# rmsd
rms reference out bb_rmsd.out :1-{reslen}@C,CA,N,O
rms reference out all_rmsd.out :1-{reslen}
atomicfluct out RESULT/rms_byres.out byres

go
quit
EOF
"""
    subprocess.run(cpptraj_input, shell=True)

def step2_Simulation(inputD):
    print(f'\nStep2: openMM Simulation is running...(duration : {inputD["md_duration_in_ns"]}ns / save every 2ps)\n')
    run_openmm(inputD)
    save_dry_traj_ensemble_pdbs()
    trajectory_analysis()
    print('Done\n')
    print('######################')
