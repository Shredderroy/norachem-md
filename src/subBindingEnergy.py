import json
import os
from os.path import exists
import subprocess
import pandas as pd

def get_reslen(rec_pdb):
    cmd  = f"pdb_wc -r {rec_pdb} | grep residues | cut -f2"
    tot_res_pdb = int(os.popen(cmd).read().strip())
    return tot_res_pdb

def ante_mmpbsa():
    if exists('comD.prmtop'):
        os.system('rm comD.prmtop recD.prmtop subD.prmtop')
    com_len = get_reslen('COM/com_dry.pdb')
    rec_len = com_len-1
    cmd = f'ante-MMPBSA.py -p COM/com_wat.prmtop \
-s "!(:1-{com_len})" -c comD.prmtop \
-m ":1-{rec_len}" -r recD.prmtop -l subD.prmtop'
    print(cmd)
    subprocess.run(cmd, shell=True)    

def mmpbsa_input_gen():
    mmpbsa_input = """
Input file for running GB in serial
&general
   interval = 1,
   entropy = 0,
   keep_files = 0,
/
&gb
  igb=2,
  saltcon=0.100,
/
"""
    with open("MMPBSA.in", 'w') as f:
        f.write(mmpbsa_input)

    mmpbsa_run_cmd = """MMPBSA.py -O \
-i MMPBSA.in \
-lp  subD.prmtop -rp  recD.prmtop -cp comD.prmtop \
-y   RESULT/com_dry.crd \
-o   FINAL_RESULTS_MMPBSA_ave.dat  \
-eo  FINAL_RESULTS_MMPBSA_perframe.dat 
"""
    with open("MMPBSA.sh", 'w') as f:
        f.write(mmpbsa_run_cmd)


def getRmsd():
    bb_cmd = "cat bb_rmsd.out | grep -v RMSD | awk '{print $2}'"
    bb_list = os.popen(bb_cmd).read().strip()
    bb_finalList = [float(m) for m in bb_list.split('\n')]

    all_cmd = "cat all_rmsd.out | grep -v RMSD | awk '{print $2}'"
    all_list = os.popen(all_cmd).read().strip()
    all_finalList = [float(m) for m in all_list.split('\n')]
    return bb_finalList, all_finalList


def Esummary():
    if exists('reference.frc'):
        os.system('rm reference.frc')
    with open(f"FINAL_RESULTS_MMPBSA_perframe.dat", 'r') as f:
        lines = f.readlines()
        myDelta = False
        with open(f"tmp.csv", 'w') as g:
            for line in lines:
                if line.startswith("DELTA Energy Terms"):
                    myDelta = True
                if myDelta and not line.startswith("DELTA Energy Terms"):
                    g.write(line)
    df = pd.read_csv('tmp.csv')
    df = df[['Frame #','DELTA TOTAL']]
    df.columns = ['Frame', 'substate_bindingE']
    df['Time_nanosecond']  = [(i+1)*2E-3 for i in df['Frame']]

    #cmd = "cat rmsd.out | grep -v RMSD | awk '{print $2}'"
    #subprocess.run(cmd, shell=True)
    df["bb_rmsd"], df["all_rmsd"] = getRmsd()

    df.to_csv('RESULT/Esummary_Rmsd.csv', index=False)
    if exists('tmp.csv'):
        os.system('rm tmp.csv')
    

def step3_Energy():
    print('Step3: Substrate Binding Energy Calculation')
    ante_mmpbsa()
    mmpbsa_input_gen()
    subprocess.run(f"bash MMPBSA.sh", shell=True)    
    Esummary()
    print('Done: see Esummary_Rmsd.csv')
    print('######################')
    print()


