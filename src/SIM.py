import json
import os
from os.path import exists
import subprocess
import re
import random
from prepSystems import *
from openmmSimulations import *
from subBindingEnergy import *

optionD = {1:"normalMD"}

""" input.json 

{
    "random_seed": 42,
    "md_duration_in_ns": 0.1,              
    "time_step_in_fs": 2.0, 
    "target_pdb": "enzyme.pdb",           
    "small_molecule_pdb": "res.pdb",
}
"""

def enzyme_resNum_atomList(enzFile):
    resNum_atomList = dict()
    with open(enzFile, 'r') as f:
        for line in f.readlines():
            atomType = line[13:16].replace(' ', '')
            resName = line[17:20].replace(' ', '')
            resNum = line[23:26].replace(' ', '')
            key = f"{resName}{resNum}"
            if not resNum_atomList.get(key):
                resNum_atomList[key] = list()
            resNum_atomList[key].append(atomType)
    return resNum_atomList
    # enzD = {'LEU402':['HG', CD1', ...}

def substrate_atomList(subFile):
    subAtomList = list()
    with open(subFile, 'r') as f:
        for line in f.readlines():
            atomType = line[13:16].replace(' ', '')
            subAtomList.append(atomType)
    return subAtomList
    # subAtomList = ['HG', CD1', ...]
            
def check_enzyme_catRes(enzFile, catResList):
    enzD = enzyme_resNum_atomList(enzFile)
    resList = enzD.keys()
    for catRes in catResList.split(','):
        if not catRes in resList:
            return False
    return True

def check_enzyme_constAtom(enzFile, constResAtom):
    enzD = enzyme_resNum_atomList(enzFile)
    resList = enzD.keys()
    constRes, constAtom = constResAtom.split('@')[:]
    if not constRes in resList:
        return False
    if not constAtom in enzD[constRes]:
        return False
    return True
    
def check_substrate_constAtom(subFile, constSubAtom):
    subAtomList = substrate_atomList(subFile)
    if not constSubAtom in subAtomList:
        return False
    return True

def assert_option(data):
         
    args_list = data
    inputD = dict()


    enz_pdb = args_list['target_pdb']
    if not exists(enz_pdb):
        print(f"\nError: can't find {target_pdb}\n")
        return False

    sub_pdb = args_list['small_molecule_pdb']
    if not exists(sub_pdb):
        print(f"\nError: can't find {small_molecule_pdb}\n")
        return False
    
    const_enz_atom = ""
    const_sub_atom =  ""

    inputD['option'] = "normalMD"
    inputD['md_duration_in_ns'] = args_list['md_duration_in_ns']
    inputD['time_step_in_fs'] = args_list['time_step_in_fs']
    inputD['enz_pdb'] = enz_pdb
    inputD['sub_pdb'] = sub_pdb
    inputD['const_enz_atom'] = const_enz_atom
    inputD['const_sub_atom'] = const_sub_atom
    inputD['random_seed'] = args_list['random_seed']
    inputD['substrate_net_charge'] = args_list['substrate_net_charge']


    ##############
    return inputD
    #   => [option, md_duration_in_ns, time_step_in_fs, enz_pdb, sub_pdb, const_enz_atom, const_sub_atom, random_seed, substrate_net_charge]


def assert_input_json():
    with open('input.json', 'r') as f:
        data = json.load(f)
    
    # options = [data[f"{opt}"]["run_this"] for opt in optionD.values()]

    # if sum(options) > 1:
    #     print()
    #     print("Error: must select only one option ")
    #     print("       by setting run_this field 1 in input.json")
    #     print()
    #     return False 

    # sel_option = options.index(1)
    inputOK = assert_option(data)
    return inputOK



if __name__ == '__main__':

    # check input.json
    inputD = assert_input_json()
    print(inputD)
    step1_Prep(inputD)
    step2_Simulation(inputD)
    step3_Energy()
