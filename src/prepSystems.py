import subprocess
import sys
import os
from os.path import exists
import time


#{'option': 'constMD', 'md_duration_in_ns': 1, 'enz_pdb': 'enzyme.pdb', 'sub_pdb': 'res.pdb', 
# 'const_enz_atom': 'ASP170@OD2', 'const_sub_atom': 'C3'}


aa31 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
        'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

aa13 = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS', 
        'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
        'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
        'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

def substrate_residue_rename_to_SUB(docked_substrate_pdb):
    with open(docked_substrate_pdb, 'r') as f:
        lines = f.readlines()
    with open("SUB.pdb", 'w') as f:
        for line in lines:
            if line[12:13] != ' ':
                f.write(line[:12] + " " + line[12:13] + line[13:14].lower() + line[14:15] + " " + "SUB" + line[20:])
            elif line.startswith('CONECT'):
                continue
            else:
                f.write(line[:17] + "SUB" + line[20:])

def prep_substrate_params(docked_substrate_pdb, net_charge):
    # mkdir and chdir to SUB
    os.makedirs('SUB', exist_ok=True)
    os.chdir('SUB')

    # rename substrate residue to 'SUB'
    substrate_residue_rename_to_SUB('../' + docked_substrate_pdb)

    if exists('leap.log'):
        subprocess.run('rm leap.log', shell=True)

    sub_input = f"""
# step1_antechamber.sh
# input = docked_substrate_pdb
# output = SUB.prepi
antechamber -fi pdb -fo prepi -i SUB.pdb -o SUB.prepi -rn SUB -c bcc -pf y -nc {net_charge} -at amber

# step2_parmchk2.sh
# input = SUB.prepi
# output = SUB.frcmod
parmchk2 -f prepi -i SUB.prepi -o SUB.frcmod 

# step3_tleap.sh
tleap -f - <<EOF
source leaprc.gaff #Source leaprc file for gaff
source leaprc.water.tip3p 
# SUB
SUB = loadpdb SUB.pdb
charge SUB
loadamberprep SUB.prepi 
loadamberparams SUB.frcmod 
saveamberparm SUB SUB.prmtop SUB.inpcrd # SUB-dry
solvateBox SUB TIP3PBOX 15.0
savepdb SUB SUB_wat.pdb
saveamberparm SUB SUB_wat.prmtop SUB_wat.inpcrd
quit
EOF

#step4_ambpdbChk.sh
ambpdb -p SUB.prmtop < SUB.inpcrd > SUB_ambpdb.pdb
"""
    subprocess.run(sub_input, shell=True)

    prep_sub_OK = False
    with open('leap.log', 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("Exiting LEaP"):
                if line.find("Errors = 0") > 0:
                    prep_sub_OK = True
    os.chdir('../')
    if not prep_sub_OK:
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("Error in making substrate paramters files")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        return False

    print("#################################################")
    print("Success: making substrate paramters files in SUB/")
    print("#################################################")
    return True



def make_complex_pdb(enzFile):

    # enzyme: del H, change HIS -> HIE
    enz_input = f"""cat ../{enzFile} | 
grep -e 'ATOM' | \
grep -v '   H' | \
sed "s/HIS/HIE/g" | \
grep -v "^END" > enzyme_noH.pdb
echo TER >> enzyme_noH.pdb
"""
    subprocess.run(enz_input, shell=True)

    # concat enzyme_noH.pdb SUB.pdb
    subprocess.run("cat enzyme_noH.pdb ../SUB/SUB.pdb > enzyme_sub.pdb", shell=True)

    if exists('leap.log'):
        subprocess.run('rm leap.log', shell=True)

    # check charge
    complex_charge = f"""
tleap -f - <<EOF
source leaprc.protein.ff14SB 
source leaprc.gaff #Source leaprc file for gaff

loadAmberPrep    ../SUB/SUB.prepi
loadAmberParams  ../SUB/SUB.frcmod

# complexM
complexM = loadpdb enzyme_sub.pdb
charge complexM

quit ##### <<< KEY
EOF
"""
    subprocess.run(complex_charge, shell=True)
    
    with open('leap.log', 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("Total unperturbed charge"):
                charge = float(line.split(':')[1].replace(" ", ""))
                break
    if charge > 0:
        return 'Cl-'
    else:
        return 'Na+'
            
        
    
def prep_complex_params(enzFile):
    # mkdir and chdir to COM
    os.makedirs('COM', exist_ok=True)
    os.chdir('COM')

    # make COM/enzyem_sub.pdb
    # compute overall charge of complex
    # and return counter ion 
    counter_ion = make_complex_pdb(enzFile)
    
    if exists('leap.log'):
        subprocess.run('rm leap.log', shell=True)

    complex_input = f"""
# step1_tleap_complex.in
tleap -f - <<EOF
source leaprc.protein.ff14SB 
source leaprc.water.tip3p 
source leaprc.gaff #Source leaprc file for gaff

loadAmberPrep    ../SUB/SUB.prepi
loadAmberParams  ../SUB/SUB.frcmod

# complexM
complexM = loadpdb enzyme_sub.pdb
saveamberparm complexM com_dry.prmtop com_dry.inpcrd # complexM-dry  
savepdb complexM com_dry.pdb
charge complexM
addions complexM {counter_ion} 0
solvateOct complexM TIP3PBOX 15.0
savepdb complexM com_wat.pdb
saveamberparm complexM com_wat.prmtop com_wat.inpcrd

quit ##### <<< KEY
EOF
#
"""
    subprocess.run(complex_input, shell=True)

    prep_complex_OK = False
    with open('leap.log', 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("Exiting LEaP"):
                if line.find("Errors = 0") > 0:
                    prep_complex_OK = True
    os.chdir('../')
    if not prep_complex_OK:
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("Error in making complex paramters files")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        return False

    print("#################################################")
    print("Success: making complex paramters files in COM/")
    print("#################################################")
    return True

def prep_md_parameter_files(enzFile, subFile, net_charge):
    # prep substrate parameter files in SUB/
    #sub_OK = prep_substrate_params(subFile)
    sub_OK = prep_substrate_params(subFile, net_charge)
    if not sub_OK:
        return False

    # prep complex parameter files in COM/
    com_OK = prep_complex_params(enzFile)
    if not com_OK:        
        return False

    return True

def step1_Prep(inputD):
    print()
    print('====================================================')
    for k, v in inputD.items():
        print(f"{k:<{20}} : {v}")
    print('====================================================')
    print()

    print("\nStep1: Preparation of parameter files from PDBs\n")
    enz_pdb = inputD['enz_pdb']
    sub_pdb = inputD['sub_pdb']
    net_charge = inputD['substrate_net_charge']
    prep_md_parameter_files(enz_pdb, sub_pdb, net_charge)
    print('Done\n')
    print('######################')
