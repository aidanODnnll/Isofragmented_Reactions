# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 19:04:52 2020

@author: bol216
"""

from rdkit import Chem
from model.training_model_GA import *
from numpy.linalg import inv
import numpy.linalg as lin


def H_count(SMILES):
    m = Chem.MolFromSmiles(SMILES)
    num1 = m.GetNumAtoms(True)
    mh = Chem.AddHs(m)
    num2 = mh.GetNumAtoms(True)
    H = num2 - num1
    return H
    
def atom_count(SMILES,pattern):
    name = SMILES
    mol = Chem.MolFromSmiles(name)
    rd_count = []
    for sub in pattern:
        if sub != 'H':  
            sub = sub.replace('R','')
            regular = len(mol.GetSubstructMatches(Chem.MolFromSmarts(sub)))
            aroma = len(mol.GetSubstructMatches(Chem.MolFromSmarts(sub.lower())))
            rd_count.append(regular + aroma)
    rd_count.append(H_count(SMILES))    
    return rd_count    

def loading_bonds(mode = 'regular'):
    if mode == 'regular':    
        pattern = ['CO','C=O','CC','C=C','C#C','C[H]','c[H]','cC','O[H]','cc','cO','co']
    elif mode == 'hybridization':
        pattern = ['[C^3][C^3]', '[C^3][C^2]', '[C^2][C^2]', '[C^1][C^1]', '[C^2][C^1]', '[C^3][C^1]', 'c[C^3]', 'c[C^2]', 'c[C^1]', 'cc', '[C^2]=[C^2]', '[C^1]#[C^1]', '[C^3][H]', '[C^2][H]', '[C^1][H]', 'c[H]', '[C^3][O^3]', '[C^2][O^3]', 'c[O^3]', 'co', '[C^2]=[O^2]', '[O^3][H]']
    elif mode == 'pathway fps':
        pattern = []
        f = open('COH_patterns_fps.txt','r').readlines()
        for line in f:
            lis = line.split()
            pattern.append(lis[0])
    else:
        print('the type of '+mode+' is not included')
    return pattern
    
    
    
def bond_count(SMILES,pattern):
    name = SMILES
    mol = Chem.MolFromSmiles(name)
    rd_count = []
    for sub in pattern:
        # print(sub)        
        sub = sub.replace('R','')
        rd_count.append(len(mol.GetSubstructMatches(Chem.MolFromSmarts(sub))))
    # rd_count.append(H_count(SMILES))    
    return rd_count   
 

        

# =============================================================================
# parameters/fragments loading
# =============================================================================
pattern_transform('QM8_patterns.txt')
# contribution,intercept = pattern_import('pattern_contributions_v2.txt')


contribution,intercept = pattern_import('pattern_contributions_OLS.txt')



    

files = ['uncertainty_information0.txt',
        'uncertainty_information1.txt',
        'uncertainty_information2.txt',
        'uncertainty_information3.txt',
        'uncertainty_information4.txt',
        'uncertainty_information5.txt']
UncertFile = ''
for file in files:
  UncertFile += open(file,'r').read()
sigma,x_ld = uncertainty_loading(UncertFile)

# =============================================================================
# uncertainty
# =============================================================================
fidx = list(np.nonzero(contribution)[0])


def file_transform(file,atoms,mode = 'regular'):
    # print(mode)
    bonds = loading_bonds(mode)
    pattern_name = atoms + bonds
    print(pattern_name)
    flag = 0
    f = open(file,'r').readlines()
    # f2 = open(file.replace('.txt','_')+'counts_'+mode+'.txt','w')
    results = ''
    strs = '{:80}{:20.2f}{:20.2f}'
    for i,line in enumerate(f):
        # print(line)
        lis = line.split()
        smiles = lis[0]
        # print(smiles)
        if len(lis) == 1:
            energy = predict_energy(smiles,contribution,intercept)
            lower = lower_matrix(smiles,fidx)
            uncertainty = uncertainty_prediction(sigma,x_ld,lower)
        else:
            energy = float(lis[1])
            uncertainty = 0
            flag = 1
        # print(energy)
        atom = atom_count(smiles,atoms)
        bond = bond_count(smiles,bonds)
        if i == 0:
            count_length = '{:15}'*(len(atom)+len(bond))
            strs += count_length
            if flag == 0:
                results += strs.replace('.2f','').replace(':1',':>1').format('SMILES','atomization_energy','uncertainty',*pattern_name)+'\n'
                # f2.write(strs.replace('.2f','').replace(':1',':>1').format('SMILES','atomization_energy','uncertainty',*pattern_name)+'\n')
            if flag == 1:
                # f2.write(strs.replace('.2f','').replace(':1',':>1').format('SMILES','formation_energy','uncertainty',*pattern_name)+'\n') 
                results += strs.replace('.2f','').replace(':1',':>1').format('SMILES','atomization_energy','uncertainty',*pattern_name)+'\n'               
        # f2.write(strs.format(smiles,energy,uncertainty,*atom,*bond)+'\n')
    # f2.close()    
    return strs.format(smiles,energy,uncertainty,*atom,*bond)+'\n'

atoms = ['C','O','H']

# mode either be 1) regular, 2) hybridization or 3) pathway fps
file_transform('test_smiles.txt',atoms, mode = 'regular')





