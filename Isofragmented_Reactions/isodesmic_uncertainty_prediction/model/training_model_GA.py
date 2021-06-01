# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 15:16:17 2020

@author: bol216
"""

from rdkit import Chem
from numpy.linalg import inv
import numpy as np
import numpy.linalg as lin

def length(mol):
    # return the number of heavy atoms in the molecule
    size = 0
    for letter in mol:
      if letter == 'C' or letter == 'c':
        size += 1
      if letter == 'F' or letter == 'f':
        size += 1
      if letter == 'O' or letter == 'o':
        size += 1
      if letter == 'N' or letter == 'n':
        size += 1
    return size  


def common_member(a, b):  
    # return common elements in two vectors
    a_set = set(a)
    b_set = set(b)
    if len(a_set.intersection(b_set)) > 0:
        return len((a_set.intersection(b_set)))
    else:
        return 0

def ringselect_rdkit(e,name):
 # return the sizes of two ring fused together as additional feature information. 
 # e.g. the correction term for C1CC2C1CCC2 is fused4+fused5.  
 numbers = []
 for i in range(len(name)):
        letter = name[i]
        if letter.isdigit() == True:
            numbers.append(letter)
 fused_comb = []
 for i in range(3,9):
    for k in range(3,i+1):
        fused_comb.append(str(i)+str(k))  
        
 fused_type = [0]*len(fused_comb)
 if len(numbers) <= 2:
  e += fused_type
 else:
   fused = 0
   m = Chem.MolFromSmiles(name)
   ssr = Chem.GetSymmSSSR(m)
   ring = []
   for ele in ssr:  
     ring.append(list(ele))     
   for i in range(1,len(ring)):
    for k in range(i):
        leng = []
        if common_member(ring[i],ring[k]) == 2 and len(ring[i]) < 9 and len(ring[k]) < 9:
            fused += 1
            leng.append(len(ring[i]))
            leng.append(len(ring[k]))
            comb = str(max(leng))+str(min(leng))
            fused_type[fused_comb.index(comb)] += 1   
   e += fused_type
 return e
       
def ringcorrection(name,AAsize,ARsize):
    # correction term for molecules with ring fragments. the format is the aromaticity + the size.
    # e.g. the term for 'c1ccccc1' is aromatic + size 6.
    Ring = np.zeros(len(ARsize))
    Aring = np.zeros(len(AAsize))
    m = Chem.MolFromSmiles(name)
    ssr = Chem.GetSymmSSSR(m)    
    for ele in ssr:
        size = len(ele)
        if m.GetAtomWithIdx(ele[0]).GetIsAromatic() == True:
            if size < 7:
               Aring[AAsize.index(size)] += 1
        else:
            if size < 9:
                Ring[ARsize.index(size)] += 1
    return list(Ring) + list(Aring)


def H_count(SMILES):
#    strs = ['C','O','N','S']
    m = Chem.MolFromSmiles(SMILES)
    num1 = m.GetNumAtoms(True)
    mh = Chem.AddHs(m)
    num2 = mh.GetNumAtoms(True)
    H = num2 - num1
#    S = SMILES.upper()
#    C,O,N,F = [S.count(ele) for ele in strs]
    return H
    

def rdkit_count(pattern,SMILES):
    name = SMILES
    mol = Chem.MolFromSmiles(name)
#    if mol == None:
#        print('error')
#        return 1
    rd_count = []
    for sub in pattern:
        sub = sub.replace('R','')
#        print(sub)
    #    print(len(mol.GetSubstructMatches(Chem.MolFromSmarts(sub))))
        rd_count.append(len(mol.GetSubstructMatches(Chem.MolFromSmarts(sub))))
    rd_count.append(H_count(SMILES))    
    return rd_count    


def pattern_transform(file):
    global AAsize
    global ARsize
    global pattern
    global fname
    fname = []
    pattern = []
    pat = open(file,'r').readlines()
    for line in pat:
       lis = line.split()
       if lis[0] == 'f':
           continue
       fname.append(lis[0])
       pattern.append(lis[0])
    fname += ['H']       
    fused_comb = []
    for i in range(3,9):
      for k in range(3,i+1):
          fused_comb.append(str(i)+str(k))   
    fname += fused_comb
    idx = [i for i in range(len(fname)) if fname[i].find('R1') != -1]
    ARsize = []
    AAsize = []
#    print('ARsize',ARsize)
    ARsize.append(8)         
#    print('AAsize',AAsize)      
    AAsize.append(4)    
    for ele in [fname[i] for i in idx]:
      if ele.find('c') == -1:
          if length(ele) not in ARsize:
              ARsize.append(length(ele))
      else:
          if length(ele) not in AAsize:
              AAsize.append(length(ele))
    for i in ARsize:
          fname += ['ARsize' + str(i)]
    for i in AAsize:
          fname += ['AAsize' + str(i)]          
    return 1
#
def sanity_check(name):
  mol = Chem.MolFromSmiles(name)
  if mol == None:
      print('error in '+name, ' not a rdkit structure')
      return 1
  else:
      return 0
    
def readmol_rdkit(name):
  # include the molecules fingerprints, properties as a dictionary.
  frag = rdkit_count(pattern,name)
  frag = ringselect_rdkit(frag,name)  
  frag = np.hstack((frag,ringcorrection(name,AAsize,ARsize)))
  return frag

#

def pattern_import(file):
    f = open(file,'r').readlines()
    fname = []
    a = []
    for line in f:
        lis = line.split()
        if line.startswith('inter'):
            intercept = float(lis[1])
        else:    
            fname.append(lis[0])
            a.append(float(lis[1]))
    return a,intercept

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


def atomization_heat(mol,ene):
    atoms = ['C','O','H']
    H_ene = 217.998
    O_ene = 249.229
    C_ene = 716.884
    C,O,H = atom_count(mol,atoms)
    return ene + np.array([C,O,H])@np.array([C_ene,O_ene,H_ene])/4.184


def predict_energy(name,contribution,intercept):
    # name = 'C1CCCCC1' # change molecules by changing the SMILES for the name    
    # the following code uses rdkit module to count each fragments and group additivity theory to sum them up to the energies.
    # the contributions for each fragement is stored in the contribution file.
    flag = sanity_check(name)
    if flag == 1:
        print('error')
    else:    
        # pattern_transform('QM8_patterns.txt')  
        # fname,contribution,intercept = pattern_import('pattern_contributions_v2.txt')
        fragments = readmol_rdkit(name)
        # print(fragments)
        energy = np.dot(fragments,np.transpose(contribution))+intercept
        # print('the predicted atomization for '+name+' is ' +str(energy) + ' kcal/mol')
    # return atomization_heat(name,energy) 
    return energy   

def select_matrix(vec,idx):
    Nvec = []
    for i in range(len(vec)):
        row = [vec[i][j] for j in idx]
        Nvec.append(row)
    return Nvec

def lower_matrix(name,fidx):
    flag = sanity_check(name)
    if flag == 1:
        print('error')    
    else:
        fragments = readmol_rdkit(name)   
    return np.array(fragments[fidx])
    
def uncertainty_loading(files):   
    # f = open(file,'r').readlines()
    f = file.splitlines()
    matrix = []
    for line in f:
        lis = line.split()
        if line.startswith('sigma'):
            sigma = float(lis[1])
        else:    
            row = []
            for ele in lis:
                row.append(float(ele))
            matrix.append(row)    
    return sigma,np.array(matrix)

def prediction_variance(V_t,x):
    return np.array(x).dot(V_t).dot(np.array(x).T)  

def uncertainty_prediction(sigma,x_ld,x_test):
    dispersion_matrix = lin.inv(np.array(x_ld).T.dot(x_ld))
    return sigma*np.sqrt(prediction_variance(dispersion_matrix,x_test))


    

