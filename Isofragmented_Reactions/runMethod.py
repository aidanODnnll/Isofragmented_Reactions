import pandas as pd
import scipy.stats
import numpy as np
import os
from pyomo.environ import *
from rdkit import Chem

# os.environ['NEOS_EMAIL'] = 'user@email' # Replace with your own email, if you want to run CPLEX through NEOS 


def justSMILES(filename):
    lines = open(filename, 'r').read().splitlines()
    SMILES = ''
    for l in lines:
        l = l.split()
        SMILES += l[0] + '\n'
    return SMILES


def keyMaker(filename):
    f = np.loadtxt(filename, dtype=object, comments = None)
    key = {}
    for l in f:
        key.update({l[0]:float(l[1])})
    return key

def CreateData(string):
    l = string.splitlines()
    z = 0 
    res = []
    #This for loop puts the data into a list to be used later 
    for x in l:
        z += 1
        # res += [delspace(x)]
        res += [x.split()]
    MolToFactors = MakeDataDict(res[0],res[1:])
    return MolToFactors

def MakeDataDict(features,data):
    res = {}
    newf = features[1:]
    mols = []
    z = 0
    #Goes through each index of data
    for l in data:
        mol = l[0] 
        mols += [mol]
        point = 1 
        # assigns the value to its corresponding feature for each molecule      
        for i in newf:
            res.update({(mol,i):float(l[point])})
            point += 1
    return res

def StoichRule1(m, i):
    return -1*m.n*m.y[i] + m.z[i] <= m.s[i]

def StoichRule2(m, i):
    return   -1*m.y[i] + 1*m.n*m.z[i] >= m.s[i]

def OneSRule(m, i):
    return m.y[i] + m.z[i] <= 1

def slack(m,k):
    sums = sum(m.a[i,k]*m.s[i] for i in m.mol)
    return sums <= m.c[k]

def slack2(m,k):
    sums = sum(m.a[i,k]*m.s[i] for i in m.mol)
    return sums >= -1*m.c[k]

def makeEta(m,i):
    return -1*m.eta[i] <= m.s[i]

def makeEta2(m,i):
    return m.eta[i] >= m.s[i]

def MakeIso(m):
    sums = sum(m.c[i]/m.w[i] for i in m.factors if  m.w[i] != 0 and m.w[i] <= 100  )
    return sums <= m.Lambda

def findRxns(UnknownMoleculesToFeatures, Weights, NumberOfRxns, maxStoich, KnownMolecules, KnownMoleculesToFeatures, features, mol):   
    d = KnownMoleculesToFeatures
    newfactors = {}
    #Adds Target molecule to the dictionary of features and list of molecules
    for x in features:
        z = UnknownMoleculesToFeatures[(mol,x)]
        d.update({('target', x):z})
    KnownMolecules += ['target']
    targetMol = 'target'
    m = ConcreteModel()
    # sets
    m.mol = Set(initialize = KnownMolecules) # the molecules
    m.factors = Set(initialize = features)# factors
    #parameters
    m.a = Param(m.mol, m.factors, initialize = d)# gives values of factors for molecules
    m.w = Param(m.factors, initialize = Weights) # gives weights of factors
    m.n = Param(initialize = maxStoich, mutable = True) # Max stoich coeff
    m.Lambda = Param(initialize = 0, mutable = True)

    #variables 
    m.y = Var(m.mol, initialize = 0, within = Binary)# Whether it is a reactant
    m.z = Var(m.mol, initialize = 0, within = Binary)# whether it is a product
    m.s = Var(m.mol, initialize = 0, within = Reals)# stoich coeff
    m.c = Var(m.factors, initialize = 0, within =NonNegativeReals) # This is the amount of extra bonds needed
    m.eta = Var(m.mol, initialize = 0, within =NonNegativeReals) # Positive value of S
    # con
    m.chosen = Constraint( expr = m.z[targetMol] >= 1) # Ensures target molecule is in all reactions 
    m.StoichLim1 = Constraint(m.mol, rule = StoichRule1) # Makes sure we don't go over the stoich limit
    m.StoichLim2 = Constraint(m.mol, rule = StoichRule2)
    m.OneSide = Constraint(m.mol, rule = OneSRule) # Make sure the molecule only apears on one side
    m.slack = Constraint(m.factors, rule = slack)# Determines m.c
    m.slack2 = Constraint(m.factors, rule = slack2) 
    m.isoD = Constraint(rule = MakeIso)#Makes the reaction isofragmented
    m.posStoich =  Constraint(m.mol, rule = makeEta) # Determines m.eta
    m.posStoich2 =  Constraint(m.mol, rule = makeEta2)
    m.constL = ConstraintList()
    # Ensuring that elements must be conserved 
    expr = 0
    for j in m.w:
        if value(m.w[j]) == 0:
            expr += m.c[j]
    m.constL.add( expr == 0 )
    
    #objective
    m.Objective = Objective(expr=sum(m.eta[i]*m.a[i,'uncertainty'] for i in m.mol), sense = minimize)
    
    #Solver options 
    # Uncomment to use NEOS
    solver = SolverManagerFactory('neos')

    # Uncomment to use Gurobi
    # solver = SolverFactory('gurobi')
    # solver.options['MIPGapAbs'] = (.1)
    # solver.options['TimeLimit'] = 250
    
    # Uncomment to use CPLEX
    # solver = SolverFactory("cplex")
    # solver.options['mip_tolerances_absmipgap'] = .1
    # solver.options['timelimit'] = 250
    
    sol = ''
    t = 1
    # iterate through reactions 
    for i in range(NumberOfRxns):
        # Uncomment if using any solver but NEOS
        results = solver.solve(m,load_solutions=False)

        # Uncomment if using NEOS
        # results = solver.solve(m, opt = 'cplex', load_solutions=False)

        # Check solver status 
        if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition == TerminationCondition.optimal):
            m.solutions.load_from(results)

            ## Options to see the results as displayed by the solver
            # print(results) 
            # print(m.display())
            
            # Turn the solution into a more understandable format
            res = ''
            resC = 0
            res2C = 0
            res2 = ''
            for v in m.component_objects(Var, active=True):
                if str(v) == 's':
                    for i in v:
                        if abs(value(v[i])) >= .01:
                            if value(v[i]) <= 0:
                                if resC >= 1:
                                    res += ' + '
                                res += str(round(abs(value(v[i])),3)) + '*' + str(i) 
                                resC += 1
                            if value(v[i]) >= 0:
                                if res2C >= 1:
                                    res2 += ' + '
                                res2 += str(round(abs(value(v[i])),3)) + '*' + str(i) 
                                res2C += 1
            sol += 'Reaction ' + str(t) + ':' + '      ' +  res + ' <--> ' + res2 + ' \n'
       
        # Possible failure cases 
        elif (results.solver.termination_condition == TerminationCondition.infeasible):
            break
        elif results.solver.status == SolverStatus.aborted:
            break
        # print(t) # Allows user to keep track of what reaction number the model is on (not necessary for use)
        t += 1
        expr = 0
        ### Adjusted Integer cuts, uncomment else statements for traditional integer cuts(not recommended)
        for j in m.y:
            if value(m.y[j]) == 1:
                expr += (1-m.y[j])
            #else:
                #expr += m.y[j]
        for r in m.z:
            if value(m.z[r]) == 1:
                expr += (1-m.z[r])
            #else:
                #expr += m.z[r]
        m.constL.add( expr >= 1 )
    return(sol,t-1)

def determineEnthalpy(rxns, n, KnownEnthalpy, KnownMoleculesToFeatures, UnknownMoleculesToFeatures):
    lines = rxns.splitlines()
    rxnCount = 0
    listR = []
    res = [] 
    dHList = []
    foundRxn = False
    for l in lines: 
        if l[:6] == 'target':
            foundRxn = True
            # target = delspace(l[8:])[0]
            # listR += delspace(l[8:])
            target = l[8:].split()[0]
            listR += l[8:].split()
            pass
        elif foundRxn == True:
            # listR += [delspace(l[13:])]
            listR += [l[13:].split()]
            rxnCount += 1
            if rxnCount == n:
                Estimate, dHList = Calculator(KnownEnthalpy,listR, n, \
                                              KnownMoleculesToFeatures, UnknownMoleculesToFeatures)
                foundRxn = False 
    return Estimate, dHList

def calcHrxn(listMol, stoich,chosenStoich, target, d, UnknownMoleculesToFeatures):
    Hrxn = sum(stoich[m]*d[(m,'atomization_energy')]  for m in listMol)\
     + chosenStoich*UnknownMoleculesToFeatures[(target,'atomization_energy')]
    return Hrxn  

def Calculator(KnownEnthalpy, listR, numRxn, KnownMoleculesToFeatures, UnknownMoleculesToFeatures):
    d = KnownMoleculesToFeatures
    target = listR[0]
    Estimate = []
    HrxnSum = 0
    dHrxnList = []
    for l in range(numRxn+1):
        #First make list of all molecules involved and a dictionary linking them to their coeff
        listMol = []
        stoich = {}
        reactant = True
        targetStoich = 0
        if l == 0:
            pass
        else:
            stillLooking = True
            for i in listR[l]:
                coeff = ''
                coeffStatus = True
                currentMol = ''  
                if i == '<-->':
                    reactant = False
                else: 
                    for k in range(len(i)):
                        if stillLooking == False:
                            break
                        if (i[k].isnumeric() or i[k] == '.') and coeffStatus:
                            coeff += i[k] 
                        elif i[k] == '*':
                            coeffStatus = False
                        elif i[k] == '+':
                            break
                        elif i[k:] == 'target':
                            if reactant:
                                targetStoich = float(coeff)*-1
                            else:
                                targetStoich = float(coeff)
                            stillLooking = False
                            break
                        else:
                            currentMol += i[k]    
                    if currentMol != '':
                        listMol += [currentMol]
                        if reactant:
                            stoich.update({currentMol:float(coeff)*-1})
                        else:
                            stoich.update({currentMol:float(coeff)})
            if  targetStoich != 0:
                Hrxn = calcHrxn(listMol, stoich,targetStoich, target, d, UnknownMoleculesToFeatures)
                Estimate += [(Hrxn - sum(KnownEnthalpy[m]*stoich[m] for m in listMol)/targetStoich)]
                dHrxnsqr = sum(abs(stoich[m])*d[(m,'uncertainty')]**2 for m in listMol) \
                + abs(targetStoich)*UnknownMoleculesToFeatures[(target,'uncertainty')]**2
                dHrxn = dHrxnsqr**.5/abs(targetStoich)
                dHrxnList += [dHrxn]
    return Estimate, dHrxnList
Estimate, dHList = determineEnthalpy(results, n, EnthalpyOfFormation,\
                                     KnownMoleculesToFeatures, UnknownMoleculesToFeatures)
def findAverages(Estimate, dHlist):
    MultH = 0
    n = len(Estimate)
    for l in range(n):
        MultH += np.log10(float(Estimate[l]) + 10000)
    logM = MultH/n
    standardDev = np.std(Estimate)
    mean = 10**logM - 10000
    correctedMult = 0
    cN = 0 
    if len(Estimate) == 1:
        cMean = listH[0]
    else:
        for x in Estimate:
            if mean - 2*standardDev <= x <= mean + 2*standardDev:
                correctedMult += np.log10(x + 10000) 
                cN += 1
            elif standardDev == 0:
                correctedMult += np.log10(x + 10000) 
                cN += 1
        clogM = correctedMult/cN
        cMean = 10**clogM - 10000
        meanE = np.mean(dHlist)
        stdevE = np.std(dHlist)
        df = n-1
        lower, upper = scipy.stats.t.interval(.90, df, meanE, stdevE)
    return(cMean, upper)



def run(UnknownMoleculesTxt, KnownMoleculesTxt = [], UseNIST = True):
    NumberOfRxns = 40
    StoichLimit = 5
    f = open(KnownMoleculesTxt)
    if len(f.readlines()) >= 1:
        EnthalpyOfFormation = keyMaker(KnownMoleculesTxt)
        KnownSMILES = justSMILES(KnownMoleculesTxt)
        KnownMolecules = KnownSMILES.splitlines()
        String_KMtF, values = file_transform(KnownSMILES,atoms, mode = 'pathway fps')
        KnownMoleculesToFeatures = CreateData(String_KMtF)
    else: 
        EnthalpyOfFormation = {}
        KnownSMILES = []
        KnownMolecules = []
        KnownMoleculesToFeatures = {}
    f.close()
    f = open(UnknownMoleculesTxt)
    if len(f.readlines()) < 1:
        return 'Error: No Molecules with Unknown Enthalpy Provided'
    UnknownMolecules = open(UnknownMoleculesTxt, 'r').readlines()
    String_UMtF, values = file_transform(open(UnknownMoleculesTxt).read(),atoms, mode = 'pathway fps') 
    UnknownMoleculesToFeatures = CreateData(String_UMtF)
    features = values[2:]
    GATestimates = values[:2]
    if UseNIST: 
        NISTKnownFile = 'NIST_ene.txt'
        NISTEnthalpy = keyMaker(NISTKnownFile)
        KnownMolecules += justSMILES(NISTKnownFile).splitlines()
        KnownMolecules = list(set(KnownMolecules))
        f = open(NISTKnownFile)
        NISTString, values = file_transform(f.read(),atoms, mode = 'pathway fps')
        f.close()
        print(CreateData(NISTString))
        KnownMoleculesToFeatures.update(CreateData(NISTString))
        EnthalpyOfFormation.update(NISTEnthalpy)
    Weights = {}
    for f in features: 
        if len(f) == 1: 
            Weights.update({f:0})
        else:
            Weights.update({f:1})
    for e in GATEstimates:
        Weights.update({e:1000})
    solutions = 'Molecule \t Estimated Enthalpy of Formation (kcal/mol) \t Uncertainty \n'
    for mol in UnknownMolecules:
        if mol in KnownMolecules:
            solutions += str(mol) +'\t' + str(EnthalpyOfFormation[mol]) + '\t Enthalpy already known \n'
        else:
            results = 'chosen = '+ str(mol) + '\n'
            rxns,n = findRxns(UnknownMoleculesToFeatures, Weights, NumberOfRxns, StoichLimit, KnownMolecules, KnownMoleculesToFeatures, values, mol)
            results += rxns
            Estimate, dHList = determineEnthalpy(results, n, EnthalpyOfFormation, KnownMoleculesToFeatures, UnknownMoleculesToFeatures)
            Mean, Error = findAverages(Estimate, dHList)
            solutions += str(mol) +'\t' + str(Mean) + '\t' + str(Error) + '\n'
    return solutions




