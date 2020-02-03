import pandas as pd
import numpy as np
import cplex
import random
from sphinx.addnodes import index

filenName = 'southeast-asia.xlsx'
threshHold = 0
numberCandidates = 20

def readData(fileName,threshHold,numberCandidates):
    dataset = pd.read_excel(fileName)
    '''
    R[i][j] is skill-depth of candidates i-th in skill j-th
    '''
    R = dataset.iloc[:,1:].values
    '''
    Sorting based on skill-depth by descending order
    '''
    R2 = -np.sort(-np.array(R),axis=0)
    '''
    E is sum of h(number of selected candidates) largest skill-depth
    '''
    E = []
    '''
    skillScore is sum of all skill-depth of a candidates
    '''
    skillScore = np.sum(R,axis=0)
    ''''''
    for i in range(len(R2[0])):
        sum_skill_largest = 0
        for j in range(numberCandidates):
            sum_skill_largest += R2[j,i]
        E.append(sum_skill_largest)
        ##
    return (E,R,skillScore)

E, R, skillScore = readData(fileName= filenName, threshHold= threshHold, numberCandidates= numberCandidates)
tau = 0.000001
z = np.zeros(len(R[0]))
totalCandidates = len(R)
numberSkill = len(R[0])

class MDSB:
    def __init__(self):
        self.model = cplex.Cplex()
    def buildModel(self):
        self.model.set_sense(self.model.objective.sense.minimize)
        return 0
    def solving(self):
        return 0

class PMDSB:
    def __init__(self):
        self.model = cplex.Cplex()
    ##
    def buildModel(self , previous_step):
        self.model.set_sense(self.model.objective.sense.minimize)
        '''
        Number of variables equals to number of candidates
        '''
        variables = []
        for i in range(totalCandidates):
            variables.append('x'+str(i))
        '''
        Add variables and set upper bounds & lower bounds
        lower bound of all variables is 0 
        upper bound of all variables is 1
        '''
        ub = []
        lb = []
        for i in range(totalCandidates):
            ub.append((variables[i],1))
            lb.append((variables[i],0))
        self.model.variables.add(names=variables,lb=lb,ub=ub)
        '''
        Set up model function, after decomposition
        Quadratics coefficient for candidate i-th will be: 
            R[i,1] ^2 + R[i,2] ^2 + ... + R[i,m] ^ 2
        Linear coefficient for candidate i-th will be:
            -(2*(E[1] * R[i,1] + E[2] * R[i,2] + .... + E[m] * R[i,m]) + tau * ( 2x previous_step[i] - 1 ))
        '''
        ## Calculating quadratic coefficients
        quad_coefficients = np.sum(np.multiply(R,R),axis=1)
        quad_variables = []
        for i in range(totalCandidates):
            quad_variables.append((i,i,quad_coefficients[i]))
        ## Linear Coefficients
        R_copy = np.copy(R)
        for i in range(len(R_copy[0])):
            R_copy[:,i] = [2*element * E[i] for element in R_copy[:,i]]
        linear_coefficients = np.sum(R_copy,axis=1)
        linear_coefficients = [(linear_coefficients[index] + tau * (2*previous_step[index] - 1)) * -1 for index in range(len(linear_coefficients))]
        linear_variables = []
        for i in range(totalCandidates):
            linear_variables.append(i, linear_coefficients[i])
        self.model.objective.set_linear(linear_variables)
        '''
        Setting ups constraints 
        1. SUM of all variables is number of candidates
        2. Minimum Skill Score 
        3. X[i] is either 0 or 1 .
        '''
        # Constraint #1
        linear_constraints = [[index for index in range(totalCandidates)]]
        rhs_linear_constraints = [totalCandidates]
        senses = ['E']
        # Constraint #2
        '''
        senses = ["E", "L", "G", "R"]
        Each entry must be one of 'G', 'L', 'E', and 'R', indicating greater-than or equal, less-than or equal, equality, and ranged constraints,respectively.
        '''
        for i in range(numberSkill):
            row = [R[index, i] for index in range(totalCandidates)]
            linear_constraints.append(row)
            rhs_linear_constraints.append(z[i])
            senses.append('G')
        self.model.linear_constraints.add(lin_expr=linear_constraints,rhs=rhs_linear_constraints,senses=senses)
        # constraints #3 x(x-1) = 0
        linear_constraints_of_quadratics = cplex.SparsePair(ind=variables,val=np.ones(totalCandidates))
        quad_constraints = cplex.SparseTriple(ind1=variables,ind2=variables,val=np.ones(totalCandidates))
        quad_sense = 'E'
        quad_contraints_rhs = np.zeros(totalCandidates)
        self.model.quadratic_constraints.add(lin_expr=linear_constraints_of_quadratics,
                                             quad_expr=quad_constraints,
                                             sense=quad_sense,
                                             rhs=quad_contraints_rhs)

def solving(option):
    # option 1 is MDSB
    # option 2 is PMDSB
    if option == 1 :
        print("hello")
    else:
        ## random pick
        indices = [i for i in range(totalCandidates)]
        X_0 = np.zeros(totalCandidates)
        for i in range(numberCandidates):
            rand_index = random.choice(indices)
            X_0[rand_index] = 1
            indices.remove(rand_index)
        ##
        epsilon = 0.0000001
        while True:
            print(0)
