import numpy as np
import cplex
import random
import pandas as pd
import math
import sys
from sphinx.addnodes import index

filenName = 'southeast-asia.xlsx'
threshHold = 0
numberCandidates = 3

def readData(fileName,threshHold,numberCandidates):
    dataset = pd.read_excel(fileName)
    '''
    R[i][j] is skill-depth of candidates i-th in skill j-th
    '''
    R = dataset.iloc[:500,1:].values
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
tau = 0.001
totalCandidates = len(R)
# z = [elem*0.4 for elem in E]
numberSkill = len(R[0])
z = np.zeros(numberSkill)
# c = [random.randrange(1, 50, 1) for i in range(totalCandidates)]
c = np.zeros(totalCandidates)
C = 1000000000

class MDSB:
    def __init__(self):
        self.model = cplex.Cplex()
    def buildModel(self):
        self.model.objective.set_sense(self.model.objective.sense.minimize)
        variables = []
        for i in range(totalCandidates):
            variables.append('x' + str(i))
        '''
        Add variables and set upper bounds & lower bounds
        lower bound of all variables is 0 
        upper bound of all variables is 1
        '''
        ub = [1 for i in range(totalCandidates)]
        lb = [0 for i in range(totalCandidates)]
        self.model.variables.add(names=variables,lb=lb,ub=ub)
        '''
        Set up variables types is binary
        ['C', 'I', 'B', 'S', 'N'] are integer, binary, semi_continuous, semi_integer
        '''
        self.model.variables.set_types([(elem,self.model.variables.type.binary) for elem in variables])
        '''
        Set up model function, after decomposition
        Quadratics coefficient for candidate i-th will be: 
            R[i,1] ^2  + R[i,2] ^2 + ... + R[i,m] ^ 2 -tau
        Linear coefficient for candidate i-th will be:
            -(2*(E[1] * R[i,1] + E[2] * R[i,2] + .... + E[m] * R[i,m]) + tau 
        '''
        ## Calculating quadratic coefficients
        # quad_coefficients = np.sum(np.multiply(R, R), axis=1)
        quad_variables = []
        # for i in range(totalCandidates):
        #     quad_variables.append((i, i, quad_coefficients[i]-tau))

        for skill_Index in range(numberSkill):
            for canIndex1 in range(numberCandidates):
                for canIndex2 in range(numberCandidates):
                    quad_variables.append((canIndex1,canIndex2,R[canIndex1][skill_Index] * R[canIndex2][skill_Index]))
        self.model.objective.set_quadratic_coefficients(quad_variables)
        ## Linear Coefficients
        R_copy = np.copy(R)
        for i in range(len(R_copy[0])):
            R_copy[:, i] = [2 * element * E[i] for element in R_copy[:, i]]
        linear_coefficients = np.sum(R_copy, axis=1)

        linear_coefficients = [(linear_coefficients[index] - tau) * -1 for index in
                               range(len(linear_coefficients))]
        linear_variables = []
        for i in range(totalCandidates):
            linear_variables.append((i, linear_coefficients[i]))
        self.model.objective.set_linear(linear_variables)
        '''
        Setting ups constraints 
        1. SUM of all variables is number of candidates
        2. Minimum Skill Score 
        3. X[i] is either 0 or 1 .
        4. Constraint sum(c[i]*x[i]) <= C
        '''
        # Constraint #1
        linear_constraints = [cplex.SparsePair(ind=[index for index in range(totalCandidates)],val=[1 for i in range(totalCandidates)])]
        rhs_linear_constraints = [numberCandidates]
        senses = ['E']
        # Constraint #2
        '''
        senses = ["E", "L", "G", "R"]
        Each entry must be one of 'G', 'L', 'E', and 'R', indicating greater-than or equal, less-than or equal, equality, and ranged constraints,respectively.
        '''
        for i in range(numberSkill):
            row = [R[index, i] for index in range(totalCandidates)]
            linear_constraints.append(cplex.SparsePair(ind=variables.copy(),val=row))
            rhs_linear_constraints.append(z[i])
            senses.append('G')
        self.model.linear_constraints.add(lin_expr=linear_constraints, rhs=rhs_linear_constraints, senses=senses)
        # Constraint #4
        self.model.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=variables.copy(), val=c)],
            senses=['L'],
            rhs=[C])

class PMDSB:
    def __init__(self):
        self.model = cplex.Cplex()
    ##
    def buildModel(self , previous_step):
        self.model.objective.set_sense(self.model.objective.sense.minimize)
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
        ub = [1 for i in range(numberCandidates)]
        lb = [0 for i in range(numberCandidates)]
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
        self.model.objective.set_quadratic_coefficients(quad_variables)
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
        4. Constraint sum(c[i]*x[i]) <= C
        '''
        # Constraint #1
        linear_constraints = [
            cplex.SparsePair(ind=[index for index in range(totalCandidates)], val=[1 for i in range(totalCandidates)])]
        rhs_linear_constraints = [numberCandidates]
        senses = ['E']
        # Constraint #2
        '''
        senses = ["E", "L", "G", "R"]
        Each entry must be one of 'G', 'L', 'E', and 'R', indicating greater-than or equal, less-than or equal, equality, and ranged constraints,respectively.
        '''
        for i in range(numberSkill):
            row = [R[index, i] for index in range(totalCandidates)]
            linear_constraints.append(cplex.SparsePair(ind=variables.copy(), val=row))
            rhs_linear_constraints.append(z[i])
            senses.append('G')
        self.model.linear_constraints.add(lin_expr=linear_constraints, rhs=rhs_linear_constraints, senses=senses)
        # constraints #3 x(x-1) = 0
        for i in range(totalCandidates):
            linear_constraints_of_quadratics = cplex.SparsePair(ind=[variables[i]], val=[1])
            quad_constraints = cplex.SparseTriple(ind1=[variables[i]], ind2=[variables[i]], val=[-1])
            quad_sense = 'E'
            quad_contraints_rhs = 0
            self.model.quadratic_constraints.add(lin_expr=linear_constraints_of_quadratics,
                                                 quad_expr=quad_constraints,
                                                 sense=quad_sense,
                                                 rhs=quad_contraints_rhs)
        # Constraint #4
        self.model.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=variables.copy(), val=c)],
            senses=['L'],
            rhs=[C])

def compute_Objective(previousState,X):
    sum = 0
    for j in range(numberSkill):
        tmp = 0
        for i in range(numberCandidates):
            tmp+=R[i][j] * X[i]
        sum += (E[j] - tmp) * (E[j] - tmp)
    ##
    for i in range(numberCandidates):
        sum -= tau * (2*previousState[i]-1) * X[i]
    return sum

def solving(option):
    # option 1 is MDSB
    # option 2 is PMDSB
    if option == 1:
        modelMDSB = MDSB()
        modelMDSB.buildModel()
        modelMDSB.model.solve()
        if (modelMDSB.model.solution.status.infeasible):
            print("Infeasible")
            return
        print("Solution for problems using MDSB model is: ")
        print(modelMDSB.model.solution.get_values())

    else:
        ## random pick
        indices = [i for i in range(totalCandidates)]
        global X_0
        global isSolvedtotalCandidates
        global previousValue
        '''
        Assign initial variables
        '''
        isSolved = True
        X_0 = np.zeros(totalCandidates)
        previousValue = sys.float_info.max
        '''
        '''
        for i in range(numberCandidates):
            rand_index = random.choice(indices)
            X_0[rand_index] = 1
            indices.remove(rand_index)
        ##
        epsilon = 0.00001
        while True:
            modelPMDSB = PMDSB()
            modelPMDSB.buildModel(X_0)
            modelPMDSB.model.solve()
            if (modelPMDSB.model.solution.status.infeasible):
                print("Infeasible")
                isSolved = False
                break
            X_1 = modelPMDSB.model.solution.get_values()
            currentValue = compute_Objective(X_0,X_1)
            if (math.fabs(previousValue - currentValue) < epsilon):
                X_0 = X_1
                break
            X_0 = X_1
            previousValue = currentValue

        if (isSolved):
            print("Solution for problems using PMDSB is: ")
            print(X_0)

solving(1)

# A = [[1,2],[3,4]]
#
# print(np.sum(A,axis=0))


