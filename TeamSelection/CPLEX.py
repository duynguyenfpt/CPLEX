from builtins import print

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

def readData(fileName,threshHold,numberCandidates,numberSkill):
    dataset = pd.read_excel(fileName)
    '''
    R[i][j] is skill-depth of candidates i-th in skill j-th
    '''
    R = dataset.iloc[:500,1:numberSkill].values
    nickNames = dataset.iloc[:500,0].values
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
            sum_skill_largest += R2[j][i]
        E.append(sum_skill_largest)
        ##
    tau = 0.001
    totalCandidates = len(R)
    z = [elem * 0.3 for elem in E]
    numberSkill = len(R[0])
    # z = np.zeros(numberSkill)
    # c = [random.randrange(1, 50, 1) for i in range(totalCandidates)]
    c = np.zeros(totalCandidates)
    C = cplex.infinity
    return (E,R,skillScore,nickNames,tau,z,c,C,numberSkill,totalCandidates)

class MDSB:
    def __init__(self):
        self.model = cplex.Cplex()
    def buildModel(self,E, R,tau,z,c,C,numberSkill,totalCandidates):
        self.model.objective.set_sense(self.model.objective.sense.minimize)
        variables = []
        for i in range(totalCandidates):
            variables.append('x' + str(i))
        '''
        Add variables and set upper bounds & lower bounds
        '''
        self.model.variables.add(names=variables)
        '''
        Set up variables types is binary
        ['C', 'I', 'B', 'S', 'N'] are integer, binary, semi_continuous, semi_integer
        '''
        self.model.variables.set_types([(elem,self.model.variables.type.binary) for elem in variables])
        '''
        Set up model function, after decomposition
        Quadratics coefficient for candidate i-th will be: 
            
        Linear coefficient for candidate i-th will be:
            -(2*(E[1] * R[i,1] + E[2] * R[i,2] + .... + E[m] * R[i,m]) + tau 
        '''
        ## Calculating quadratic coefficients
        quad_variables = []
        for canIndex1 in range(totalCandidates):
            for canIndex2 in range(canIndex1,totalCandidates):
                coefficients = 0
                for skill_Index in range(numberSkill):
                    if (canIndex1 != canIndex2):
                        coefficients += 2 * R[canIndex1][skill_Index] * R[canIndex2][skill_Index]
                    else:
                        coefficients += R[canIndex1][skill_Index] * R[canIndex2][skill_Index] - tau
                quad_variables.append((canIndex1,canIndex2,coefficients))

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
    def buildModel(self , previous_step,E, R,tau,z,c,C,numberSkill,totalCandidates):
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
        # ub = [1 for i in range(numberCandidates)]
        # lb = [0 for i in range(numberCandidates)]
        # self.model.variables.add(names=variables,lb=lb,ub=ub)
        self.model.variables.add(names=variables)
        '''
        Set up variables types is binary
        ['C', 'I', 'B', 'S', 'N'] are integer, binary, semi_continuous, semi_integer
        '''
        self.model.variables.set_types([(elem, self.model.variables.type.binary) for elem in variables])
        '''
        Set up model function, after decomposition
        Quadratics coefficient for candidate i-th will be: 
            R[i,1] ^2 + R[i,2] ^2 + ... + R[i,m] ^ 2
        Linear coefficient for candidate i-th will be:
            -(2*(E[1] * R[i,1] + E[2] * R[i,2] + .... + E[m] * R[i,m]) + tau * ( 2x previous_step[i] - 1 ))
        '''
        ## Calculating quadratic coefficients
        quad_variables = []
        for canIndex1 in range(totalCandidates):
            for canIndex2 in range(canIndex1,totalCandidates):
                coefficients = 0
                for skill_Index in range(numberSkill):
                    if (canIndex1 != canIndex2):
                        coefficients += 2 * R[canIndex1][skill_Index] * R[canIndex2][skill_Index]
                    else:
                        coefficients += R[canIndex1][skill_Index] * R[canIndex2][skill_Index] - tau
                quad_variables.append((canIndex1,canIndex2,coefficients))
        self.model.objective.set_quadratic_coefficients(quad_variables)
        ## Linear Coefficients
        R_copy = np.copy(R)
        for i in range(len(R_copy[0])):
            R_copy[:,i] = [2*element * E[i] for element in R_copy[:,i]]
        linear_coefficients = np.sum(R_copy,axis=1)
        linear_coefficients = [(linear_coefficients[index] + tau * (2*previous_step[index] - 1)) * -1 for index in range(len(linear_coefficients))]
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
        # for i in range(totalCandidates):
        #     linear_constraints_of_quadratics = cplex.SparsePair(ind=[variables[i]], val=[1])
        #     quad_constraints = cplex.SparseTriple(ind1=[variables[i]], ind2=[variables[i]], val=[-1])
        #     quad_sense = 'E'
        #     quad_contraints_rhs = 0
        #     self.model.quadratic_constraints.add(lin_expr=linear_constraints_of_quadratics,
        #                                          quad_expr=quad_constraints,
        #                                          sense=quad_sense,
        #                                          rhs=quad_contraints_rhs)
        # Constraint #4
        self.model.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=variables.copy(), val=c)],
            senses=['L'],
            rhs=[C])

def compute_ObjectiveMDSB(X,numberSkill,E,R):
    sum = 0
    for index in range(numberSkill):
        tmp = 0
        for j in X:
            tmp+= R[j][index]
        sum += (E[index]-tmp) **2
    return sum

def compute_Objective(previousState,X,numberSkill,totalCandidates,tau,E,R):
    sum = 0
    indcies = [index for index in range(len(X)) if (X[index] != 0)]
    for j in range(numberSkill):
        tmp = 0
        for i in range(len(indcies)):
            tmp+=R[indcies[i]][j]
        sum += (E[j] - tmp) * (E[j] - tmp)
    ##
    for index in indcies:
        sum -= tau * (2*previousState[index]-1)
    # print(sum)
    return sum

def solving(option,number_Skill):
    # option 1 is MDSB
    # option 2 is PMDSB
    if option == 1:
        E, R, skillScore, nickNames, tau, z, c, C, numberSkill, totalCandidates = readData(fileName=filenName,
                                                                                           threshHold=threshHold,
                                                                                           numberCandidates=numberCandidates,
                                                                                           numberSkill=number_Skill)
        modelMDSB = MDSB()
        modelMDSB.buildModel( E=E, R=R,tau=tau,z= z,c= c,C= C,numberSkill= numberSkill,totalCandidates= totalCandidates)
        modelMDSB.model.set_log_stream(None)
        modelMDSB.model.set_error_stream(None)
        modelMDSB.model.set_warning_stream(None)
        modelMDSB.model.set_results_stream(None)
        startTime = modelMDSB.model.get_time()
        modelMDSB.model.solve()
        endTime = modelMDSB.model.get_time()
        # if (modelMDSB.model.solution.status.infeasible):
        #     print("Infeasible")
        #     return
        print("Solution for problems using MDSB model is: ")
        solution = modelMDSB.model.solution.get_values()
        print([nickNames[element] for element in range(len(solution)) if (solution[element] != 0)])
        print("Objective Function Values is: ")
        result = modelMDSB.model.solution.get_objective_value()
        for index in range(numberSkill):
            result += E[index]**2
        # print(math.sqrt(result))
        # result = math.sqrt(result)
        # print(result)
        indecies = [element for element in range(len(solution)) if (solution[element] != 0)]
        print(indecies)
        print(math.sqrt(compute_ObjectiveMDSB(indecies,numberSkill,E,R)))
        # DuongsValue = [6,330,456]
        # print(math.sqrt(compute_ObjectiveMDSB(DuongsValue)))
        print("Time Execution")
        print(endTime-startTime)
        print("Number skill: ",numberSkill)

    else:
        E, R, skillScore, nickNames, tau, z, c, C, numberSkill, totalCandidates = readData(fileName=filenName,
                                                                                           threshHold=threshHold,
                                                                                           numberCandidates=numberCandidates,
                                                                                           numberSkill=number_Skill)
        ## random pick
        indices = [i for i in range(totalCandidates)]
        global X_0
        global isSolvedtotalCandidates
        global previousValue
        global objectValue
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
            modelPMDSB.buildModel(X_0,E=E, R=R,tau=tau,z= z,c= c,C= C,numberSkill= numberSkill,totalCandidates= totalCandidates)
            modelPMDSB.model.set_log_stream(None)
            modelPMDSB.model.set_error_stream(None)
            modelPMDSB.model.set_warning_stream(None)
            modelPMDSB.model.set_results_stream(None)
            modelPMDSB.model.solve()
            # if (modelPMDSB.model.solution.status.infeasible):
            #     print("Infeasible")
            #     isSolved = False
            #     break
            X_1 = modelPMDSB.model.solution.get_values()
            currentValue = compute_Objective(X_0,X_1,numberSkill,totalCandidates,tau,E,R)
            objectValue = currentValue
            if (math.fabs(previousValue - currentValue) < epsilon):
                X_0 = X_1
                break
            X_0 = X_1
            previousValue = currentValue

        if (isSolved):
            print("Solution for problems using PMDSB is: ")
            print([index for index in range(len(X_0)) if X_0[index] != 0])
            print("Objective Value is: ")
            print(math.sqrt(objectValue))
            print("Number Skill: ")
            print(numberSkill)

for skillNumber in range(3,39):
    solving(2,skillNumber)

# solving(2,37)

