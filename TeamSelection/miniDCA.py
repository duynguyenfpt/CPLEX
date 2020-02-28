import cplex
import numpy as np
import  pandas as pd

filenName = 'southeast-asia.xlsx'
threshHold = 0
numberCandidates = 3

def readData(fileName,threshHold,numberCandidates,numberSkill):
    dataset = pd.read_excel(fileName)
    '''
    R[i][j] is skill-depth of candidates i-th in skill j-th
    '''
    R = dataset.iloc[:4,1:numberSkill].values
    nickNames = dataset.iloc[:4,0].values
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

E,R,skillScore,nickNames,tau,z,c,C,numberSkill,totalCandidates = readData(fileName=filenName,threshHold=threshHold,numberCandidates=numberCandidates,numberSkill=2)