#pragma once

#include <ilcplex/ilocplex.h>
#include "dataInput.h"
#include <vector>
#include <stdlib.h>   
#include <limits>

using namespace std;

class PMDSB
{
public:
	IloNumArray solvingModel(IloNumArray previousState, dataInput input) {
		IloNumArray defaultValue;
		for (int index = 0; index < input.totalCandidates; index++) {
			defaultValue.add(0);
		}
		//
		IloEnv env;
		try {
			IloModel model(env);
			IloNumVarArray vars(env);
			// set up parameters is boolean
			for (int index = 0; index < input.totalCandidates; index++) {
				vars.add(IloNumVar(env, ILOBOOL));
			}
			// constraint sum of variables equals to numberCandidates
			IloExpr constraint1(env);
			for (int index = 0; index < input.totalCandidates; index++) {
				constraint1 += vars[index];
			}
			model.add(constraint1 == input.numberCandidates);
			//
			// IloObjective obj = IloMinimize(env, objExpr);
			// model.add(obj);
			//

			IloExpr expressions;
			expressions += vars[1];
			model.add(expressions >= 0);
			//
			for (int index = 0; index < input.numberSkill; index++) {
				IloExpr expressions;
				for (int member = 0; member < input.totalCandidates - 1; member++) {
					expressions += vars[member] * input.R[member][index];
				}
				model.add(expressions >= input.z[index]);
			}
			// quadratics expression
			IloExpr objExpr;
			for (int index = 0; index < input.totalCandidates; index++) {
				for (int index2 = index; index < input.totalCandidates; index2++) {
					int coefficients = 0;
					for (int skill = 0; skill < input.numberSkill; skill++) {
						if (index == index2) {
							coefficients += (input.R[index][skill] * input.R[index][skill]);
						}
						else {
							coefficients += 2 * input.R[index][skill] * input.R[index2][skill];
						}
					}
					objExpr += coefficients * vars[index] * vars[index2];
				}
			}
			//
			for (int index = 0; index < input.totalCandidates; index++) {
				int coefficient = 0;
				for (int skill = 0; skill < input.numberSkill; skill++) {
					coefficient += -2 * input.E[skill] * input.R[index][skill];
				}
				objExpr += coefficient * vars[index] - input.tau * (2 * previousState[index] - 1);
			}
			//
			IloObjective obj = IloMinimize(env, objExpr);
			model.add(obj);
			//
			IloCplex cplex(model);
			if (!cplex.solve()) {
				env.error() << "Failed to optimize LP." << endl;
				throw(-1);
			}
			IloNumArray vals(env);
			env.out() << "Solution status = " << cplex.getStatus() << endl;
			env.out() << "Solution value = " << cplex.getObjValue() << endl;
			cplex.getValues(vals, vars);
			env.out() << "Values = " << vals << endl;
			return vals;
		}
		catch (IloException & e) {
			cerr << "Concert exception caught: " << e << endl;
			return defaultValue;
		}
		catch (...) {
			cerr << "Unknown exception caught" << endl;
			return defaultValue;
		}
		env.end();
		return defaultValue;
	}

	bool isUnsolve(IloNumArray aArray) {
		for (int index = 0; index < aArray.getSize(); index++) {
			if (aArray[index] != 0) {
				return false;
			}
		}
		return true;
	}

	double getObjectivePMDSB(IloNumArray previousSate, IloNumArray currentSate, dataInput input) {
		double result = 0;
		for (int skill; skill < input.numberSkill; skill++) {
			double tmp = 0;
			for (int can = 0; can < input.totalCandidates; can++) {
				tmp += input.R[can][skill] * currentSate[can];
			}
			result += (input.E[skill] - tmp) * (input.E[skill] - tmp);
		}
		//
		for (int can = 0; can < input.totalCandidates; can++) {
			result += input.tau * (2 * previousSate[can] - 1) * currentSate[can];
		}
		return result;
	}

	void solving(dataInput input) {
		// create a random state
		bool isSolved = true;
		std::vector<int> indicies;
		for (int index = 0; index < input.totalCandidates; index++) {
			indicies.push_back(index);
		}
		//
		IloNumArray previousSate;
		for (int index = 0; index < input.totalCandidates; index++) {
			previousSate.add(0);
		}
		for (int index = 0; index < input.numberCandidates; index++) {
			int randomIndex = rand() * indicies.size;
			previousSate[randomIndex] = 1;
			indicies.erase(indicies.begin() + randomIndex - 1);
		}
		// looping for solving
		// print out data
		double epsilon = 0.00001;
		double previousValue = numeric_limits<double>::max();
		double currentValue = 0;
		while (true) {
			IloNumArray currentSate = solvingModel(previousSate, input);
			if (isUnsolve(currentSate)) {
				std::cout << "No Solution" << std::endl;
				isSolved = false;
			}
			else {
				currentValue = getObjectivePMDSB(previousSate, currentSate, input);
				previousSate = currentSate;
				if (fabs(currentValue - previousValue) < epsilon) {
					previousValue = currentValue;
					break;
				}
				previousValue = currentValue;
			}
		}
		//
		if (isSolved) {
			previousSate.getEnv().out << "Values: " << previousSate;
		}
	}
};
