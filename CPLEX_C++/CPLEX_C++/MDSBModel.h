#pragma once
#include <ilcplex/ilocplex.h>
#include "dataInput.h"
#include <vector>
using namespace std;

class MDSBModel
{
	public:
		void buildModel() {
			IloEnv env;
			try {
				IloModel model(env);
				IloNumVarArray vars(env);
				for (int i = 0; i < 2; i++) {
					vars.add(IloNumVar(env));
				}
				IloExpr linearExpr(env);
				/*linearExpr += vars[0] + vars[1];*/
				for (int i = 0; i < 2; i++) {
					linearExpr += vars[i];
				}
				model.add(linearExpr == 10);
				//
				IloExpr linearExpr2(env);
				linearExpr2 += vars[0] - vars[1];
				model.add(linearExpr2 == 1);
				IloExpr objExpr(env);
				objExpr += vars[0] * vars[0] + vars[1]*vars[1] + vars[0]*vars[1] + vars[1]*vars[0];
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
			}
			catch (IloException & e) {
				cerr << "Concert exception caught: " << e << endl;
			}
			catch (...) {
				cerr << "Unknown exception caught" << endl;
			}
				env.end();
		}

		void buildModel(dataInput input) {
			IloEnv env;
			try {
				IloModel model(env);
				IloNumVarArray vars(env);
				// set up parameters is boolean
				for (int index = 0; index < input.totalCandidates; index++) {
					vars.add(IloNumVar(env,ILOBOOL));
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
				expressions+= vars[1];
				model.add(expressions >= 0);
				//
				for (int index = 0; index < input.numberSkill; index++) {
					IloExpr expressions;
					for (int member = 0; member < input.totalCandidates -1 ; member++) {
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
						coefficient += -2* input.E[skill]* input.R[index][skill];
					}
					objExpr += coefficient * vars[index];
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
			}
			catch (IloException & e) {
				cerr << "Concert exception caught: " << e << endl;
			}
			catch (...) {
				cerr << "Unknown exception caught" << endl;
			}
			env.end();
		}
};

