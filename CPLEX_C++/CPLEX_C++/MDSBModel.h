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
			objExpr += vars[0] * vars[0] + vars[1] * vars[1] + vars[0] * vars[1] + vars[1] * vars[0];
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
			IloBoolVarArray vars(env);
			// set up parameters is boolean
			for (int index = 0; index < input.totalCandidates; index++) {
				//IloInt n = 2;
				//IloNumArray possibleValues(env, 2);
				//possibleValues.add(IloInt(0));
				//possibleValues.add(IloInt(1));
				vars.add(IloBoolVar(env));
				/*std::cout <<"type: " << vars[index].getType() << endl;*/
			}
			//
			//std::cout <<"type: " << vars[0] << endl;
			//std::cout <<"type: " << vars[1] << endl;
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

			//
			for (int index = 0; index < input.numberSkill; index++) {
				IloExpr expressions(env);
				for (int member = 0; member < input.totalCandidates - 1; member++) {
					expressions += vars[member] * input.R[member][index];
				}
				model.add(expressions >= input.z[index]);
			}
			// quadratics expression
			IloExpr objExpr(env);
			for (int index = 0; index < input.totalCandidates; index++) {
				for (int index2 = index; index < input.totalCandidates; index2++) {
					if (index2 == input.totalCandidates) {
						break;
					}
					int coefficients = 0;
					for (int skill = 0; skill < input.numberSkill; skill++) {
						if (index == index2) {
							coefficients += (input.R[index][skill] * input.R[index][skill]);
						}
						else {
							coefficients += 2 * input.R[index][skill] * input.R[index2][skill];
						}
					}
					/*std::cout << index << "-" << index2 << "-" << coefficients << std::endl;*/
					objExpr += coefficients * vars[index] * vars[index2];
				}
			}
			//
			for (int index = 0; index < input.totalCandidates; index++) {
				int coefficient = 0;
				for (int skill = 0; skill < input.numberSkill; skill++) {
					coefficient += -2 * input.E[skill] * input.R[index][skill];
				}
				/*std::cout << coefficient << " ";*/
				objExpr += coefficient * vars[index];
			}
			//
			IloObjective obj = IloMinimize(env, objExpr);
			model.add(obj);
			//
			IloCplex cplex(model);
			cplex.setParam(IloCplex::Param::OptimalityTarget, 3);
			cplex.setOut(env.getNullStream());
			std::cout << "start solving" << std::endl;
			IloTimer timer(env);
			IloNum timeStart = timer.getTime();
			if (!cplex.solve()) {
				IloNumArray vals(env);
				cplex.getValues(vals, vars);
				std::cout << "end solving" << std::endl;
				throw(-1);
			}
			IloNum timeEnd = timer.getTime();
			std::cout << "end solving" << std::endl;
			std::cout << "time for solving: " << timeEnd - timeStart << std::endl;
			IloNumArray vals(env);
			IloNum objValuie = cplex.getObjValue();
			env.out() << "Solution status = " << cplex.getStatus() << endl;
			env.out() << "Solution value = " << objValuie << endl;
			cplex.getValues(vals, vars);
			env.out() << "Values = " << vals << endl;
		}
		catch (IloAlgorithm::CannotExtractException & e) {
			IloExtractableArray& failed = e.getExtractables();
			std::cerr << "Failed to extract:" << std::endl;
			for (IloInt i = 0; i < failed.getSize(); ++i)
				std::cerr << "\t" << failed[i] << std::endl;
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

