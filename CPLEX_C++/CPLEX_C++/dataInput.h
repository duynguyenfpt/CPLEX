#pragma once
#include <string>
class dataInput
{
	public:
		// (E,E_before_normalize,R,R_before_normalize,skillScore,nickNames,tau,z,c,C,numberSkill,totalCandidates)
		double* E_before_normalize;
		int* E;
		double** R_before_normalize;
		int** R;
		std::wstring* nickNames;
		double tau;
		double* z;
		double* c;
		double C;
		int numberCandidates;
		int numberSkill;
		int totalCandidates;
		//
		dataInput(){}
		//
		dataInput(int numberSkillCon, int totalCandidatesCon, int numberCandidatesCon) {
			numberSkill = numberSkillCon;
			totalCandidates = totalCandidatesCon;
			numberCandidates = numberCandidatesCon;
			E_before_normalize = new double[numberSkill];
			E = new int[numberSkill];		
			//
			R_before_normalize = new double*[totalCandidates];
			for (int index = 0; index < totalCandidates; index++) {
				R_before_normalize[index] = new double[numberSkill];
			}
			//
			R = new int* [totalCandidates];
			for (int index = 0; index < totalCandidates; index++) {
				R[index] = new int[numberSkill];
			}
			//
			nickNames = new std::wstring[totalCandidates];
			z = new double[numberSkill];
			c = new double[totalCandidates];
		}
};

