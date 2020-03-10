#pragma once
#include <string>
#include <vector>
class dataInput
{
	public:
		// (E,E_before_normalize,R,R_before_normalize,skillScore,nickNames,tau,z,c,C,numberSkill,totalCandidates)
		std::vector<double> E_before_normalize;
		std::vector<int> E;
		std::vector<std::vector<double>> R_before_normalize;
		std::vector<std::vector<int>> R;
		std::vector<std::string> nickNames;
		double tau;
		std::vector<double> z_before_normalize;
		std::vector<double> z;
		double* c;
		double C;
		int numberCandidates;
		int numberSkill;
		int totalCandidates;
		//
		dataInput(){}
		//
		dataInput(int numberSkillCon, int numberCandidatesCon, int totalCandidatesCon) {
			numberSkill = numberSkillCon;
			totalCandidates = totalCandidatesCon;
			numberCandidates = numberCandidatesCon;	
			//
			
			R_before_normalize = std::vector<std::vector<double>>(totalCandidates,std::vector<double>(numberSkill));
			R = std::vector<std::vector<int>>(totalCandidates,std::vector<int>(numberSkill));
			//
		}
};

