#pragma once
#include <string>
#include <iostream>
#include "libxl.h"
#include <vector>
#include <map>
#include <windows.h>
#include <string>
#include <libxl.h>
#include <iostream>
#include <WinNls.h>
#include "dataInput.h"

using namespace libxl;
class dataInputHandler
{
public:
	std::string fileName;
	int numberSkill;
	int numberCandidates;
	int totalCandidates;
	dataInput data;

	dataInputHandler(std::string fileNameCon, int numberSkillCon, int numberCandidatesCon, int totalCandidatesCon) {
		fileName = fileNameCon;
		numberCandidates = numberCandidatesCon;
		numberSkill = numberSkillCon;
		totalCandidates = totalCandidatesCon;
		data = dataInput(numberSkill, totalCandidates, numberCandidates);
	}
	template <typename T>
	void printOutMatrix(T** R) {
		for (int index = 0; index < totalCandidates; index++) {
			for (int index2 = 0; index2 < numberSkill; index2++) {
				std::cout << R[index][index2] << "   ";
			}
			std::cout << std::endl;
		}
	}
	// generic for both double and c++
	template <typename T>
	void sortByCollumn(T** R) {
		for (int skill_ith = 0; skill_ith < numberSkill; skill_ith++) {
			// sorting each row using selection sort
			for (int candidates_i = 0; candidates_i < totalCandidates - 1; candidates_i++) {
				for (int candidates_j = candidates_i + 1; candidates_j < totalCandidates; candidates_j++) {
					if (R[candidates_i][skill_ith] > R[candidates_j][skill_ith]) {
						double tmp = R[candidates_i][skill_ith];
						R[candidates_i][skill_ith] = R[candidates_j][skill_ith];
						R[candidates_j][skill_ith] = tmp;
					}
				}
			}
		}
	}
	template <typename T>
	std::map<std::string, int>* normalizingData(T** R) {
		// normalizing equal values
		std::map<std::string, int>* listMap = new std::map<std::string, int>[numberSkill];
		for (int skill_ith = 0; skill_ith < numberSkill; skill_ith++) {
			T value = R[0][skill_ith];
			int ranking = 1;
			listMap[skill_ith][std::to_string(value)] = ranking;
			for (int candidates_i = 1; candidates_i < totalCandidates; candidates_i++) {
				if (R[candidates_i][skill_ith] > value) {
					value = R[candidates_i][skill_ith];
					ranking++;
					listMap[skill_ith][std::to_string(value)] = ranking;
				}
			}
		}
		//
		printOutMatrix(data.R_before_normalize);
		//
		for (int index = 0; index < totalCandidates; index++) {
			for (int index2 = 0; index2 < numberSkill; index2++) {
				data.R[index][index2] = listMap[index2][std::to_string(data.R_before_normalize[index][index2])];
			}
		}
		// return ranking
		printOutMatrix<int>(data.R);
		//
		return listMap;
	}

	void readData() {
		Book* book = xlCreateBook();
		std::wstring wide_string = std::wstring(fileName.begin(), fileName.end());
		const wchar_t* result = wide_string.c_str();
		//
		if (book->load(L"southeast-asia.xls"))
		{
			Sheet* sheet = book->getSheet(0);
			if (sheet)
			{
				//std::cout << sheet->firstRow() << std::endl;
				//std::cout << sheet->lastRow();

				for (int row = sheet->firstRow(); row < sheet->lastRow(); ++row)
				{
					if (row == totalCandidates + 1) { break; };
					// first collumn is string nickname
					int firstCol = sheet->firstCol();
					const wchar_t* s = sheet->readStr(row, firstCol);
					/*std::wcout << s << std::endl;*/
					data.nickNames[row - 1] = std::wstring(s);
					//
					for (int col = firstCol + 1; col < sheet->lastCol(); ++col)
					{
						if (col == numberSkill + 1) { break; };
						CellType cellType = sheet->cellType(row, col);
						double d = sheet->readNum(row, col);
						data.R_before_normalize[row - 1][col - 1] = d;
					}
				}
				// int R_size = sizeof(data.R)/sizeof(data.R[0]);
				// create a new copy of R_before_normalize for normalizing
				double** R_before_normalize_copy = new double* [totalCandidates];
				for (int index = 0; index < totalCandidates; index++) {
					R_before_normalize_copy[index] = new double[numberSkill];
				}
				//
				//std::cout << sizeof(data.R_before_normalize) << std::endl;
				//std::cout << sizeof(data.R_before_normalize[0]) << std::endl;
				memcpy(R_before_normalize_copy, data.R_before_normalize, sizeof(data.R_before_normalize) * sizeof(data.R_before_normalize[0]));
				std::cout << data.R_before_normalize << std::endl;
				std::cout << R_before_normalize_copy << std::endl;
				printOutMatrix<double>(data.R_before_normalize);
				sortByCollumn<double>(R_before_normalize_copy);
				//check if address of two R_before_normalize and its copy are the same
				// printOutMatrix(R_before_normalize_copy);
				std::map< std::string, int>* listMap = normalizingData(R_before_normalize_copy);

			}
		}
		else {
			std::cout << "Unsuccessful";
		}

		book->release();
	}
};


