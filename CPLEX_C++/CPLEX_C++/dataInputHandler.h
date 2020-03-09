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

using namespace libxl;
class dataInputHandler
{
public:
	std::string fileName;
	int numberSkill;
	int numberCandidates;

	dataInputHandler(std::string fileNameCon, int numberSkillCon, int numberCandidatesCon) {
		this->fileName = fileNameCon;
		this->numberCandidates = numberCandidatesCon;
		this->numberSkill = numberSkillCon;
	}

	std::map<std::string, int>* normalizingData(double R[500][37]) {
		// indexing
		std::vector<int> normailzedData[500][37];
		// ranking score of each collumn
		for (int skill_ith = 0; skill_ith < numberSkill; skill_ith++) {
			// sorting each row using selection sort
			for (int candidates_i = 0; candidates_i < 499; candidates_i++) {
				for (int candidates_j = 0; candidates_j < 500; candidates_j) {
					if (R[candidates_i][skill_ith] > R[candidates_j][skill_ith]) {
						double tmp = R[candidates_i][skill_ith];
						R[candidates_i][skill_ith] = R[candidates_j][skill_ith];
						R[candidates_j][skill_ith] = tmp;
					}
				}
			}
		}
		// normalizing equal values
		std::map<std::string, int> listMap[37];
		for (int skill_ith = 0; skill_ith < numberSkill; skill_ith++) {
			double value = R[0][skill_ith];
			int ranking = 1;
			listMap[skill_ith][std::to_string(value)] = ranking;
			for (int candidates_i = 1; candidates_i < 500; candidates_i++) {
				if (R[candidates_i][skill_ith] > value) {
					value = R[candidates_i][skill_ith];
					ranking++;
					listMap[skill_ith][std::to_string(value)] = ranking;
				}
			}
		}
		// return ranking
		return listMap;
	}

	void readData() {
		Book* book = xlCreateBook();
		std::wstring wide_string = std::wstring(fileName.begin(), fileName.end());
		const wchar_t* result = wide_string.c_str();
		if (book->load(L"southeast-asia.xls"))
		{
			Sheet* sheet = book->getSheet(0);
			if (sheet)
			{
				//std::cout << sheet->firstRow() << std::endl;
				//std::cout << sheet->lastRow();

				for (int row = sheet->firstRow(); row < sheet->lastRow(); ++row)
				{

					if (row == 500) { break; };
					for (int col = sheet->firstCol(); col < sheet->lastCol(); ++col)
					{
						CellType cellType = sheet->cellType(row, col);
						std::wcout << "(" << row << ", " << col << ") = ";
						if (sheet->isFormula(row, col))
						{
							const wchar_t* s = sheet->readFormula(row, col);
							std::wcout << (s ? s : L"null") << " [formula]";
						}
						else
						{
							switch (cellType)
							{
							case CELLTYPE_EMPTY: std::wcout << "[empty]"; break;
							case CELLTYPE_NUMBER:
							{
								double d = sheet->readNum(row, col);
								std::wcout << d << " [number]";
								break;
							}
							case CELLTYPE_STRING:
							{
								const wchar_t* s = sheet->readStr(row, col);
								std::wcout << (s ? s : L"null") << " [string]";
								break;
							}
							case CELLTYPE_BOOLEAN:
							{
								bool b = sheet->readBool(row, col);
								std::wcout << (b ? "true" : "false") << " [boolean]";
								break;
							}
							case CELLTYPE_BLANK: std::wcout << "[blank]"; break;
							case CELLTYPE_ERROR: std::wcout << "[error]"; break;
							}
						}
						std::wcout << std::endl;
					}
				}
			}
		}
		else {
			std::cout << "Unsuccessful";
		}

		book->release();
	}
};


