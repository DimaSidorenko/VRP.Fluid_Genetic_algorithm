#include <iostream>
#include <string>
#include "Utils.h"
#include "SolverFluidGeneticfirst.h"
#include "SolverFluidGeneticSecond.h"
#include "BruteAlgorithm.h"
#include "InputData.h"
#include "SolverGenetic.h"

#pragma warning(disable : 4996)

using namespace std;

const char* fileName = "input.txt";
//InputFormat inputFormat = InputFormat::file;


double SolverFluidGeneticFirst::diversityRate = 0.0005;
double SolverFluidGeneticFirst::individualLR = 0.05;
double SolverFluidGeneticFirst::globalLR = 0.05;

double SolverFluidGeneticSecond::diversityRate = 0.0005;
double SolverFluidGeneticSecond::individualLR = 0.05;
double SolverFluidGeneticSecond::globalLR = 0.05;


void test() {
	InputData input(80);	
	input.WriteInFile(fileName);

	//SolverFluidGeneticFirst::Solve(input, 100, 1000);
	SolverFluidGeneticSecond::Solve(input, 60, 1000);
	SolverGenetic::Solve(input, 10);
	
	//BruteAlgorithm::Solve(input);
}


int main() {
	srand(time(0));
	cout.setf(ios::fixed);
	cout.precision(8);
	//read();
	
	double startT = clock();	
	
	bool test_mode = true;

	if (test_mode) {
		test();
	}

	cout << "All Time Work  = " << (clock() - startT) / CLOCKS_PER_SEC << '\n';
}




