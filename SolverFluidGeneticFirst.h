#pragma once
#include "Utils.h"
#include "InputData.h"


vector<vector<int>> SplitPaths(vector<int>& path, InputData& input);

class SolverFluidGeneticFirst
{
public:
	static double diversityRate;
	static double individualLR;
	static double globalLR;

	static bool Solve(InputData& input, int populationSize, int cntIteration);
};

