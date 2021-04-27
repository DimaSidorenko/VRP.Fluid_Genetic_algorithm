#pragma once
#include "Utils.h"
#include "InputData.h"
#include <map>
#include <algorithm>
#include "Optimizations.h"
#include "SolverFluidGeneticSecond.h"

class SolverGenetic {
public:
	static vector<vector<int>> Solve(InputData &input, int populationSize);
};

