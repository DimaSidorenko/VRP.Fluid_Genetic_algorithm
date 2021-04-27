#pragma once

#include "Utils.h"
#include "InputData.h"

bool LocalReverseOptimization(std::vector<int>& path, InputData& input, bool pathHasFirst0);
bool StringCrossOptimization(std::vector<std::vector<int>>& paths, InputData& input, bool pathsHaveFirst0 = true);