#include "Optimizations.h"

bool LocalReverseOptimization(std::vector<int>& path, InputData& input, bool pathHasFirst0) {
	path.push_back(0);
	if (!pathHasFirst0) {
		path.push_back(0);
		rotate(path.begin(), path.begin() + path.size() - 1, path.end());
	}
	// we have path: 0 4 7 1 3 0

	bool improvedSolution = false;

	while (true) {
		bool improved = false;

		// try to reverse [i, j]
		for (size_t i = 1; i + 1 < path.size(); ++i) {
			for (size_t j = i + 1; j + 1 < path.size(); ++j) {
				double curDist = input.Distance(path[i - 1], path[i]) +
					input.Distance(path[j], path[j + 1]);
				double anotherDist = input.Distance(path[i - 1], path[j]) +
					input.Distance(path[i], path[j + 1]);
				if (anotherDist < curDist) {
					reverse(path.begin() + i, path.begin() + j + 1);
					improved = true;
					improvedSolution = true;
				}
			}
		}

		if (!improved)
			break;
	}
	if (!pathHasFirst0) {
		rotate(path.begin(), path.begin() + 1, path.end());
		assert(path.back() == 0);
		path.pop_back();
	}
	assert(path.back() == 0);
	path.pop_back();

	if (!pathHasFirst0) {
		assert(find(path.begin(), path.end(), 0) == path.end());
	}

	return improvedSolution;
}

// Recalculates prefix and sufix sums of distances of path
void recalcPS(std::vector<double>& pref, std::vector<double>& suf, std::vector<int>& path, InputData& input) {
	if (!path.size()) {
		pref.resize(0);
		suf.resize(0);
		return;
	}
	pref.resize(path.size());
	suf.resize(path.size());
	pref[0] = 0;
	for (int i = 1; i < (int)pref.size(); i++) {
		pref[i] = pref[i - 1] + input.Distance(path[i - 1], path[i]);
	}
	suf[(int)path.size() - 1] = 0;
	for (int i = (int)path.size() - 2; i >= 0; i--) {
		suf[i] = suf[i + 1] + input.Distance(path[i], path[i + 1]);
	}
}


void MakeStringCross(std::vector<std::vector<int>>& paths, std::vector<std::vector<double>>& pref,
	std::vector<std::vector<double>>& suf, int i1, int i2, int p1, int p2, InputData& input) {
	auto tmp = paths[p1];
	paths[p1].resize(i1 - i2 + (int)paths[p2].size());
	for (int j = i2 + 1; j < (int)paths[p2].size(); ++j)
		paths[p1][j - i2 + i1] = paths[p2][j];
	paths[p2].resize(i2 - i1 + (int)tmp.size());
	for (int j = i1 + 1; j < (int)tmp.size(); ++j)
		paths[p2][j - i1 + i2] = tmp[j];
	recalcPS(pref[p1], suf[p1], paths[p1], input);
	recalcPS(pref[p2], suf[p2], paths[p2], input);
}


// Paths have to be in format {{0, 1, 2, 5}, {0, 6, 4}} or {{1, 2, 5}, {6, 4}}
bool StringCrossOptimization(std::vector<std::vector<int>>& paths, InputData& input, bool pathsHaveFirst0) {
	for (auto& path : paths)
		path.emplace_back(0);
	if (!pathsHaveFirst0) {
		for (auto& path : paths)
			path.insert(path.begin(), { 0 });
	}

	std::vector<std::vector<double>> pref(paths.size());
	std::vector<std::vector<double>> suf(paths.size());

	for (int p = 0; p < (int)paths.size(); ++p) {
		recalcPS(pref[p], suf[p], paths[p], input);
	}

	bool improvedSolution = false;

	for (int p1 = 0; p1 < (int)paths.size(); ++p1) {
		for (int p2 = p1 + 1; p2 < (int)paths.size(); ++p2) {
			for (int i1 = 1; i1 < (int)paths[p1].size() - 2; ++i1) {
				for (int i2 = 1; i2 < (int)paths[p2].size() - 2; ++i2) {
					//assert(paths[p1].size() == pref[p1].size() && paths[p1].size() == suf[p1].size());
					//assert(paths[p2].size() == pref[p2].size() && paths[p2].size() == suf[p2].size());
					double pref1 = pref[p1][i1];
					double suf1 = suf[p1][i1 + 1];
					double pref2 = pref[p2][i2];
					double suf2 = suf[p2][i2 + 1];
					double edge1 = input.Distance(paths[p1][i1], paths[p1][i1 + 1]);
					double edge2 = input.Distance(paths[p2][i2], paths[p2][i2 + 1]);
					double crossEdge1 = input.Distance(paths[p1][i1], paths[p2][i2 + 1]);
					double crossEdge2 = input.Distance(paths[p2][i2], paths[p1][i1 + 1]);

					if (edge1 + edge2 > crossEdge1 + crossEdge2) {
							MakeStringCross(paths, pref, suf, i1, i2, p1, p2, input);
							improvedSolution = true;
					}
				}
			}
		}
	}

	if (!pathsHaveFirst0) {
		for (auto& path : paths) {
			assert(path[0] == 0);
			path.erase(path.begin());
		}
	}

	for (auto& path : paths) {
		assert(path.back() == 0);
		path.pop_back();
	}

	if (!pathsHaveFirst0) {
		for (auto& path : paths) {
			assert(find(path.begin(), path.end(), 0) == path.end());
		}
	}

	return improvedSolution;
}