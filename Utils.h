#pragma once
#include <random>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <numeric>


using namespace std;

class Math {
public:
	int static GenInt(int l, int r);
	double static GenDouble(int l, int r);
	double static GenAllele();
	int static GenAllele(vector<double> &probabilities);
	int static GenRandom(vector<double> &probabilities);
};


enum class InputFormat
{
	file,
	console
};


// It may be use as type of graph
class IntMatrix{
public:
	int N;
	vector<vector<int>> Matrix;

	IntMatrix(vector<vector<int>> Matrix) {
		this->Matrix = Matrix;
		N = Matrix.size();
	}
};
