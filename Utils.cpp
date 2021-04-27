#include "Utils.h"

static mt19937 randGen;

int Math::GenInt(int l, int r) {
	uniform_int_distribution<int> uid(l, r);
	return uid(randGen);
}

double Math::GenDouble(int l, int r) {
	uniform_real_distribution<double> urd(l, r);
	return urd(randGen);
}



double Math::GenAllele() {
	uniform_real_distribution<double> urd(0.0, 1.0);
	return urd(randGen);
}


int Math::GenRandom(vector<double>& probabilities) {
	random_device device;
	mt19937 engine(device()); // Seed the random number engine;
	discrete_distribution<> dist(probabilities.begin(), probabilities.end()); // Create the distribution.

	int result = dist(engine);

	return result;
}


int Math::GenAllele(vector<double> &probabilities) {
	for (auto prob : probabilities) {
		assert(prob <= 1.0);
	}

	int allele = GenRandom(probabilities) + 1;

	return allele;
}



