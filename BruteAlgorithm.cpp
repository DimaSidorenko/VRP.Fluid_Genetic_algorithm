#include "BruteAlgorithm.h"
#include "SolverFluidGeneticSecond.h"


static double CalcFitness(vector<int>& seq, InputData& input) {
	if (seq.size() + 1 != input.Size()) {
		return -1;
	}

	seq.insert(seq.begin(), { 0 });
	int n = (int)(seq.size());

	double* prefLen = new double[n];
	fill(prefLen, prefLen + n, 0);

	for (int i = 1; i < n; ++i) {
		prefLen[i] = prefLen[i - 1] + input.Distance(seq[i - 1], seq[i]);
	}

	double* dp = new double[n];
	fill(dp, dp + n, INF);

	dp[0] = 0;
	for (int v = 0; v + 1 < n; ++v) {
		if (dp[v] == INF) {
			continue;
		}
		double len = 0;
		for (int u = v + 1; u < n; ++u) {
			// relaxing using edge (v, u)
			len = prefLen[u] - prefLen[v + 1] + input.Distance(0, seq[v + 1]) + input.Distance(seq[u], 0);
			dp[u] = min(dp[u], dp[v] + len);
		}
	}

	delete[] prefLen;
	auto answer = dp[n - 1];
	delete[] dp;
	seq.erase(seq.begin());
	return answer;
}


bool BruteAlgorithm::Solve(InputData &input) {
	int n = input.Size() - 1;

	vector<int> seq(n);
	for (int i = 0; i < n; i++) {
		seq[i] = i + 1;
	}

	vector<int> best_seq;
	int best_seq_fitness = INT_MAX;

	do {
		int fitness = CalcFitness(seq, input);

		if (fitness < best_seq_fitness) {
			best_seq_fitness = fitness;
			best_seq = seq;
		}
	} while (next_permutation(seq.begin(), seq.end()));

	cout << "Brute Algorithm Answer Fitness = " << best_seq_fitness << '\n';
	
	return false;
}