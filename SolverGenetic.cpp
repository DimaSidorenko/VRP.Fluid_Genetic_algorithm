#include "SolverGenetic.h"



struct Chromosome2 {
	std::vector<int> Seq;
	double Fitness;

	Chromosome2() {}
	Chromosome2(std::vector<int>& seq, double fitness) :
		Seq(seq),
		Fitness(fitness)
	{}

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

	Chromosome2(std::vector<int>& seq, InputData& input) : Seq(seq) {
		Fitness = CalcFitness(seq, input);
	}

	Chromosome2(std::vector<std::vector<int>>& paths, InputData& input) {
		for (auto& path : paths)
			for (auto x : path)
				Seq.push_back(x);
		if (Seq.size() != input.Size() - 1) {
			Fitness = 0;
			cout << "not correct path in genetic algorithm" << '\n';
		}
		else
			Fitness = CalcFitness(Seq, input);
	}

};

std::mt19937 tmpGen = std::mt19937(1);

Chromosome2 GenRandomChromosome(int n, InputData& input) {
	std::vector<int> seq(n);
	for (int i = 1; i <= n; i++)
		seq[i - 1] = i;
	std::shuffle(seq.begin(), seq.end(), tmpGen);
	return Chromosome2(seq, input);
}

struct Population2 {
	vector<Chromosome2> Chromosomes;
	map<int, int> FitByDelta;
	double Delta;

	Population2(double delta) : Delta(delta) {
		Chromosomes = {};
		FitByDelta = {};
	}

	bool Add(Chromosome2 x) {
		if (Delta == 0) {
			Chromosomes.push_back(x);
			return true;
		}
		if (FitByDelta.count((int)(x.Fitness / Delta)) && FitByDelta[(int)(x.Fitness / Delta)] > 0) {
			return false;
		}
		if (Delta) {
			++FitByDelta[(int)(x.Fitness / Delta)];
		}
		Chromosomes.push_back(x);
		for (int i = (int)Chromosomes.size() - 1; i >= 1; i--) {
			if (Chromosomes[i].Fitness < Chromosomes[i - 1].Fitness) {
				swap(Chromosomes[i], Chromosomes[i - 1]);
			}
		}
		return true;
	}

	bool IsAddible(Chromosome2& x) {
		return (Delta == 0 || !FitByDelta.count((int)(x.Fitness / Delta)) || FitByDelta[(int)(x.Fitness / Delta)] == 0);
	}

	void Del(Chromosome2& x) {
		for (auto it = Chromosomes.begin(); it != Chromosomes.end(); ++it) {
			if (it->Seq == x.Seq) {
				if (Delta) {
					int arg = (int)(x.Fitness / Delta);
					FitByDelta[arg]--;
					if (FitByDelta[arg] == 0)
						FitByDelta.erase(arg);
				}

				Chromosomes.erase(it);
				break;
			}
		}
	}

	void Del(int index) {
		auto it = Chromosomes.begin() + index;
		if (Delta) {
			int arg = (int)(it->Fitness / Delta);
			FitByDelta[arg]--;
			if (FitByDelta[arg] == 0)
				FitByDelta.erase(arg);
		}
		Chromosomes.erase(it);
	}

	int Size() {
		return Chromosomes.size();
	}
};

const int N = 201;

bool used[N];
double prefLen[N];
double dp[N];
int prev[N];

Chromosome2 Crossover(Chromosome2& parent1, Chromosome2& parent2, InputData& input) {
	auto par1 = &parent1;
	auto par2 = &parent2;
	if (Math::GenInt(0, 1)) {
		std::swap(par1, par2);
	}
	int n = par1->Seq.size();
	int tl = Math::GenInt(0, n - 1);
	int tr = Math::GenInt(0, n - 1);
	if (tr < tl) std::swap(tr, tl);

	std::fill(used, used + n + 1, false);
	std::vector<int> child(n, 0);

	for (int i = tl; i <= tr; ++i) {
		child[i] = par1->Seq[i];
		used[par1->Seq[i]] = true;
	}

	int curChildInd = (tr + 1) % n;
	for (int i1 = (tr + 1); i1 < tr + n + 1; i1++) {
		int i = i1;
		while (i >= n)
			i -= n;
		if (!used[par2->Seq[i]]) {
			used[par2->Seq[i]] = true;
			child[curChildInd] = par2->Seq[i];
			curChildInd++;
			if (curChildInd >= n) curChildInd -= n;
		}
	}

	return Chromosome2(child, input);
}



Chromosome2 Mutation(Chromosome2& c, double p, InputData& input) {
	double x = Math::GenDouble(0, 1);
	if (x > p) return c;

	auto seq = c.Seq;
	auto paths = SplitPaths(seq, input);

	bool improved = true;

	while (improved) {
		improved = false;
		improved |= StringCrossOptimization(paths, input, false);
		for (auto& path : paths)
			improved |= LocalReverseOptimization(path, input, false);
	}

	return Chromosome2(paths, input);
}

Population2 InitPopulation(InputData& input, int populationSize, double delta) {
	Population2 population(delta);

	while (population.Size() < populationSize) {

		auto chromosomeCur = GenRandomChromosome(input.Size() - 1, input);

		population.Add(chromosomeCur);
	}

	return population;
}

// Args are: {n, alpha, betta, delta, p, timeLimit}
// Suggested args: {30, 30000, 10000, 0.5, 0.05, 0}
vector<vector<int>> SolverGenetic::Solve(InputData& input, int populationSize) {

	// 2 near chromosomes' fitnesses must differ at list on delta
	double delta = 0.5;

	// if timeLimit is not 0, these 2 variables are useless
	int maxProductiveIters = 30000; // alpha
	int maxNonImproveIters = 10000; // betta

	double mutationProb = 0.05; // p
	// if timeLimit is 0, then cycle is limited with maxProductiveIters and maxNonImproveIters
	int timeLimit = 5;

	double startT = clock();

	Population2 population = InitPopulation(input, populationSize, delta);

	int productiveIters = 0;
	int nonImproveIters = 0;
	int numIters = 0;

	int n = population.Size();

	while (timeLimit ? clock() - startT < timeLimit * CLOCKS_PER_SEC :
		productiveIters < maxProductiveIters && nonImproveIters < maxNonImproveIters) {

		numIters++;

		assert(population.Size() == n);

		int parent1ind = std::min(Math::GenInt(0, n - 1), Math::GenInt(0, n - 1));
		int parent2ind = std::min(Math::GenInt(0, n - 1), Math::GenInt(0, n - 1));
		while (parent2ind == parent1ind) {
			parent2ind = std::min(Math::GenInt(0, n - 1), Math::GenInt(0, n - 1));
		}
		if (parent2ind < parent1ind) std::swap(parent1ind, parent2ind);

		auto child = Crossover(population.Chromosomes[parent1ind], population.Chromosomes[parent2ind], input);
		auto mutatedChild = Mutation(child, mutationProb, input);

		int replacedInd = Math::GenInt(n / 2, n - 1);

		auto deletedChromosome = population.Chromosomes[replacedInd];

		population.Del(replacedInd);

		if (!population.IsAddible(mutatedChild) && population.IsAddible(child)) {
			mutatedChild = child;
		}

		nonImproveIters++;

		if (nonImproveIters % 50 == 0) {
			//std::cout << "iters: " << nonImproveIters << " " << productiveIters << std::endl;
		}

		if (population.IsAddible(mutatedChild)) {
			if (mutatedChild.Fitness < population.Chromosomes[0].Fitness) {
				nonImproveIters = 0;
				//std::cout << "New best fitness: " << mutatedChild.Fitness << std::endl;
			}
			productiveIters++;

			population.Add(mutatedChild);
		}
		else {
			population.Add(deletedChromosome);
		}
	}

	cout << "Genetic :\n";
	cout << "Fitness = " << population.Chromosomes[0].Fitness << '\n';
	cout << "Iteration Count = " << numIters << "\n\n";

	return SplitPaths(population.Chromosomes[0].Seq, input);
}
