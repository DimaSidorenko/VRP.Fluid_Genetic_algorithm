#include "SolverFluidGeneticFirst.h"


bool enable_blueprint1 = true;

struct Gene1 {
	int N; // cnt alleles
	vector<double> probabilities; // probabilities for gene to become corresponding allel 

	Gene1() {
		N = rand();
		probabilities.resize(N);
		for (size_t i = 0; i < probabilities.size(); i++)
		{
			probabilities[i] = 1 / double(N);
		}
	}

	Gene1(int size) {
		N = size;
		probabilities.resize(N);
		for (size_t i = 0; i < probabilities.size(); i++)
		{
			probabilities[i] = 1 / double(N);
		}
	}

	int GetAllel() {
		return Math::GenAllele(probabilities);
	}

	int Size() {
		return N;
	}
};


struct Chromosome1 {
	int N; // size of chromosome
	vector<Gene1> Genes; // probabilities for sequences

	Chromosome1() {
		N = 0;
		Genes = {};
	}

	Chromosome1(int N) {
		this->N = N;

		Genes.resize(N);
		for (int i = 0; i < N; i++) {
			Genes[i] = Gene1(N - i);
		}
	}

	int Size() {
		return N;
	}
};

//it may be last chromosome of popupation or have some flag in population to detect that it is blueprint
Chromosome1 blueprint;


struct Individual1 {
	Chromosome1 Chr; // probabilities for Sequences
	vector<int> EncodedSeq;
	vector<int> Seq; // sequences
	int Fitness;


	Individual1() {
		Chr = Chromosome1(0);
		EncodedSeq = {};
		Seq = {};
		Fitness = INF;
	}

	Individual1(int size) {
		Chr = Chromosome1(size);
		EncodedSeq.resize(size);
		Seq.resize(size);
		Fitness = INF;
	}

	Individual1(Chromosome1& chr, InputData& input) {
		Chr = chr;
		EncodedSeq = BornAnIndividual(chr);
		Seq = Decode(EncodedSeq);
		Fitness = CalcFitness(Seq, input);
	}

	static int CalcFitness(vector<int>& seq, InputData& input) {
		if (seq.size() + 1 != input.Size()) {
			return -1;
		}

		seq.insert(seq.begin(), { 0 });
		int n = (int)(seq.size());

		int* prefLen = new int[n];
		fill(prefLen, prefLen + n, 0);

		for (int i = 1; i < n; ++i) {
			prefLen[i] = prefLen[i - 1] + input.Distance(seq[i - 1], seq[i]);
		}

		int* dp = new int[n];
		fill(dp, dp + n, INF);

		dp[0] = 0;
		for (int v = 0; v + 1 < n; ++v) {
			if (dp[v] == INF) {
				continue;
			}
			int len = 0;
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

	vector<int> BornAnIndividual(Chromosome1& chr) {
		int n = chr.Size();
		vector<int> seq(n, 0);


		double gLR = SolverFluidGeneticFirst::globalLR;


		for (int i = 0; i < n; i++) {
			if (enable_blueprint1) {
				for (int j = 0; j < chr.Genes[i].Size(); j++) {
					chr.Genes[i].probabilities[j] = gLR * blueprint.Genes[i].probabilities[j] + (1 - gLR) * chr.Genes[i].probabilities[j];
				}
			}

			seq[i] = chr.Genes[i].GetAllel();
		}

		return seq;
	}

	vector<int> Decode(vector<int>& seq) {
		int n = seq.size();
		vector<bool> used(n, false);
		vector<int> res;


		//iota(res.begin(), res.end(), 1);

		//It may be possible to solve this in O(n * logn) time.
		for (int i = 0; i < n; i++) {
			int cnt = 0;
			int cur_vert = seq[i] - 1;
			for (int j = 0; j < n; j++) {
				if (!used[j]) {
					if (cnt == cur_vert) {
						res.push_back(j + 1);
						used[j] = true;
						break;
					}
					else {
						cnt++;
					}
				}
			}
		}

		return res;
	}
};

//update certain gene probability convert to some allel using individual LR and diversity rate 
void UpdateProbability1(double& probability, double iLR, double DR) {
	probability += iLR;
	probability = min(probability, 1.0 - DR);
	probability = max(probability, DR);
}


Chromosome1 Crossover(Individual1& parent1, Individual1& parent2) {
	auto par1 = &parent1;
	auto par2 = &parent2;
	if (Math::GenInt(0, 1)) {
		swap(par1, par2);
	}

	int n = par1->Chr.Size();
	int crossoverPoint = Math::GenInt(0, n);

	Individual1 result(n);
	for (int i = 0; i < crossoverPoint; i++) {
		result.EncodedSeq[i] = par1->EncodedSeq[i];
		result.Chr.Genes[i] = par1->Chr.Genes[i];
	}

	for (int i = crossoverPoint; i < n; i++) {
		result.EncodedSeq[i] = par2->EncodedSeq[i];
		result.Chr.Genes[i] = par2->Chr.Genes[i];
	}


	//Recalc chromosome probabilities using individual LR


	//for (int i = 0; i < n; i++) {
	//	int seqValue = result.EncodedSeq[i] - 1;

	//	double iLR = SolverFluidGenetic::individualLR;
	//	double DR = SolverFluidGenetic::diversityRate;

	//	UpdateProbability(result.Chr.Genes[i].probabilities[seqValue], sqrt(n) * iLR, DR);
	//}


	for (int i = 0; i < n; i++) {
		int seqValue = result.EncodedSeq[i] - 1;
		auto prob = result.Chr.Genes[i].probabilities;

		assert(n >= 2);

		double iLR = SolverFluidGeneticFirst::individualLR;
		double DR = SolverFluidGeneticFirst::diversityRate;

		UpdateProbability1(prob[seqValue], iLR, DR);
		double decreaseValue = -iLR / double(n - 1);
		for (int j = 0; j < prob.size(); j++) {
			if (j != seqValue) {
				UpdateProbability1(prob[j], decreaseValue, DR);
			}
		}

		result.Chr.Genes[i].probabilities = prob;
	}

	return result.Chr;
}


vector<vector<int>> SplitPaths1(vector<int>& path, InputData& input) {
	if (path.empty()) {
		return {};
	}
	path.insert(path.begin(), { 0 });
	int n = static_cast<int>(path.size());

	vector<double> prefLen(n, 0);
	for (int i = 1; i < n; i++) {
		prefLen[i] = prefLen[i - 1] + input.Distance(path[i - 1], path[i]);
	}

	vector<double> dp(n, INF);
	vector<int> prev_dp(n, INF);

	dp[0] = 0;

	for (int v = 0; v + 1 < n; ++v) {
		assert(dp[v] != INF);

		double len = 0;
		for (int u = v + 1; u < n; ++u) {
			// relaxing using edge (v, u)
			len = prefLen[u] - prefLen[v + 1] + input.Distance(0, path[v + 1]) + input.Distance(path[u], 0);
			if (dp[v] + len < dp[u]) {
				dp[u] = dp[v] + len;
				prev_dp[u] = v;
			}
		}
	}

	vector<vector<int>> paths;
	// min() to fix warning
	for (size_t i = n - 1; i != 0; i = prev_dp[min((size_t)n - 1, i)]) {
		paths.push_back({});
		for (int j = (int)i; j > max(-1, prev_dp[i]); --j) {
			paths.back().push_back(path[j]);
		}
		reverse(paths.back().begin(), paths.back().end());
	}

	path.erase(path.begin());


	return paths;
}


bool fitness_comparator1(Individual1& ind1, Individual1& ind2) {
	return ind1.Fitness > ind2.Fitness;
}


void recalcBlueprint(vector<Chromosome1>& population) {
	int sizePopulation = population.size();

	//iterate on every Gene
	for (int itGene = 0; itGene < population[0].Size(); itGene++) {
		int cntAllelesInGene = population[0].Genes[itGene].Size();
		vector<double> sumProbabilities(cntAllelesInGene, 0.0);

		for (int itChr = 0; itChr < sizePopulation; itChr++) {
			for (int itAlelle = 0; itAlelle < cntAllelesInGene; itAlelle++) {
				sumProbabilities[itAlelle] += population[itChr].Genes[itGene].probabilities[itAlelle];
			}
		}

		//take mean
		for (auto& sum : sumProbabilities) {
			sum /= sizePopulation;
		}

		blueprint.Genes[itGene].probabilities = sumProbabilities;
	}
}


void printPopulation(vector<Chromosome1>& chrs, InputData& input) {
	for (auto& chr : chrs) {
		Individual1 ind = Individual1(chr, input);

		cout << "Seq : ";
		for (auto& to : ind.Seq) {
			cout << to << ' ';
		}
		cout << '\n';

		cout << "Fitness : ";
		cout << ind.Fitness << '\n';

		auto path = SplitPaths(ind.Seq, input);

		cout << "Path : " << '\n';
		cout << path.size() << '\n';
		for (auto& p : path) {
			for (auto vert : p) {
				cout << vert << ' ';
			}
			cout << '\n';
		}
	}
}


bool SolverFluidGeneticFirst::Solve(InputData& input, int populationSize, int cntIteration) {
	int N = populationSize;
	int cnt_vertices = input.Size() - 1;

	//inizialization
	vector<Chromosome1> chrs(N, Chromosome1(cnt_vertices));
	blueprint = Chromosome1(cnt_vertices);

	//cout << "-------------------Init Population-----------------------\n";
	//printPopulation(chrs, input);

	//probability depends on individual's fitness
	vector<double> probabilities(N);

	Individual1 answer;

	int timeLimit = 10;

	double startT = clock();
	int numIteration = 0;

	ofstream mean_first("C:/Users/dimas/Desktop/Jupiter/mean_first.txt");
	ofstream best_first("C:/Users/dimas/Desktop/Jupiter/best_first.txt");


	//while (clock() - startT < timeLimit * CLOCKS_PER_SEC) {
	while (numIteration < cntIteration) {
		numIteration++;

		//calc population
		vector<Individual1> population;
		for (int i = 0; i < N; i++) {
			population.push_back(Individual1(chrs[i], input));
		}

		//calc probability for every individuals to particapate in replication
		//using inverse coefficients because we need to minimize individuals fitnesses

		double meanFitness = 0.0;
		int bestFitness = INT_MAX;
		for (int i = 0; i < N; i++) {
			Individual1 ind = population[i];
			meanFitness += 1.0 / double(ind.Fitness);

			bestFitness = min(bestFitness, ind.Fitness);

			if (ind.Fitness < answer.Fitness) {
				answer = ind;
			}
		}

		if (numIteration % 10 == 0) {
			//fout << numIteration << ' ' << bestFitness << '\n';
			//cout << numIteration << ' ' << bestFitness << '\n';
		}
		//cout << numIteration << ' ' << meanFitness << ' ' << bestFitness << '\n';
		mean_first << numIteration << ' ' << meanFitness << '\n';
		best_first << numIteration << ' ' << bestFitness << '\n';
		cout << numIteration << ' ' << meanFitness << '\n';


		for (int i = 0; i < N; i++) {
			probabilities[i] = (1.0 / double(population[i].Fitness)) / meanFitness;
		}

		chrs.clear();
		//make chromosome for new population
		while (chrs.size() < N) {
			int it_par1 = Math::GenRandom(probabilities);
			//probabilities[it_par1] = max(0.01, probabilities[it_par1] - 1.0);

			int it_par2 = Math::GenRandom(probabilities);
			//probabilities[it_par2] = max(0.01, probabilities[it_par2] - 1.0);

			chrs.push_back(Crossover(population[it_par1], population[it_par2]));
		}

		if (enable_blueprint1) {
			recalcBlueprint(chrs);
		}
	}


	mean_first.close();
	best_first.close();

	//cout << "-------------------Finish Population----------------------------\n";
	//printPopulation(chrs, input);
	//cout << "Fluid Genetic First:\n";
	//cout << "Fitness = " << answer.Fitness << '\n';
	//cout << "Iteration Count = " << numIteration << "\n\n";

	return false;
}
