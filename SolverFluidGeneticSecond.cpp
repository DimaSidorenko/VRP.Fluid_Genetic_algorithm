#include "SolverFluidGeneticSecond.h"


bool enable_blueprint = true;

struct Gene {
	int N; // cnt alleles
	vector<double> probabilities; // probabilities for gene to become corresponding allel 

	Gene() {
		N = rand();
		probabilities.resize(N);
		for (size_t i = 0; i < probabilities.size(); i++)
		{
			probabilities[i] = 1 / double(N);
		}
	}
	
	Gene(int size) {
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


struct Chromosome {
	int N; // size of chromosome
	vector<Gene> Genes; // probabilities for sequences

	Chromosome() {
		N = 0;
		Genes = {};
	}

	Chromosome(int N) {
		this->N = N;

		Genes.resize(N);
		for (size_t i = 0; i < N; i++){
			Genes[i] = Gene(N);
		}
	}

	int Size() {
		return N;
	}
};

//it may be last chromosome of popupation or have some flag in population to detect that it is blueprint
Chromosome blueprint;


struct Individual {
	Chromosome Chr; // probabilities for Sequences
	vector<int> Seq; // sequences
	int Fitness;


	Individual() {
		Chr = Chromosome(0);
		Seq = {};
		Fitness = INF;
	}
	
	Individual(int size) {
		Chr = Chromosome(size);
		Seq.resize(size);
		Fitness = INF;
	}

	Individual(Chromosome &chr, InputData& input) {
		Chr = chr;
		Seq = BornAnIndividual(chr);
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

	vector<int> BornAnIndividual(Chromosome &chr) {
		int n = chr.Size();
		vector<int> seq(n);
		vector<int> used(n + 1, false);

		double gLR = SolverFluidGeneticSecond::globalLR;


		for (int i = 0; i < n; i++) {
			if (enable_blueprint) {
				for (int j = 0; j < chr.Genes[i].Size(); j++) {
					chr.Genes[i].probabilities[j] = gLR * blueprint.Genes[i].probabilities[j] + (1 - gLR) * chr.Genes[i].probabilities[j];
				}
			}

			vector<double> prob;
			vector<int> key;

			for (int j = 0; j < chr.Genes[i].Size(); j++) {
				if (!used[j + 1]) {
					key.push_back(j + 1);
					prob.push_back(chr.Genes[i].probabilities[j]);
				}
			}

			int randVal = Math::GenRandom(prob);
			seq[i] = key[randVal];
			used[key[randVal]] = true;
		}

		return seq;
	}

	vector<int> Decode(vector<int> &seq) {
		int n = seq.size();
		vector<bool> used(n, false);
		vector<int> res;

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
void UpdateProbability(double &probability, double iLR, double DR) {
	probability += iLR;
	probability = min(probability, 1.0 - DR);
	probability = max(probability, DR);
}


Chromosome Crossover(Individual& parent1, Individual& parent2) {
	auto par1 = &parent1;
	auto par2 = &parent2;
	if (Math::GenInt(0, 1)) {
		swap(par1, par2);
	}

	int n = par1->Chr.Size();
	int crossoverPoint = Math::GenInt(0, n);

	Individual result(n);
	for (int i = 0; i < crossoverPoint; i++) {
		result.Chr.Genes[i] = par1->Chr.Genes[i];
	}

	for (int i = crossoverPoint; i < n; i++) {
		result.Chr.Genes[i] = par2->Chr.Genes[i];
	}

	return result.Chr;
}

vector<vector<int>> SplitPaths(vector<int>& path, InputData& input) {
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


bool fitness_comparator(Individual& ind1, Individual& ind2) {
	return ind1.Fitness < ind2.Fitness;
}


void recalcBlueprint(vector<Chromosome> &population) {		
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
		for (auto &sum : sumProbabilities) {
			sum /= sizePopulation;
		}

		blueprint.Genes[itGene].probabilities = sumProbabilities;
	}
}


void printPopulation(vector<Chromosome>& chrs, InputData &input) {
	for (auto& chr : chrs) {
		Individual ind = Individual(chr, input);

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


bool SolverFluidGeneticSecond::Solve(InputData& input, int populationSize, int cntIteration) {
	int N = populationSize;
	int cnt_vertices = input.Size() - 1;

	//inizialization
	vector<Chromosome> chrs(N, Chromosome(cnt_vertices));
	blueprint = Chromosome(cnt_vertices);

	double iLR = SolverFluidGeneticSecond::individualLR;
	double DR = SolverFluidGeneticSecond::diversityRate;

	int numIteration = 0;

	ofstream mean_second("mean_second.txt");
	ofstream best_second("best_second.txt");

	Individual best_individual;

	while (numIteration < cntIteration) {
		numIteration++;

		//calc population
		vector<Individual> population;
		for (int i = 0; i < N; i++) {
			population.push_back(Individual(chrs[i], input));
		}

		sort(population.begin(), population.end(), fitness_comparator);

		int n = population[0].Chr.Size();
		for (int itInd = 0; itInd < population.size() / 4; itInd++) {
			
			for (int itGene = 0; itGene < n; itGene++) {

				int seqValue = population[itInd].Seq[itGene] - 1;

				UpdateProbability(population[itInd].Chr.Genes[itGene].probabilities[seqValue], iLR, DR);
				double decreaseValue = -iLR / double(n - 1);
				for (int j = 0; j < n; j++) {
					if (j != seqValue) {
						UpdateProbability(population[itInd].Chr.Genes[itGene].probabilities[j], decreaseValue, DR);
					}
				}
			}
		}

		for (int itInd = population.size() - population.size() / 4; itInd < population.size(); itInd++) {

			for (int itGene = 0; itGene < n; itGene++) {

				int seqValue = population[itInd].Seq[itGene] - 1;

				UpdateProbability(population[itInd].Chr.Genes[itGene].probabilities[seqValue], -iLR, DR);
				double decreaseValue = iLR / double(n - 1);
				for (int j = 0; j < n; j++) {
					if (j != seqValue) {
						UpdateProbability(population[itInd].Chr.Genes[itGene].probabilities[j], decreaseValue, DR);
					}
				}
			}
		}

		//calc probability for every individuals to particapate in replication
		//using inverse coefficients because we need to minimize individuals fitnesses
		double meanFitness = 0.0;
		int bestFitness = INT_MAX;
		for (int i = 0; i < N; i++) {
			Individual ind = population[i];
			meanFitness += 1.0 / double(ind.Fitness);

			bestFitness = min(bestFitness, ind.Fitness);

			if (ind.Fitness < best_individual.Fitness) {
				best_individual = ind;
			}
		}


		if (numIteration % 50 == 0) {
			cout << numIteration << ' ' << meanFitness << ' ' << bestFitness << '\n';
		}
		
		mean_second << numIteration << ' ' << meanFitness << '\n';
		best_second << numIteration << ' ' << bestFitness << '\n';
		
		//cout << numIteration << ' ' << meanFitness << ' ' << bestFitness << '\n';
		//cout << numIteration << ' ' << meanFitness << '\n';


		//probabilities to become a parent
		//depends on individual's fitness
		vector<double> probabilities(N);
		for (int i = 0; i < N; i++) {
			probabilities[i] = (1.0 / double(population[i].Fitness))  / meanFitness;
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

		if (enable_blueprint) {
			recalcBlueprint(chrs);
		}
	}

	//Printing results of algorithm
	cout << "Fluid Genetic Second:\n";
	cout << "Fitness = " << best_individual.Fitness << '\n';
	cout << "Iteration Count = " << numIteration << "\n\n";

	mean_second.close();
	best_second.close();


	return false;
}
