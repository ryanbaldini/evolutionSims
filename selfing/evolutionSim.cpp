#include <vector>
#include <iostream>
#include <ctime>
#include <numeric>
#include <cmath>
#include <fstream>

using namespace std;

//I've included NR's Ran and RanQ1 here. The latter is faster, shorter.
//They say RanQ1 is good up to about 10^12 samples
//The model uses O(t*n*l), number of sims times pop size times number of loci
//So, e.g., a big simulation would be 10000*10000*100 = 10^10. 
//So we're probably fine. But if anything gets enormous, switch to the other one. 
//(The speed difference is not big)
struct Ran
{
	// unsigned long long u, v, w;
	unsigned long long v;
	
	// Ran(unsigned long long j): v(4101842887655102017LL), w(1)
	// {
	// 	u = j ^ v; int64();
	// 	v = u; int64();
	// 	w = v; int64();
	// }
	Ran(unsigned long long j): v(4101842887655102017LL)
	{
		v ^= j;
		v = int64();
	}
	// inline unsigned long long int64()
	// {
	// 	u = u * 286293555777941757LL + 7046029254386353087LL;
	// 	v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
	// 	w = 4294957665U*(w & 0xffffffff) + (w >> 32);
	// 	unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
	// 	return (x + v) ^ w;
	// }
	inline unsigned long long int64()
	{
		v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
		return v*2685821657736338717LL;
	}
	inline double doub()
	{
		return 5.42101086242752217E-20 * int64();
	}
	inline vector<double> doub(int n)
	{
		vector<double> x(n);
		for(int i=0; i<n; i++) x[i] = doub();
		return x;
	}
	inline bool bern()
	{
		return (int64() % 2);
	}
	inline bool bern(double p)
	{
		if(doub() < p) return true;
		else return false;
	}
	inline vector<bool> bern(int n)
	{
		vector<bool> x(n);
		for(int i=0; i<n; i++) x[i] = bern();
		return x;
	}
	inline vector<bool> bern(int n, double p)
	{
		vector<bool> x(n);
		for(int i=0; i<n; i++) x[i] = bern(p);
		return x;
	}
	
	vector<int> sampleWithReplacement(int times, vector<double> &p)
	{
		int nElements = p.size();
		vector<double> cumSum(nElements);
		partial_sum(p.begin(), p.end(), cumSum.begin());
		double pSum = cumSum[nElements-1];
		
		vector<int> out(times);
		double x;
		int index;
		for(int i=0; i<times; i++)
		{
			x = pSum*doub();
			index = 0;
			while(x > cumSum[index]) index++;
			out[i] = index;
		}
		return out;
	}

	void sampleWithReplacement(int times, vector<double> &p, vector<double> &cumSumP, vector<int> &out)
	{
		int nElements = p.size();
		double pSum = cumSumP[nElements-1];
		
		double x;
		int index;
		for(int i=0; i<times; i++)
		{
			x = pSum*doub();
			index = 0;
			while(x > cumSumP[index]) index++;
			out[i] = index;
		}
	}
};

double gammln(const double xx)
{
	int j;
	double x, tmp, y, ser;
	static const double cof[14]={57.1562356658629235,-59.5979603554754912, 14.1360979747417471,-0.491913816097620199,.339946499848118887e-4, .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3, -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
	if(xx <= 0) throw("bad arg in gammln");
	y = x = xx;
	tmp = x+5.24218750000000000;
	tmp = (x+0.5)*log(tmp)-tmp;
	ser = 0.999999999999997092;
	for (j=0;j<14;j++) ser += cof[j]/++y;
	return tmp+log(2.5066282746310005*ser/x);
}

struct Binomialdev: Ran {
	double pp, p, pb;
	int n;
	double cdf[64];
	
	Binomialdev(int nn, double ppp, unsigned long long int i): Ran(i), pp(ppp), n(nn) {
		int j;
		pb = p = (pp <= 0.5 ? pp : 1.0-pp);
		cdf[0] = exp(n*log(1-p));
		for(j=1; j<64; j++)
			cdf[j] = cdf[j-1] + exp(gammln(n+1.) - gammln(j+1.) - gammln(n-j+1.) + j*log(p) + (n-j)*log(1.-p));
	}
	
	int dev()
	{
		int j, k, kl, km;
		double y;
		y = doub();
		kl = -1;
		k = 64;
		while(k-kl>1)
		{
			km = (kl+k)/2;
			if(y < cdf[km]) k = km;
			else kl = km;
		}
		if(p != pp) k = n - k;
		return k;
	}
};

struct Genome
{
	int nLoci;
	vector< vector<bool> > chromosomes;
	
	Genome(int l): nLoci(l), chromosomes(2, vector<bool>(l, 0))
	{}
};

struct EvolveSim
{
	//parameters
	int n;	//number diploid individuals
	int l;	//number loci
	double s; 	//selection coefficient against deleterious alleles
	double muGood;	//beneficial mutation rate
	double muBad;	//deleterious mutation rate
	double muSelf;
	bool recessive;
	bool selfing;
	
	//state
	Ran ran;	//random number generator
	Binomialdev binomGood;	//binomial generator for positive mutations
	Binomialdev binomBad;	//binomial generator for positive mutations
	vector<Genome> genePool;	//gene pool of alleles
	vector<double> load;	//mutational load at each time interval
	vector<double> fitness;	//fitness at present
	vector<double> alleleFrequencies;
	vector< vector<bool> > selfingLocus;
	vector<double> selfingLocusFrequency;
	int iter;	//simulation iteration
	
	//constructor
	EvolveSim(int nn, int ll, double ss, double muuGood, double muuBad, double muuSelf, int seed, bool recessive_, bool selfing_): n(nn), l(ll), s(ss), muGood(muuGood), muBad(muuBad), muSelf(muuSelf), recessive(recessive_), selfing(selfing_), ran(seed), binomGood(2*ll, muuGood, seed+1), binomBad(2*ll, muuBad, seed+2), fitness(nn), alleleFrequencies(ll, 0.0), selfingLocus(nn, vector<bool>(2, 0)), iter(0)
	{
		for(int i=0; i<n; i++)
		{ 
			genePool.push_back(Genome(l));
			fitness[i] = 1.0;
		}
	}

	inline void updateAlleleFrequencies()
	{
		double twoN = 2.0*n;
		for(int i=0; i<l; i++)
		{
			alleleFrequencies[i] = 0;
			for(int j=0; j<n; j++)
			{
				alleleFrequencies[i] += (genePool[j].chromosomes[0][i] + genePool[j].chromosomes[1][i]);
			}
			alleleFrequencies[i] /= twoN;
		}
		double f = 0.0;
		for(int i=0; i<n; i++)
		{
			f += selfingLocus[i][0] + selfingLocus[i][1];
		}
		selfingLocusFrequency.push_back(f/twoN);
	}
	
	inline void updateFitness()
	{
		int cost;
		if(recessive)
		{
			for(int i=0; i<n; i++)
			{
				cost = 0;
				for(int j=0; j<l; j++)
				{
					cost += genePool[i].chromosomes[0][j]*genePool[i].chromosomes[1][j];
				}
				fitness[i] = 1.0 - s*cost;
				if(fitness[i] < 1e-8) fitness[i] = 1e-8;
			}
		} else {
			//additive
			for(int i=0; i<n; i++)
			{
				cost = accumulate((genePool[i].chromosomes[0]).begin(), (genePool[i].chromosomes[0]).end(), 0.0);
				cost = accumulate((genePool[i].chromosomes[1]).begin(), (genePool[i].chromosomes[1]).end(), cost);
				fitness[i] = 1.0 - s*cost;
				if(fitness[i] < 1e-8) fitness[i] = 1e-8;
			}
		}
	}
	
	inline void recordLoad()
	{
		load.push_back(1.0 - accumulate(fitness.begin(), fitness.end(), 0.0)/n);
	}
	
	inline void makeGamete(Genome &genome, vector<bool> &newGamete)
	{
		for(int i=0; i<l; i++)
			newGamete[i] = genome.chromosomes[ran.bern()][i];
	}
	
	inline void mutate()
	{
		int mutGood, mutBad, ind1, ind2;
		for(int j=0; j<n; j++)
		{
			//genome mutations
			mutGood = binomGood.dev();
			if(mutGood > 0)
			{
				for(int k=0; k<mutGood; k++)
				{
					ind1 = ran.bern();
					ind2 = ran.int64() % l;
					genePool[j].chromosomes[ind1][ind2] = 0;
				}
			}
			mutBad = binomBad.dev();
			if(mutBad > 0)
			{
				for(int k=0; k<mutBad; k++)
				{
					ind1 = ran.bern();
					ind2 = ran.int64() % l;
					genePool[j].chromosomes[ind1][ind2] = 1;
				}
			}
			//selfing locus mutations
			if(iter > 1000)	//need time for inbreeding depression to accumulate
			{
				selfingLocus[j][0] = (selfingLocus[j][0] ^ ran.bern(muSelf));
				selfingLocus[j][1] = (selfingLocus[j][1] ^ ran.bern(muSelf));
			}
		}
	}
	
	//some evolve functions...
	void evolve(int t)
	{
		vector<Genome> newGenePool(n, Genome(l));
		vector< vector<bool> > newSelfingLocus(n, vector<bool>(2, 0));
		vector<int> parentIndices(2);
		vector<double> cumSumFitness(n);
		int mutGood, mutBad, ind1, ind2;
		for(int i=0; i<t; i++)
		{			
			iter++;
			//recordkeeping
			updateAlleleFrequencies();
			updateFitness();
			recordLoad();
			//reproduce with selection
			partial_sum(fitness.begin(), fitness.end(), cumSumFitness.begin());
			// if(selfing)
			// {
			// 	for(int j=0; j<n; j++)
			// 	{
			// 		ran.sampleWithReplacement(1, fitness, cumSumFitness, parentIndices);
			// 		makeGamete(genePool[parentIndices[0]], newGenePool[j].chromosomes[0]);
			// 		makeGamete(genePool[parentIndices[0]], newGenePool[j].chromosomes[1]);
			// 	}
			// }
			// else
			// {
			// 	for(int j=0; j<n; j++)
			// 	{
			// 		ran.sampleWithReplacement(2, fitness, cumSumFitness, parentIndices);
			// 		makeGamete(genePool[parentIndices[0]], newGenePool[j].chromosomes[0]);
			// 		makeGamete(genePool[parentIndices[1]], newGenePool[j].chromosomes[1]);
			// 	}
			// }
			for(int j=0; j<n; j++)
			{
				//select first parent, make gamete
				ran.sampleWithReplacement(1, fitness, cumSumFitness, parentIndices);
				makeGamete(genePool[parentIndices[0]], newGenePool[j].chromosomes[0]);
				newSelfingLocus[j][0] = selfingLocus[parentIndices[0]][ran.bern()];
				//if first parent is outbreeding, select second parent. make sure second parent is not a selfer.
				//otherwise, sample from same parent again: selfing
				if(selfingLocus[parentIndices[0]][0] + selfingLocus[parentIndices[0]][1] == 0)
				{
					ran.sampleWithReplacement(1, fitness, cumSumFitness, parentIndices);
					while(selfingLocus[parentIndices[0]][0] + selfingLocus[parentIndices[0]][1] > 0)
						ran.sampleWithReplacement(1, fitness, cumSumFitness, parentIndices);
					//this could take forever if selfing is very common (but not fixed)
					//better solution would be to sample from the known outbreeding population
				}
				makeGamete(genePool[parentIndices[0]], newGenePool[j].chromosomes[1]);
				newSelfingLocus[j][1] = selfingLocus[parentIndices[0]][ran.bern()];
			}
			genePool = newGenePool;
			selfingLocus = newSelfingLocus;
			
			//mutate
			//ryan mutates so he can jump higher than baby brutal, hisfoe and friend and mentor
			mutate();
		}
	}
};

int main()
{
	int n = 1000;						//population size
	int l = 1000;						//number of loci
	int t = 10000;						//number of generations to simulate
	double s = 0.1;						//selection coefficient
	double muGood = 1e-6;				//beneficial mutation rate
	double muBad = 1e-5;				//deleterious mutation rate
	double muSelf = 1e-5;				//mutation rate of selfing locus
	unsigned long long int seed = 5;	//seed for random generator
	bool recessive = true;				//is the deleterious allele recessive?
	// bool selfing = false;				//defunct; set to false below
	EvolveSim evolveSim(n, l, s, muGood, muBad, muSelf, seed, recessive, false);	
	
	clock_t begin = clock();
	evolveSim.evolve(t);
	clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "time: " << elapsed_secs << '\n';
			
	//write results out to csv
	ofstream output("output.csv", fstream::out);
	//make headers
	output << "load" << ',' << "selfingLocusFrequency" << '\n';
	for(int i=0; i<t; i++)
	{
		output << evolveSim.load[i] << ',';
		output << evolveSim.selfingLocusFrequency[i] << '\n';
	}
	output.close();
}
