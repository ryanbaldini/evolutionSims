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

	vector<int> sampleWithReplacement(int times, vector<int> &p)
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
	
	double normaldev(double mu, double sig)
	{
		double u, v, x, y, q;
		do
		{
			u = doub();
			v = 1.7156*(doub() - 0.5);
			x = u - 0.449871;
			double v_abs = (v > 0.0) ? v : -v;
			y = v_abs + 0.386595;
			q = x*x + y*(0.19600*y - 0.25472*x);
		}
		while(q > 0.27597 && (q > 0.27846 || v*v > -4.0*log(u)*u*u ));
		return mu + sig*v/u;
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

int sum(vector<bool> x)
{
	return accumulate(x.begin(), x.end(), 0);
}

struct Individual
{
	bool male;
	// int nLoci;	//probably store this elsewhere, once
	//genotype
	vector< vector<bool> > G_D;
	vector< vector<bool> > G_P;
	//phenotype
	double D;
	double P;
	
	void calculatePhenotype(double nLoci, double eps_D, double eps_P, Ran &ran)
	{
		D = (sum(G_D[0]) + sum(G_D[1]))/(2.0*nLoci) + ran.normaldev(0, eps_D);
		D = (D < 0.0) ? 0.0 : D;
		P = 10.0*(sum(G_P[0]) + sum(G_P[1]) - nLoci)/(2.0*nLoci) + ran.normaldev(0, eps_P);
	}
	
	//constructor that initializes genes randomly
	Individual(bool male_, int nLoci, double freq_D, double eps_D, double freq_P, double eps_P, Ran &ran): 
	male(male_), 
	G_D(2, vector<bool>(nLoci, 0)),
	G_P(2, vector<bool>(nLoci, 0))
	{
		//set initial frequencies
		for(int i=0; i<nLoci; i++)
		{
			G_D[0][i] = ran.bern(freq_D);
			G_D[1][i] = ran.bern(freq_D);
			G_P[0][i] = ran.bern(freq_P);
			G_P[1][i] = ran.bern(freq_P);
			
		}
		//get phenoype
		// D = (sum(G_D[0]) + sum(G_D[1]))/(2.0*nLoci) + ran.normaldev(0, eps_D);
		// D = (D < 0.0) ? 0.0 : D;
		// P = 10.0*(sum(G_P[0]) + sum(G_P[1]) - nLoci)/(2.0*nLoci) + ran.normaldev(0, eps_P);
		calculatePhenotype(nLoci, eps_D, eps_P, ran);
	}

	//constructuor that just makes blank values
	Individual(bool male_, int nLoci): 
	male(male_), 
	G_D(2, vector<bool>(nLoci, 0)),
	G_P(2, vector<bool>(nLoci, 0))
	{}
		
};


// struct Genome
// {
// 	int nLoci;
// 	vector< vector<bool> > chromosomes;
//
// 	Genome(int l): nLoci(l), chromosomes(2, vector<bool>(l, 0))
// 	{}
// };
//
struct EvolveSim
{
	//parameters
	int n;	//number diploid individuals
	int l;	//number loci
	double c;	//exponential cost of display trait
	int m;	//number of mates compared by female
	double mu;	//mutation rate
	double eps_D;	//environmental noise on display trait
	double eps_P;	//environmental noise on display trait
	// double s; 	//selection coefficient against deleterious alleles
	// double muGood;	//beneficial mutation rate
	// double muBad;	//deleterious mutation rate
	// double muSelf;
	// bool recessive;
	// bool selfing;

	//state
	Ran ran;	//random number generator
	Binomialdev binom;	//binomial generator for # of mutations per chromosome
	vector<Individual> males;
	vector<Individual> females;
	// vector<double> load;	//mutational load at each time interval
	//vector<double> fitness;	//fitness at present
	vector<double> mean_D;
	vector<double> mean_P;
	//vector<double> alleleFrequencies;
	// vector< vector<bool> > selfingLocus;
	// vector<double> selfingLocusFrequency;
	// int iter;	//simulation iteration

	//constructor
	EvolveSim(int nn, int ll, double cc, int mm, double muu, double eps_D_, double eps_P_, double freq_D, double freq_P): 
	n(nn), l(ll), c(cc), m(mm), mu(muu), eps_D(eps_D_), eps_P(eps_P_), ran(0), binom(ll, muu, 1)
	{
		for(int i=0; i<(n/2); i++)
		{
			males.push_back(Individual(true, l, freq_D, eps_D, freq_P, eps_P, ran));
			females.push_back(Individual(false, l, freq_D, eps_D, freq_P, eps_P, ran));
		}
	}

// 	//some evolve functions...
	void evolve(int t)
	{	
		vector<Individual> newMales(n/2, Individual(true,l));
		vector<Individual> newFemales(n/2, Individual(false,l));
		vector<int> aliveMales(n/2);
		//data to be collected
		mean_D.resize(t, 0.0);
		mean_P.resize(t, 0.0);
		for(int i=0; i<t; i++)
		{
			// iter++;
			//record keeping
			for(int j=0; j<n/2; j++)
			{
				mean_D[i] += males[j].D + females[j].D;
				mean_P[i] += males[j].P + females[j].P;
			}
			mean_D[i] /= n;
			mean_P[i] /= n;
			cout << "mean D: " << mean_D[i] << '\n';
			cout << "mean P: " << mean_P[i] << '\n';
			
			// selection on male survival
			for(int j=0; j<n/2; j++)
			{
				aliveMales[j] = ran.bern(exp(-c*males[j].D));
			}
			
			//produce new males
			for(int j=0; j<n/2; j++)
			{
				//randomly select mother
				Individual& mother = females[ran.int64() % (n/2)];
				//have mother compare m males and select 1
				vector<int> suitorIndices = ran.sampleWithReplacement(m, aliveMales);	//shouldn't be with replacement...
				vector<double> probSelect(m);
				for(int k=0; k<m; k++) probSelect[k] = pow(males[suitorIndices[k]].D, mother.P);	//no need to normalize
				Individual& father = males[suitorIndices[ran.sampleWithReplacement(1, probSelect)[0]]];
				//make the baby
				for(int k=0; k<l; k++)
				{
					//might be faster not to have to call bern so many times, but use binom?
					newMales[j].G_D[0][k] = mother.G_D[ran.bern()][k];
					newMales[j].G_D[1][k] = father.G_D[ran.bern()][k];
					newMales[j].G_P[0][k] = mother.G_P[ran.bern()][k];
					newMales[j].G_P[1][k] = father.G_P[ran.bern()][k];					
				}
				//mutate
				for(int chromosome=0; chromosome<2; chromosome++)
				{
					int nMut = binom.dev();
					if(nMut > 0)
					{
						for(int k=0; k<nMut; k++)
						{
							int ind = ran.int64() % l;
							// newMales[j].G_D[chromosome][ind] = (newMales[j].G_D[chromosome][ind] + 1) % 1;
							newMales[j].G_D[chromosome][ind] = (newMales[j].G_D[chromosome][ind]) ? 0 : 1;
						}
					}
					nMut = binom.dev();
					if(nMut > 0)
					{
						for(int k=0; k<nMut; k++)
						{
							int ind = ran.int64() % l;
							// newMales[j].G_P[chromosome][ind] = (newMales[j].G_P[chromosome][ind] + 1) % 1;
							newMales[j].G_P[chromosome][ind] = (newMales[j].G_P[chromosome][ind]) ? 0 : 1;
						}
					}					
				}
				//calculate phenotype
				newMales[j].calculatePhenotype(l, eps_D, eps_P, ran);
			}
			males = newMales;
			
			//produce new females
			for(int j=0; j<n/2; j++)
			{
				//randomly select mother
				Individual& mother = females[ran.int64() % (n/2)];
				//have mother compare m males and select 1
				vector<int> suitorIndices = ran.sampleWithReplacement(m, aliveMales);	//shouldn't be with replacement...
				vector<double> probSelect(m);
				for(int k=0; k<m; k++) probSelect[k] = pow(males[suitorIndices[k]].D, mother.P);	//no need to normalize
				Individual& father = males[suitorIndices[ran.sampleWithReplacement(1, probSelect)[0]]];
				//make the baby
				for(int k=0; k<l; k++)
				{
					//might be faster not to have to call bern so many times, but use binom?
					newFemales[j].G_D[0][k] = mother.G_D[ran.bern()][k];
					newFemales[j].G_D[1][k] = father.G_D[ran.bern()][k];
					newFemales[j].G_P[0][k] = mother.G_P[ran.bern()][k];
					newFemales[j].G_P[1][k] = father.G_P[ran.bern()][k];					
				}
				//mutate
				for(int chromosome=0; chromosome<2; chromosome++)
				{
					int nMut = binom.dev();
					if(nMut > 0)
					{
						for(int k=0; k<nMut; k++)
						{
							int ind = ran.int64() % l;
							// newFemales[j].G_D[chromosome][ind] = (newFemales[j].G_D[chromosome][ind] + 1) % 1;
							newFemales[j].G_D[chromosome][ind] = (newFemales[j].G_D[chromosome][ind]) ? 0 : 1;
						}
					}
					nMut = binom.dev();
					if(nMut > 0)
					{
						for(int k=0; k<nMut; k++)
						{
							int ind = ran.int64() % l;
							// newFemales[j].G_P[chromosome][ind] = (newFemales[j].G_P[chromosome][ind] + 1) % 1;
							newFemales[j].G_P[chromosome][ind] = (newFemales[j].G_P[chromosome][ind]) ? 0 : 1;
						}
					}					
				}
				//calculate phenotype
				newFemales[j].calculatePhenotype(l, eps_D, eps_P, ran);
			}
			females = newFemales;
		}
	}
};
//
int main()
{
	int t = 10000;	//number of generations
	int n = 5000;
	int l = 100;
	double c = 0.1;
	int m = 20;
	double mu = 1e-5;
	double eps_D = 0.01;
	double eps_P = 0.01;
	double freq_D = 0.1;
	double freq_P = 0.5;
		
	// unsigned long long int seed = 5;	//seed for random generator
	// bool recessive = true;				//is the deleterious allele recessive?
	// bool selfing = false;				//defunct; set to false below
	EvolveSim evolveSim(n, l, c, m, mu, eps_D, eps_P, freq_D, freq_P);

	// cout << evolveSim.males.size() << '\n';
	// cout << sum(evolveSim.males[1].G_D[0]) + sum(evolveSim.males[1].G_D[1]) << '\n';
	// cout << evolveSim.males[1].D << '\n';
	//
	clock_t begin = clock();
	evolveSim.evolve(t);
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "time: " << elapsed_secs << '\n';
	//
	// //write results out to csv
	// ofstream output("output.csv", fstream::out);
	// //make headers
	// output << "load" << ',' << "selfingLocusFrequency" << '\n';
	// for(int i=0; i<t; i++)
	// {
	// 	output << evolveSim.load[i] << ',';
	// 	output << evolveSim.selfingLocusFrequency[i] << '\n';
	// }
	// output.close();
}
