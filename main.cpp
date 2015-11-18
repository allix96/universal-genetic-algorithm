#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <ctime>
#include <string>
#include <climits>
#include <iomanip>

using std::cout;
using std::cin;
using std::endl;
using std::left;
using std::setw;
using std::vector;
using std::string;
using std::max;
using std::min;
using std::sort;
using std::ofstream;
using std::ios_base;


typedef unsigned int uint;

class Task
{
public:
	double valMin = -5.12;
	double valMax = 5.12;

	double operator() (vector <double> x)
	{
		int n = x.size();
		int A = 10;

		double returnValue = A * n;
		for (int i = 0; i < n; ++i){
			returnValue += x[i] * x[i] - A * cos(2 * M_PI * x[i]);
		}
		return returnValue;		
	}

} task;

double getRandNum (double randMin, double randMax)
{
	double temp = static_cast<double> (rand()) / RAND_MAX;
	temp *= (randMax - randMin);
	temp += randMin;

	return temp;	
}

class Gen
{
public:
	uint gen;

	Gen()
	{
		gen = rand();
	}

	Gen(const Gen &genIn)
	{
		gen = genIn.gen;
	}

	operator int()
	{
		return gen;
	}

	operator double()
	{
		double randMin = task.valMin;
		double randMax = task.valMax;

		double temp = static_cast<double> (gen) / RAND_MAX;
		temp *= (randMax - randMin);
		temp += randMin;
		return temp;
	}

	Gen& operator= (const int right)
	{
		gen = right;
		return *this;
	}

	Gen& operator+= (const int right)
	{
		gen += right;
		return *this;
	}

	Gen& operator^= (const uint right)
	{
		gen ^= right;
		return *this;
	}


};

class Individual
{
public:
	vector<Gen> gens;

	Individual(int sizeIn = 0): gens(sizeIn){}

	Individual(const Individual &individualIn)
	{
		gens = individualIn.gens;
	}

	Individual& operator= (const Individual &right)
	{
		gens = right.gens;
		return *this;
	}

	vector<double> toDouble()
	{
		vector<double> ret(gens.size());
		for (int i = 0; i < ret.size(); ++i){
			ret[i] = gens[i];
		}
		return ret;
	}

	uint size()
	{
		return gens.size();
	}

	void resize(uint newSize)
	{
		gens.resize(newSize);
	}

};

struct IndivCmp
{
	bool operator() (Individual first, Individual second)
	{ 
		return (task(first.toDouble()) < task(second.toDouble()));
	}

} indivCmp;
 

class Population
{
public:
	vector<Individual> individuals;
	uint gensSize;

	Population(int genSizeIn = 0, int indivSizeIn = 0): individuals(indivSizeIn), gensSize(genSizeIn)
	{
		for (int i = 0; i < indivSizeIn; ++i){
			individuals[i].resize(genSizeIn);
		}
	}

	Population(Population &populationIn)
	{
		individuals = populationIn.individuals;
	}

	void crossingIndiv(uint first, uint second)
	{
		uint bitsNum = sizeof(uint) * 8;
		uint crossBitPos = getRandNum(1, bitsNum - 1);

		vector<Individual> newIndiv(2);
		for (int i = 0; i < newIndiv.size(); ++i){
			newIndiv[i].gens.resize(gensSize);
		}

		uint temValue = 0;
		for (int i = 0; i < gensSize; ++i){
			newIndiv[0].gens[i] = (individuals[first].gens[i] >> crossBitPos) << crossBitPos;
			newIndiv[0].gens[i] += (individuals[second].gens[i] << bitsNum - crossBitPos) >> bitsNum - crossBitPos;

			newIndiv[1].gens[i] = (individuals[second].gens[i] >> crossBitPos) << crossBitPos;
			newIndiv[1].gens[i] += (individuals[first].gens[i] << bitsNum - crossBitPos) >> bitsNum - crossBitPos;
		}
		for (int i = 0; i < newIndiv.size(); ++i){
			individuals.push_back(newIndiv[i]);
		}

	}

	void mutationRand(vector<uint> &mutationChance)
	{
		uint bitsNum = sizeof(uint) * 8;
		

		uint ctrlBit = 1;

		for (int i = 0; i < individuals.size(); ++i)
			if (rand() % 100 < mutationChance[0])
			{
				for (int j = 0; j < gensSize; ++j)
					if (rand() % 100 < mutationChance[1])
					{
						for (int k = 0; k < bitsNum; ++k)
							if (rand() % 100 < mutationChance[2])
							{
								individuals[i].gens[j] ^= ctrlBit << k;
							}
					}		
			}
	}

	void sortAscend()
	{
		sort(individuals.begin(), individuals.end(), indivCmp);
	}

	uint size()
	{
		return individuals.size();
	}

	void resize(uint newSize)
	{
		individuals.resize(newSize);
	}
};



int main()
{
	ofstream outFile("out", ios_base::app);
	unsigned int start_time =  clock();


	//Parametrs
	uint POP_SIZE = 400;
	uint INDIV_SIZE = 16;
	uint STEPS_MAX = 400;
	vector<uint> mutationChance = { 20,		//Chance of individual mutations in %
									30,		//Chance of gen mutations in %
									20};	//Chance of bit mutations in %

	//srand( time(0) );
	double answer = rand();



	uint MATCHES_NUM_MAX = 20;

	Population population(INDIV_SIZE, POP_SIZE);

	double taskValueBackup = answer;

	for (uint step = 0, matches = 0; matches < MATCHES_NUM_MAX && step < STEPS_MAX; ++step)
	{
		//--------crossing
		uint crossElem = 0;
		for (int i = 0; i < POP_SIZE - 1; ++i){
			crossElem = getRandNum(i + 1, POP_SIZE - 1);
			population.crossingIndiv(i, crossElem);
		}

		//--------mutation
		population.mutationRand(mutationChance);
		
		//--------sorting
		population.sortAscend();
		
		//--------resize
		population.resize(POP_SIZE);

		answer = task(population.individuals[0].toDouble());

		if(answer == taskValueBackup)
		{
			++matches;
		}
		else
		{
			matches = 0;
			taskValueBackup = answer;
		}
	}

	//cout << "answer: " << answer << endl;

    unsigned int end_time = clock();
    unsigned int search_time = end_time - start_time; 


	uint outSize = 15;
/*
	outFile << setw(outSize) << left << "POP_SIZE" 
			<< setw(outSize) << left << "Dim"
			<< setw(outSize) << left << "mChance_1" 
			<< setw(outSize) << left << "mChance_2" 
			<< setw(outSize) << left << "mChance_3" << endl;
*/	
	outFile << setw(outSize) << left << answer
			<< setw(outSize) << left << POP_SIZE 
			<< setw(outSize) << left << INDIV_SIZE 
			<< setw(outSize) << left << mutationChance[0] 
			<< setw(outSize) << left << mutationChance[1] 
			<< setw(outSize) << left << mutationChance[2]
			<< setw(outSize) << left << float(search_time) / CLOCKS_PER_SEC << endl;

	outFile.close();
	return 0;
}
