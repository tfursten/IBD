#ifndef POP_H_INCLUDED
#define POP_H_INCLUDED

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <unistd.h>
#include <map>
#include <boost/foreach.hpp>
#include "xorshift64.h"
#include "rexp.h"
#include "disperse.h"
#include "disk.h"



struct individual
{
    unsigned int nWeight;
    int nAllele;
    int nParent_id;

    individual(unsigned int weight, int allele, int parent_id) {
		nWeight = weight;
        nAllele = allele;
        nParent_id = parent_id;
    }
};


class Population
{
private:
    int m_nMaxX;
    int m_nMaxY;
    std::string m_sBound;
    int m_nOffspring;
    double m_dSigma;
	double m_dMut;
	int m_nMutCount;
	int m_nIndividuals;
	int m_nSample;
	int m_nTransPos;
	xorshift64 m_myrand;
	Dispersal disp;
	std::ofstream & pout;
	std::ofstream & dout;
	std::ofstream & gout;
	bool verbose;
	bool m_bTorus;
	std::vector<individual> m_vPop1;
	std::vector<individual> m_vPop2;
	std::vector<int> m_vtransIndex;
	std::vector<int> m_vtransDist;
	std::vector<int> m_DistCount;
	std::vector<double> m_vAvgIBD;
	int m_nAlleleID;
	void setMutCount();
	void step(int parent);
	int mutation(int allele);
	void samplePop(int gen);

protected:
    int(Population::*disperse)(int,int);

public:
    Population(std::ofstream &p, std::ofstream &d, std::ofstream &g, bool v): pout(p), dout(d), gout(g), verbose(v) {};
    void initialize(int nMaxX, int nMaxY, int nOffspring, double dSigma,  double dMut, unsigned int seed, int nTransPos, int nSample, string dist_name, string bound, float param, bool fast);
	void evolve(int m_nGenerations, int m_nBurnIn);
};

#endif // POP_H_INCLUDED
