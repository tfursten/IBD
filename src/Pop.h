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


typedef pair<int,int> xyCoord;
struct individual
{
    unsigned int nWeight;
    int nAllele;
    int nAllele2;
    int nParent_id;
    int nParent_id2;
    xyCoord xy;


    individual(xyCoord x){
    	nWeight = 0;
    	nAllele = 0;
    	nAllele2 = 0;
    	nParent_id = 0;
    	nParent_id2 = 0;
    	xy = x;
    }

    individual(xyCoord x, unsigned int weight, int allele, int parent_id) {
		nWeight = weight;
		nAllele = allele;
		nAllele2 = 0;
		nParent_id = parent_id;
		nParent_id2 = 0;
        xy = x;
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
	int m_nPopSample;
	int m_nTransPos;
	int m_nLenTrans;
	int m_nTransIdx;
	int m_nBoxIdxMin;
	int m_nBoxIdxMax;
	int m_nDistClass;
	int m_nPairs;
	xorshift64 m_myrand;
	std::ofstream & pout;
	std::ofstream & dout;
	std::ofstream & gout;
	std::ofstream & iout;
	std::ofstream & popout;
	std::ofstream & distout;
	std::ofstream & nbout;
	std::ofstream & demout;
	bool verbose;
	Dispersal disp;
	std::vector<individual> m_vPop1;
	std::vector<individual> m_vPop2;
	std::vector<int> m_vtransIndex;
	std::vector<int> m_vtransDist;
	std::vector<int> m_DistCount;
	std::vector<double> m_vAvgIBD;
	int m_nAlleleID;
	void setMutCount();
	void step(int parent);
	int mutate(int allele);
	void sampleIBD(int gen);
	void samplePop();
	void sampleDist();
	void sampleNb();

protected:
    int(Population::*disperse)(int,int);

public:
    
    Population(std::ofstream &p, std::ofstream &d, std::ofstream &g,std::ofstream &i, std::ofstream &pop, std::ofstream &dt,std::ofstream &nb, std::ofstream &dem, bool v): pout(p), dout(d), gout(g), iout(i), popout(pop), distout(dt), nbout(nb), demout(dem), verbose(v), disp(dt) {};
    void initialize(int nMaxX, int nMaxY, int nOffspring, double dSigma,  double dMut, unsigned int seed, int nTransPos, int nSample, int nPopSample, string dist_name, string bound, float param, bool fast, int nclass, int npairs);
	void evolve(int m_nGenerations, int m_nBurnIn);
};

#endif // POP_H_INCLUDED
