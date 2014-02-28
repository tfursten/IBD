#ifndef POP_H_INCLUDED
#define POP_H_INCLUDED

#include <vector>
#include <iostream>
#include <cmath>
#include <unistd.h>
#include "xorshift64.h"
#include "rexp.h"



struct individual
{
    unsigned int nWeight;
    int nAllele;
    int nParent_id;

    individual(unsigned int weight, int allele, int parent_id)
        {
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
    int m_nOffspring;
    float m_fSigma;
	double m_dMut;
	int m_nMutCount;
	int m_nIndividuals;
	xorshift64 m_myrand;
	std::vector<individual> pop1;
	std::vector<individual> pop2;
	std::vector<int> m_vtransect;
	int m_nGenerations;
	int m_nBurnIn;
	int m_nAlleleID;
	void setMutCount();
	int dispersal(int x, int y);
	void step(int parent);
	int mutation(int allele);
	void transIndices(int nTransPos);
public:
    void initialize(int nMaxX, int nMaxY, int nOffspring, float fSigma, double dMut, unsigned int seed, int nTransPos);
	void evolve(int m_nGenerations, int m_nBurnIn);



};





#endif // POP_H_INCLUDED
