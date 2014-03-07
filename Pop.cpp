#include "Pop.h"


//urandom dev/urandom if it exists use it else use create random seed
using namespace std;
// random seed generator
inline unsigned int create_random_seed() {
    unsigned int v = static_cast<unsigned int>(getpid());
    v += ((v << 15) + (v >> 3)) + 0x6ba658b3; // Spread 5-decimal PID over 32-bit number
    v^=(v<<17);
    v^=(v>>13);
    v^=(v<<5);
    v += static_cast<unsigned int>(time(NULL));
    v^=(v<<17);
    v^=(v>>13);
    v^=(v<<5);
    return (v == 0) ? 0x6a27d958 : (v & 0x7FFFFFFF); // return at most a 31-bit seed
}



inline int mod(int a, int b)
/*the implementation of the % operator and fmod computes the remainder of division and is defined to keep the same sign as the dividend.
Ex: 2%5 = 2 but -2%5 = -2 instead of 3.
*/
{
    int r = a % b;
    return (r < 0? r+b : r);
}

void Population::initialize(int nMaxX, int nMaxY, int nOffspring, float fSigma, double dMut,unsigned int seed, int nTransPos)
{
    m_nMaxX = nMaxX;
    m_nMaxY = nMaxY;
    m_nOffspring = nOffspring;
    m_fSigma = fSigma;
    m_dMut = dMut;
    m_nIndividuals = nMaxX * nMaxY;
    m_nAlleleID = m_nIndividuals;
    transIndices(nTransPos);
    setMutCount();
    if (seed==0)
        {
        seed = create_random_seed();
        cout << "Using Generated PRNG Seed: "<< seed << endl;
        }
    m_myrand.seed(seed);

    // Initialize Population
    for(int iii=1; iii<=m_nIndividuals; iii++)
    {
        m_vPop1.push_back(individual(1,iii,iii-1));
        m_vPop2.push_back(individual(0,0,0));
    }

}



void Population::setMutCount()
{
    m_nMutCount = floor(rand_exp(m_myrand, m_dMut));
}


int Population::mutation(int allele)
{
    if (--m_nMutCount > 0)
    {
        return allele;
    }

    else
    {
        setMutCount();
        return m_nAlleleID++;
    }
}

void Population::evolve(int m_nBurnIn, int m_nGenerations)
{
    for(int ggg=0;ggg<m_nBurnIn;++ggg)
    {
        for(int parent=0; parent<m_nIndividuals;parent++)
        {
            step(parent);
        }
        std::swap(m_vPop1,m_vPop2);
    }

    for(int ggg=0;ggg<m_nGenerations;++ggg)
    {
        for(int parent=0; parent<m_nIndividuals;parent++)
        {
            step(parent);
        }
        std::swap(m_vPop1,m_vPop2);
    }
}

int Population::offDispersal(int x, int y)
{
    double a = m_myrand.get_double52() * 2.0 * M_PI;
    double r = rand_exp(m_myrand,m_fSigma);
    int newX = mod(int(floor(r*cos(a)+x+0.5)), m_nMaxX);
    int newY = mod(int(floor(r*sin(a)+y+0.5)), m_nMaxY);
    return newX * m_nMaxY + newY;
}

void Population::step(int parent)
{
    unsigned int &parentHere = m_vPop1[parent].nWeight;
    if (parentHere)
    {
        parentHere=0;
        int nX = parent/m_nMaxX;
        int nY = parent%m_nMaxY;
        for (int off=0; off<m_nOffspring; off++)
        {
            int nNewCell = offDispersal(nX,nY);
            unsigned int nSeedWeight = m_myrand.get_uint32();
            unsigned int nCellWeight = m_vPop2[nNewCell].nWeight;
            int offAllele = mutation(m_vPop1[parent].nAllele);
            if (nSeedWeight > nCellWeight)
            {
                m_vPop2[nNewCell].nAllele=offAllele;
                m_vPop2[nNewCell].nParent_id=parent;
                m_vPop2[nNewCell].nWeight=nSeedWeight;
            }
        }
    }
}

void Population::transIndices(int nTransPos)
//create a vector of the indices in the trasect based on transect postion nTransPos.
{
    int i0 = nTransPos * m_nMaxY + 0; //index of first position in row[nTransPos]
    for(int yyy=0; yyy<m_nMaxY; yyy++)
    {
        m_vtransIndex.push_back(i0+yyy);
    }


}




int Population::minDist(int a, int b)
//calculate the minimum distance between two indices.
{
    int d1 = mod((a-b),m_nMaxX);
    int d2 = mod((b-a),m_nMaxX);
    return min(d1,d2);
}


/*
void Population::samplePop()
{
    vector<int> vDistances(floor(m_nMaxX/2.0));

    for(iii=0; iii<m_nMaxX; iii++)
    {
        for(jjj=1; jjj<(m_nMaxX-iii); jjj++)
        {
            if(m_vPop2[m_vtransIndex[i]]=m_vPop2[m_vtransIndex[i+j]])
                vDistances[minDist(i,i+j)]++
        }
    }
}


*/



//create a vector to store every possible distance size: floor(xmax/2)
//for i=0 to xmax-1:
//  for j=1 to xmax-i:
//      if trans[i] == trans[i+j]
//          k=minDist(i,i+j,xmax)
//          vec[k] += 1
//then iterate through the vector and divide the count by: xmax.







