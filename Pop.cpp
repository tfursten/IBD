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
{
    int r = a % b;
    return (r < 0? r+b : r);
}

void Population::initialize(int nMaxX, int nMaxY, int nOffspring, float fSigma, double dMut,unsigned int seed)
{
    m_nMaxX = nMaxX;
    m_nMaxY = nMaxY;
    m_nOffspring = nOffspring;
    m_fSigma = fSigma;
    m_dMut = dMut;
    m_nIndividuals = nMaxX * nMaxY;
    m_nAlleleID = m_nIndividuals;
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
        pop1.push_back(individual(1,iii,iii-1));
        pop2.push_back(individual(0,0,0));
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
        std::swap(pop1,pop2);
    }

    for(int ggg=0;ggg<m_nGenerations;++ggg)
    {
        for(int parent=0; parent<m_nIndividuals;parent++)
        {
            step(parent);
        }
        std::swap(pop1,pop2);
    }
}

int Population::dispersal(int x, int y)
{
    double a = m_myrand.get_double52() * 2.0 * M_PI;
    double r = rand_exp(m_myrand,m_fSigma);
    int newX = mod(int(floor(r*cos(a)+x+0.5)), m_nMaxX);
    int newY = mod(int(floor(r*sin(a)+y+0.5)), m_nMaxY);
    return newX * m_nMaxY + newY;
}

void Population::step(int parent)
{
    unsigned int &parentHere = pop1[parent].nWeight;
    if (parentHere)
    {
        parentHere=0;
        int nX = parent/m_nMaxX;
        int nY = parent%m_nMaxY;
        for (int off=0; off<m_nOffspring; off++)
        {
            int nNewCell = dispersal(nX,nY);
            unsigned int nSeedWeight = m_myrand.get_uint32();
            unsigned int nCellWeight = pop2[nNewCell].nWeight;
            int offAllele = mutation(pop1[parent].nAllele);
            if (nSeedWeight > nCellWeight)
            {
                pop2[nNewCell].nAllele=offAllele;
                pop2[nNewCell].nParent_id=parent;
                pop2[nNewCell].nWeight=nSeedWeight;
            }
        }
    }
}









