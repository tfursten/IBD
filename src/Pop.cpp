#include <fstream>

#include "Pop.h"

//urandom dev/urandom if it exists use it else use create random seed
using namespace std;
// random seed generator
inline unsigned int create_random_seed() {
	unsigned int v;
	ifstream urandom("/dev/urandom", ios::in|ios::binary);
	if(urandom.good()) {
		urandom >> v;
	} else {
		v = static_cast<unsigned int>(getpid());
		v += ((v << 15) + (v >> 3)) + 0x6ba658b3; // Spread 5-decimal PID over 32-bit number
		v^=(v<<17);
		v^=(v>>13);
		v^=(v<<5);
		v += static_cast<unsigned int>(time(NULL));
		v^=(v<<17);
		v^=(v>>13);
		v^=(v<<5);
	}
    return (v == 0) ? 0x6a27d958 : (v & 0x7FFFFFFF); // return at most a 31-bit seed
}

inline pair<int,int> i2xy(int i, int mx, int my)
{
    return make_pair((i/my),(i%my));
}

inline int xy2i(int x, int y, int mx, int my) {
	return x*my+y;
}

inline int xy2i(pair<int,int> xy, int mx, int my) {
	return xy2i(xy.first,xy.second,mx,my);
}

int minDist(int a, int b, int max)
//calculate the minimum distance between two indices.
{
    int d1 = (a-b+max) % max;
    int d2 = (b-a+max) % max;
    return std::min(d1,d2);
}

void Population::initialize(int nMaxX, int nMaxY, int nOffspring, double dSigma,  double dMut,unsigned int seed, int nTransPos, int nSample, string dist_name)
{
    m_nMaxX = nMaxX;
    m_nMaxY = nMaxY;
    m_nOffspring = nOffspring;
    m_dSigma = dSigma;
    m_dMut = -log(1.0-dMut);
    m_nIndividuals = nMaxX * nMaxY;
    m_nAlleleID = m_nIndividuals;
    m_fAvgSig = 0.0;
    setMutCount();
    ostringstream out;
    m_nSample = nSample;
    dist.initialize(dist_name);
    out << "Dispersal distribution set to " << dist.getName() << endl;

    //set Random seed
    if (seed==0)
        {
        seed = create_random_seed();
        out << "Using Generated PRNG Seed: "<< seed << endl;
        }
    m_myrand.seed(seed);

    // Initialize Population: each individual has unique allele
    for(int iii=1; iii<=m_nIndividuals; iii++)
    {
        m_vPop1.push_back(individual(1,iii,iii-1));
        m_vPop2.push_back(individual(0,0,0));
    }

    //Make separate class for transect

    //create a vector of the indices in the trasect based on transect postion nTransPos.
    int i0 = nTransPos * m_nMaxY;
    for(int yyy=0; yyy<m_nMaxY; yyy++)
    {
        m_vtransIndex.push_back(i0+yyy);
    }


    //keep running total of pIBD for each distance bin over every sampled generation
    for(int iii=0; iii<m_nMaxY/2; iii++)
    {
        m_vAvgIBD.push_back(0.0);
    }

    //vector containing the minimum distance for each pair of locations in transect
    for(int iii=0; iii<m_nMaxY;iii++)
    {
        for(int jjj=1; jjj<(m_nMaxY-iii);jjj++)
        {
            m_vtransDist.push_back(minDist(iii,iii+jjj,m_nMaxY));
        }
    }

    //vector containing the total number of samples in each distance bin.
    for(int iii=1; iii<=m_nMaxY/2; iii++)
    {
        m_DistCount.push_back(count(m_vtransDist.begin(), m_vtransDist.end(), iii));
    }

    cout << out.str();
    mout << out.str() << "#" << endl;

}

void Population::setMutCount()
{
    m_nMutCount = floor(rand_exp(m_myrand, m_dMut));
}


//Each mutational event creates a new allele unlike any other allele currently in the population
//so that identity in state for two or more alleles is always an indication of idenity by descent.
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
    //run burn-in period
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
        mout << "pIBD Gen " << ggg << ": ";

        if (ggg % m_nSample == 0)
        {
            samplePop();
        }

        std::swap(m_vPop1,m_vPop2);
    }

    mout << "Average IBD: " ;
    for(int iii=0; iii<m_nMaxY; iii++)
    {
        mout << m_vAvgIBD[iii]/(double)m_nGenerations << "\t";
    }
    mout << "\nAverage Sigma2: " <<m_fAvgSig/(float)m_nGenerations << endl;
    cout <<"Done"<< endl;
}

int wrap_around(int x, int w) {
	return ((x % w) + w) % w;
}

int Population::dispersal(int x, int y)
{
    double a = m_myrand.get_double52() * 2.0 * M_PI;
    double r = dist(m_myrand,m_dSigma);
    int newX = wrap_around(static_cast<int>(floor(r*cos(a)+x+0.5)), m_nMaxX);
    int newY = wrap_around(static_cast<int>(floor(r*sin(a)+y+0.5)), m_nMaxY);
    
    return xy2i(newX,newY, m_nMaxX,m_nMaxY);
}

void Population::step(int parent)
{
    unsigned int &parentHere = m_vPop1[parent].nWeight;
    if (parentHere > 0)
    {
        parentHere=0;
        int nX = parent/m_nMaxX;
        int nY = parent%m_nMaxY;
        for (int off=0; off<m_nOffspring; off++)
        {
            int nNewCell = dispersal(nX,nY);
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

void Population::samplePop()
{
    vector<int> vIBD(m_nMaxY/2);
    int kkk = 0;
    for(int iii=0; iii<m_nMaxY; iii++)
    {
        if(m_vPop2[m_vtransIndex[iii]].nWeight)
        {
            for(int jjj=1; jjj<(m_nMaxY-iii); jjj++)
            {
                if(m_vPop2[m_vtransIndex[iii+jjj]].nWeight)
                {
                    if(m_vPop2[m_vtransIndex[iii]].nAllele==m_vPop2[m_vtransIndex[iii+jjj]].nAllele)
                    {
                        int dist = m_vtransDist[kkk];
                        vIBD[dist-1]++;
                        kkk++;
                    }
                    else kkk++;
                }
                else kkk++;
            }
        }
        else kkk+=(m_nMaxY-iii-1);
    }



    for(int iii=0; iii<m_nMaxY/2;iii++)
    {
        double fpIBD =  (double)vIBD[iii] / (double)m_DistCount[iii];
        mout << fpIBD << "\t";
        m_vAvgIBD[iii] += fpIBD;
    }
    mout << endl;



    int nSumDist2 = 0;
    for(int iii=0; iii<m_nMaxY; iii++)
    {
        pair<int,int> npXY = i2xy(m_vPop2[m_vtransIndex[iii]].nParent_id, m_nMaxX, m_nMaxY);
        pair<int,int> niXY = i2xy(m_vtransIndex[iii], m_nMaxX, m_nMaxY);
        nSumDist2 += (pow(minDist(npXY.first,niXY.first,m_nMaxX),2) + pow(minDist(npXY.second,niXY.second,m_nMaxY),2));
    }
    float nSigma2 = (nSumDist2/(float)m_nMaxY) * 0.5;
    mout << "sigma2: " << nSigma2 << endl << endl;
    m_fAvgSig += nSigma2;

}

