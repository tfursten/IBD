#include <fstream>

#include "Pop.h"

//urandom dev/urandom if it exists use it else use create random seed
using namespace std;
// random seed generator
inline unsigned int create_random_seed() {
	unsigned int v;
	//ifstream urandom("/dev/urandom", ios::in|ios::binary);
	//if(urandom.good()) {
		//urandom >> v;
	//} else {
	v = static_cast<unsigned int>(getpid());
	v += ((v << 15) + (v >> 3)) + 0x6ba658b3; // Spread 5-decimal PID over 32-bit number
	v^=(v<<17);
	v^=(v>>13);
	v^=(v<<5);
	v += static_cast<unsigned int>(time(NULL));
	v^=(v<<17);
	v^=(v>>13);
	v^=(v<<5);
	//}
    return (v == 0) ? 0x6a27d958 : (v & 0x7FFFFFFF); // return at most a 31-bit seed
}

inline pair<int,int> i2xy(int i, int mx, int my)
{
    assert(0 <= i && i < mx*my);
    return make_pair((i/my),(i%my));
}

inline int xy2i(int x, int y, int mx, int my) {
	assert(0 <= x && x < mx);
	assert(0 <= y && y < my);
	return x*my+y;
}

inline int xy2i(pair<int,int> xy, int mx, int my) {
	return xy2i(xy.first,xy.second,mx,my);
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
    m_nTransPos = nTransPos;
    setMutCount();
    ostringstream out;
    m_nSample = nSample;
    dist.initialize(dist_name);
    out << "Dispersal distribution set to " << dist.getName() << endl;

    //set Random seed
    if (seed==0) {
        seed = create_random_seed();
        out << "Using Generated PRNG Seed: "<< seed << endl;
    }
    m_myrand.seed(seed);

    // Initialize Population: each individual has unique allele
    for(int iii=0; iii<m_nIndividuals; iii++) {
        m_vPop1.emplace_back(1,iii+1,iii);
    }
    m_vPop2.assign(m_nIndividuals, individual(0,0,0));

    cout << out.str();
    pout << out.str();

}

void Population::setMutCount() {
    m_nMutCount = floor(rand_exp(m_myrand, m_dMut));
}


//Each mutational event creates a new allele unlike any other allele currently in the population
//so that identity in state for two or more alleles is always an indication of idenity by descent.
int Population::mutation(int allele)
{
    if (--m_nMutCount > 0)
        return allele;
    setMutCount();
    return m_nAlleleID++;
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

    dout << "#Gen\tSigma2\tKo\tKe\tf\tnIBD" << endl;
    for(int ggg=0;ggg<m_nGenerations;++ggg)
    {
        for(int parent=0; parent<m_nIndividuals;parent++)
            step(parent);
        if (ggg % m_nSample == 0)
            samplePop(ggg);
        std::swap(m_vPop1,m_vPop2);
    }
}

int wrap_around(int x, int w) {
	return ((x % w) + w) % w;
}

int Population::dispersal(int x, int y)
{
    double a = m_myrand.get_double52() * 2.0 * M_PI;
    double r = dist(m_myrand,m_dSigma);
    double dX = floor(r*cos(a)+x+0.5);
    double dY = floor(r*sin(a)+y+0.5);
    int newX = wrap_around(static_cast<int>(dX), m_nMaxX);
    int newY = wrap_around(static_cast<int>(dY), m_nMaxY);

    return xy2i(newX,newY, m_nMaxX,m_nMaxY);
}

void Population::step(int parent)
{
    unsigned int &parentHere = m_vPop1[parent].nWeight;
    if(parentHere == 0)
    	return;
    parentHere=0;
    pair<int,int> xy = i2xy(parent,m_nMaxX,m_nMaxY);
    int nX = xy.first;
    int nY = xy.second;
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

double minEuclideanDist2(int i, int j, int mx, int my) {
	auto xy1 = i2xy(i,mx,my);
	auto xy2 = i2xy(j,mx,my);
	double dx = abs(1.0*(xy1.first - xy2.first));
	double dy = abs(1.0*(xy1.second - xy2.second));
	dx = (dx < mx*0.5) ? dx : mx-dx;
	dy = (dy < my*0.5) ? dy : my-dy;
	return (dx*dx+dy*dy);
}

void Population::samplePop(int gen)
{
    vector<int> vIBD(1+m_nMaxY/2,0);
    vector<int> vN(1+m_nMaxY/2,0);
    typedef map<int,int> mapType;
    mapType alleleMap;
    int szSample = 0;
    double dSigma2 = 0.0;
    double ko = 0.0;
    double ke = 0.0;

    int i0 = m_nTransPos * m_nMaxY;
    for(int i = i0; i < i0+m_nMaxY; ++i) {
    	individual ind = m_vPop2[i];
        if(ind.nWeight == 0)
    		continue;
    	szSample += 1;
    	alleleMap[ind.nAllele] += 1;
    	int p = ind.nParent_id;
    	dSigma2 += minEuclideanDist2(i,p,m_nMaxX,m_nMaxY);

    	for(int j=i; j < i0+m_nMaxY; ++j) {
    		if(m_vPop2[j].nWeight == 0)
    			continue;
    		int k = j-i;
    		k = (k <= m_nMaxY/2) ? k : m_nMaxY-k;
       		if(ind.nAllele == m_vPop2[j].nAllele) {
    			vIBD[k] += 1;
    		}
    		vN[k] += 1;
    	}
    }
    int df = 0;
    BOOST_FOREACH(mapType::value_type & vv, alleleMap){
            df += vv.second*vv.second;
        }
    ko = (double)alleleMap.size();
    double f = df/(double)(szSample*szSample);
    ke = 1.0/f;
    cout << "Gen: " << gen << " Ko: " << ko << " Ke: " << ke << endl;

    dout << gen << "\t" << dSigma2/(2.0*szSample)<< "\t"<<ko<<"\t"<<ke<<"\t"<<f<<"\t";
    for(unsigned int k=0; k<vIBD.size();++k)
        dout << vIBD[k] << "/" << vN[k] << ((k< vIBD.size()-1) ? "\t" : "\n");
    for(int i=0;i<m_nIndividuals;i++)
    {
        if(m_vPop2[i].nWeight == 0)
            gout << -1 << " ";
        else
            gout << m_vPop2[i].nAllele << " ";
    }
    gout << endl;
}


