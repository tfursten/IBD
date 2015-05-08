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

typedef pair<int,int> xyCoord;

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

void Population::initialize(int nMaxX, int nMaxY, int nOffspring, double dSigma,  double dMut, unsigned int seed, int nTransPos, int nSample, string dist_name, string bound, float param, bool fast)
{
    m_nMaxX = nMaxX;
    m_nMaxY = nMaxY;
    m_sBound = bound;
    m_nOffspring = nOffspring;
    m_dMut = -log(1.0-dMut);
    m_nIndividuals = nMaxX * nMaxY;
    m_nAlleleID = m_nIndividuals;
    m_nTransPos = nTransPos;
    m_nLenTrans = ceil(m_nMaxX/2.0);
    m_nTransIdx = m_nTransPos*m_nMaxX + m_nMaxX/4;
    setMutCount();
    ostringstream out;
    m_nSample = nSample;
    m_dSigma = dSigma;
    disp.initialize(dist_name, m_nMaxX, m_nMaxY, fast, m_sBound, dSigma, param);
    out << "Dispersal distribution set to " << disp.getName() << ".\n" ;
    out << "Extra parameter set to " << param << ".\n";
    out << "Landscape set to " << m_sBound << ".\n";


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
    m_vPop2.assign(m_nIndividuals, individual());

    cout << out.str();
    pout << out.str();

}

void Population::setMutCount() {
    m_nMutCount = floor(rand_exp(m_myrand, m_dMut));
}


//Each mutational event creates a new allele unlike any other allele currently in the population
//so that identity in state for two or more alleles is always an indication of idenity by descent.
int Population::mutate(int allele)
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

    dout << "#Gen\tSigma2\tSigma3\tKo\tKe\tf\tIBD" << endl;
    gout << "#Gen\tSigma2\tSigma3\tKo\tKe\tf\tIBD" << endl;
    iout << "#Gen\tSigma2\tSigma3\tKo\tKe\tf\tIBD" << endl;
    for(int ggg=0;ggg<m_nGenerations;++ggg)
    {
        for(int parent=0; parent<m_nIndividuals;parent++)
            step(parent);
        if (ggg % m_nSample == 0)
            samplePop(ggg);
        std::swap(m_vPop1,m_vPop2);
    }
}


void Population::step(int parent)
{
    //unsigned int &parentHere = m_vPop1[parent].nWeight[0];
    individual &parentHere = m_vPop1[parent];
    if(parentHere.nWeight[0] == 0)
    	return;
    parentHere.nWeight[0] = 0;
    parentHere.nWeight[1] = 0;
    parentHere.threshold = 0;
    pair<int,int> xy = i2xy(parent,m_nMaxX,m_nMaxY);
    int nX = xy.first;
    int nY = xy.second;
    for (int off=0; off<m_nOffspring; off++)
    {
        int nNewCell = disp(m_myrand,nX,nY);
        if (nNewCell == -1)
        {
            mutate(parentHere.nAllele[0]);
            continue;
        }
        individual &parentThere = m_vPop2[nNewCell];
        unsigned int nSeedWeight = m_myrand.get_uint32();
        //I figure that the minimum key value (threshold) will become large
        //quickly and will be less likely to be replaced
        //so I stored it as a variable to avoid having 
        //to calulate the min each time.
        unsigned int nCellWeight = parentThere.threshold;
        int offAllele = mutate(parentHere.nAllele[0]);
        if (nSeedWeight > nCellWeight)
        {
            int r = (parentThere.nWeight[0]<=parentThere.nWeight[1])? 0 : 1;
            parentThere.nAllele[r] = offAllele;
            parentThere.nWeight[r] = nSeedWeight;
            parentThere.nParent_id[r] = parent;
            parentThere.threshold = min(parentThere.nWeight[0],parentThere.nWeight[1]);
            //m_vPop2[nNewCell].nAllele=offAllele;
            //m_vPop2[nNewCell].nParent_id=parent;
            //m_vPop2[nNewCell].nWeight=nSeedWeight;
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

double minEuclideanDist3(int i, int j, int mx, int my){
    auto xy1 = i2xy(i,mx,my);
    auto xy2 = i2xy(j,mx,my);
    double dx = abs(1.0*(xy1.first - xy2.first));
    double dy = abs(1.0*(xy1.second - xy2.second));
    dx = (dx < mx*0.5) ? dx : mx-dx;
    dy = (dy < my*0.5) ? dy : my-dy;
    return (dx*dx*dx + dy*dy*dy);
}

double euclideanDist2(int i, int j, int mx, int my) {
    auto xy1 = i2xy(i,mx,my);
    auto xy2 = i2xy(j,mx,my);
    double dx = abs(1.0*(xy1.first - xy2.first));
    double dy = abs(1.0*(xy1.second - xy2.second));
    return (dx*dx+dy*dy);
}

double euclideanDist3(int i, int j, int mx, int my){
    auto xy1 = i2xy(i,mx,my);
    auto xy2 = i2xy(j,mx,my);
    double dx = abs(1.0*(xy1.first - xy2.first));
    double dy = abs(1.0*(xy1.second - xy2.second));
    return (dx*dx*dx+dy*dy*dy);
}
/*
double minAxialDist2(int i, int j, int mx, int my) {
    auto xy1 = i2xy(i,mx,my);
    auto xy2 = i2xy(j,mx,my);
    double dy = abs(1.0*(xy1.second - xy2.second));
    dy = (dy < my*0.5) ? dy : my-dy;
    return (dy*dy);
}

double axialDist2(int i, int j, int mx, int my) {
    auto xy1 = i2xy(i,mx,my);
    auto xy2 = i2xy(j,mx,my);
    double dy = abs(1.0*(xy1.second - xy2.second));
    return dy*dy;
}
*/

void Population::samplePop(int gen)
{
    vector<int> vIBD(m_nLenTrans,0);
    vector<int> vgIBD(m_nLenTrans,0);
    vector<int> vpIBD(m_nLenTrans,0);
    vector<int> vN(m_nLenTrans,0);
    typedef map<int,int> mapType;
    mapType alleleMap;
    int szSample = 0;
    double dSigma2 = 0.0;
    //double dSigma2_1D = 0.0;
    double dSigma3 = 0.0;
    double ko = 0.0;
    double ke = 0.0;
    //loop through inner 50% of transect
    for(int i = m_nTransIdx; i < m_nTransIdx+m_nLenTrans; ++i) {
    	individual & ind = m_vPop2[i];
        if(ind.nWeight[0] == 0)
    		continue;
    	szSample += 1;
    	alleleMap[ind.nAllele[0]] += 1;
    	int & p = ind.nParent_id[0];

        if (m_sBound == "torus")
        {   
            dSigma2 += minEuclideanDist2(i,p,m_nMaxX,m_nMaxY);
            dSigma3 += minEuclideanDist3(i,p,m_nMaxX,m_nMaxY);
        }
        else
        {
            dSigma2 += euclideanDist2(i,p,m_nMaxX,m_nMaxY);
            dSigma3 += euclideanDist3(i,p,m_nMaxX,m_nMaxY);        
        }
 
        for(int j=i; j < m_nTransIdx+m_nLenTrans; ++j) {
            individual & ind2 = m_vPop2[j];
            if(ind2.nWeight[0] == 0)
                continue;
            int k = j-i;
            //if(m_sBound == "torus")
                //***NOT NECESSARY WHEN USING 50% TRANSECT
                //k = (k <= m_nMaxX/2) ? k : m_nMaxX-k;
            //within individual
            if(k==0){
                //check for presence of competition allele
                if(ind.nWeight[1] == 0)
                    continue;
                // count IIS
                if(ind.nAllele[0]==ind.nAllele[1])
                    vIBD[k] += 1;
                // count grandparental IBD
                if(m_vPop1[p].nParent_id[0] == m_vPop1[ind.nParent_id[1]].nParent_id[0])
                    vgIBD[k] += 1;
                // count parental IBD
                if(p == ind.nParent_id[1])
                    vpIBD[k] += 1;
                vN[k] += 1;
            }
            //between individuals
            else{
                // count IIS
                if(ind.nAllele[0] == ind2.nAllele[0])
                    vIBD[k] += 1;
                // count grandparental IBD
                if(m_vPop1[p].nParent_id[0] == m_vPop1[ind2.nParent_id[0]].nParent_id[0]) 
                    vgIBD[k] += 1;
                // count parental IBD
                if(p == ind2.nParent_id[0])
                    vpIBD[k] += 1;
                vN[k] += 1;
            }
            
        }


    }
    int df = 0;
    BOOST_FOREACH(mapType::value_type & vv, alleleMap){
            df += vv.second*vv.second;
        }
    ko = (double)alleleMap.size();
    double f = df/(double)(szSample*szSample);
    ke = 1.0/f;
    if(verbose)
        cout << "Gen: " << gen << " Ko: " << ko << " Ke: " << ke << endl;

    dout << gen << "\t" << dSigma2/(2.0*szSample) <<"\t" << dSigma3/(2.0*szSample)<<"\t"<<ko<<"\t"<<ke<<"\t"<<f<<"\t";
    gout << gen << "\t" << dSigma2/(2.0*szSample) <<"\t" << dSigma3/(2.0*szSample)<<"\t"<<ko<<"\t"<<ke<<"\t"<<f<<"\t";
    iout << gen << "\t" << dSigma2/(2.0*szSample) <<"\t" << dSigma3/(2.0*szSample)<<"\t"<<ko<<"\t"<<ke<<"\t"<<f<<"\t";
    for(unsigned int k=0; k<vIBD.size();++k){
        dout << vIBD[k] << "/" << vN[k] << ((k < vIBD.size()-1) ? "\t" : "\n");
        gout << vgIBD[k] << "/" << vN[k] << ((k < vgIBD.size()-1) ? "\t" : "\n");
        iout << vpIBD[k] << "/" << vN[k] << ((k < vpIBD.size()-1) ? "\t" : "\n");
    }
    //for(int i=0;i<m_nIndividuals;i++)
    //{
      //  if(m_vPop2[i].nWeight <= 0)
        //    gout << -1 << " ";
       // else
         //   gout << m_vPop2[i].nAllele << " ";
   // }
   // gout << endl;
}


