#include "Main.h"

typedef std::pair<int,int> xyCoord;


inline int xy2i( xyCoord xy, int nMaxY ) {return xy.first * nMaxY + xy.second;}

inline int mod(int a, int b)
{
    int r = a % b;
    return (r < 0? r+b : r);
}

inline xyCoord disperse(xorshift64& myrand1, int x, int y, double mu, int nMaxX, int nMaxY )
{
    using namespace std;
    double r = rand_exp(myrand1, mu);
    double a = myrand1.get_double53() *2.0*M_PI;
    return xyCoord(
    mod(int(floor(r*cos(a)+0.5+x)), nMaxX),
    mod(int(floor(r*sin(a)+0.5+y)), nMaxY));
}



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


struct Plant
{
    //Plant( int nw=0, int na=0 ) {nWeight = nw; nAllele = na;}
    unsigned int anWeight[2];
    int anAllele[2];
};

struct allelePairs
{
    int a1;
    int a2;
};

Plant* initPop(Plant anInitPop[], int length)
{
    for (int iii=0; iii< length; iii++)
	{
        anInitPop[iii] = Plant();
        anInitPop[iii].anWeight[0] = 1;
        anInitPop[iii].anAllele[0] = iii;

    }
    return anInitPop;
}

unsigned int geometric (xorshift64& myrand1, double fMutp)
    {
        double u = myrand1.get_double52();
        unsigned int k;
        if (fMutp == 1)
        {
            k=1;
        }
        else
        {
            k = log (u) / log (1-fMutp) + 1;
        }
        return k;
    }

inline bool mutate(unsigned int &gmut)
    {
        if (gmut <= 0)
            return true;
        else return false;
    }




int main(int ac, char** av)
{
	namespace po = boost::program_options;
    using namespace std;

    static int nGenerations, nMaxX, nMaxY, nOffspring;
    unsigned int seed;
    static double fMutp;
    static float mu;

    enum Distribution {EXPONENTIAL, GAUSSIAN, TRIANGULAR};
    Distribution disp;
    string dist,infile, outfile;

    //use xorshift64 for PRNG
    xorshift64 myrand;



    try
    {
        po::options_description generic("General Options");
        generic.add_options()
            ("help", "Produce help message")
            ;

        po::options_description config("Configuration");
        config.add_options()
            ("MaxX,X", po::value<int>(&nMaxX)->default_value(100),"Set X dimension")
            ("MaxY,Y", po::value<int>(&nMaxY)->default_value(100),"Set Y dimension")
            ("Generations,g", po::value<int>(&nGenerations)->default_value(10), "Set number of Generations to run")
            ("Offspring,o", po::value<int>(&nOffspring)->default_value(10), "Set number of offspring per individual")
            ("MutP,m", po::value<double>(&fMutp)->default_value(0.0), "Set mutation rate")
            ("Distribution,d", po::value<string>(&dist), "Set Dispersal Distribution")
            ("Mu", po::value<float>(&mu), "Set dispersal mean")
            ("Output_File, out", po::value<string>(&outfile),"Output File Name")
            ("Seed", po::value<unsigned int>(&seed)->default_value(0), "Set PRNG seed, 0 to create random seed")
            ;

        po::options_description hidden("Hidden Options");
        hidden.add_options()
            ("input-file, -i", po::value<string>(&infile), "input file")
            ;

        po::options_description cmdline_options;
        cmdline_options.add(generic).add(config).add(hidden);

        po::options_description config_file_options;
        config_file_options.add(config);

        po::options_description visible("Allowed Options");
        visible.add(generic).add(config);

        po::positional_options_description p;
        p.add("input-file", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(ac,av).options(cmdline_options).positional(p).run(), vm);
        po::notify(vm);


        if (vm.count("help"))
        {
            cout << visible << "\n";
            return 0;
        }


        ifstream ifs(infile.c_str());
        if (!ifs)
        {
            cout << "can not open config file: "<< infile << "\n";
            return 0;
        }
        else
        {
            po::store(parse_config_file(ifs, config_file_options), vm);
            po::notify(vm);
        }




        cout << "X dimension was set to " << vm["MaxX"].as<int>() << ".\n";
        cout << "Y dimension was set to " << vm["MaxY"].as<int>() << ".\n";
        cout << "Run for " << vm["Generation"].as<int>() << " generations.\n";
        cout << "Number of Offspring is set to " << vm["Offspring"].as<int>() << ".\n";
        cout << "Mutation rate set to " << vm["MutP"].as<double>() << ".\n";
        if (vm.count("Distribution"))
            {
            static std::pair<const char*, Distribution> pairs[] = {{"EXPONENTIAL", EXPONENTIAL},{"TRIANGULAR", TRIANGULAR},{"GAUSSIAN",GAUSSIAN}};
            static bool good = false;
            for (int iii=0; iii<sizeof pairs/sizeof pairs[0]; ++iii)
            {
                if (vm["Distribution"].as<string>()==pairs[iii].first)
                    {
                    cout << "The dispersal distribution is " << vm["Distribution"].as<string>() << ".\n";
                    disp = pairs[iii].second;
                    good = true;
                    }
            }
            if (!good)
                {
                throw std::invalid_argument("Error: Distribution not found\n");
                }
            }

        if (vm["Seed"].as<unsigned int>() == 0)
            {
            seed = create_random_seed();
            myrand.seed(seed);
            cout << "Using Generated PRNG Seed: "<< seed << endl;
            }
        else
            {
            cout << "User set PRNG seed to: " << vm["Seed"].as<unsigned int>() << ".\n";
            myrand.seed(vm["Seed"].as<unsigned int>());
            }


    }

    catch(exception& e)
    {
        cout<< e.what() << "\n";
        return 1;
    }


	//Initialize Population
    const int nCells = nMaxX * nMaxY;
    static int nNextAllele = nCells;
    Plant tmp[nCells];
    Plant * vPop;
    vPop = initPop(tmp,nCells);
    bool gen = false;
    unsigned static int gmut = geometric(myrand,fMutp);
    ofstream fout(outfile);


    //Main Generation Loop

    for (int iii=0; iii<nGenerations; iii++)
    {
        for (int jjj=0; jjj<nCells; jjj++)
        {
            unsigned int &cellIsOccupied = vPop[jjj].anWeight[gen];
            if (cellIsOccupied)
            {
                cellIsOccupied = 0;
                for (int lll=0; lll<nOffspring; lll++)
                {
                    int nX = jjj/nMaxY;
                    int nY = mod(jjj,nMaxY);
                    int nNewCell = xy2i(disperse(myrand, nX, nY, mu, nMaxX, nMaxY), nMaxY);
                    unsigned int nSeedWeight = myrand.get_uint64(); //not sure if long is necessary
                    unsigned int &nCellWeight = vPop[nNewCell].anWeight[!gen];
                    gmut --;
                    bool mflag = mutate(gmut);
                    if (nSeedWeight > nCellWeight)
                        {
                        if (!mflag)
                            vPop[nNewCell].anAllele[!gen] = vPop[jjj].anAllele[gen];
                        else
                            {
                            vPop[nNewCell].anAllele[!gen] = nNextAllele;
                            nNextAllele ++;
                            }
                        vPop[nNewCell].anWeight[!gen] = nSeedWeight;
                        }
                    if (mflag)
                        gmut = geometric(myrand, fMutp);

                }
            }
        }

        int &n = nMaxX; //later add option to change the axis for the transect

        int transect[n];
        for (int iii=0; iii<n; iii++)
        {
            xyCoord xy = xyCoord((nMaxY/2), iii);
            int index = xy2i(xy,nMaxY);
            if (vPop[index].anWeight[!gen])
                transect[iii] = vPop[xy2i(xy, nMaxY)].anAllele[!gen];
            else
                transect[iii] = -1;
        }


        //allelePairs anAllelePairs[n-1][n];
        int nIBD[n-1];


        for (int iii=1; iii<n; iii++)
        {
            int count = 0;
            for (int jjj=0; jjj<n; jjj++)
            {
                if (transect[mod(jjj,n)] == transect[mod((jjj+iii),n)])
                    count ++;
            }
            nIBD[iii-1] = count;


                //allelePairs &current = anAllelePairs[iii-1][jjj];
                //current = allelePairs();
                //current.a1 = transect[mod(jjj,n)];
                //current.a2 = transect[mod((jjj+iii),n)];

        }



        if(fout.is_open())
        {

            for (int iii=0; iii<n-1; iii++)
            {
                fout << nIBD[iii];
            }
            fout << endl;
        }
        else
        {
            cout << "Error: Output file could not be opened." << endl;
        }

        /*for (int iii=0; iii<n-1; iii++)
        {
            cout << nIBD[iii] << endl;
        }
        cout << endl;
            for (int jjj=0; jjj<n; jjj++)
            {
                cout << anAllelePairs[iii][jjj].a1 << "  " << anAllelePairs[iii][jjj].a2 << endl;
            }
            cout << endl;
        }
            */


        gen = !gen;
        //cout <<"Gen" << iii << endl;


    }


	return 0;

}
