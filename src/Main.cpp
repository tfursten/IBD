#include "Main.h"


int main(int ac, char** av)
{
	namespace po = boost::program_options;
    using namespace std;

    static int nGenerations, nMaxX, nMaxY, nOffspring, nBurnIn, nTransPos, nSample;
    unsigned int seed;
    static double dMut;
    static float fSigma, param;
    bool f;
    string dist_name, bound, infile, outfileName;
    bool verbose;

    ostringstream out;
    try
    {
        po::options_description generic("General Options");
        generic.add_options()
            ("help", "Produce help message")
            ;

        po::options_description config("Configuration");
        config.add_options()
            ("maxX,x", po::value<int>(&nMaxX)->default_value(100),"Set X dimension")
            ("maxY,y", po::value<int>(&nMaxY)->default_value(100),"Set Y dimension")
            ("generations,g", po::value<int>(&nGenerations)->default_value(10), "Set number of Generations to run after burn-in")
            ("offspring,o", po::value<int>(&nOffspring)->default_value(10), "Set number of offspring per individual")
            ("mut,m", po::value<double>(&dMut)->default_value(0.0), "Set mutation rate")
            ("distribution,d", po::value<string>(&dist_name)->default_value("triangular"), "Set Dispersal Distribution")
            ("sigma,s", po::value<float>(&fSigma)->default_value(2.0), "Set dispersal parameter")
            ("burn,b", po::value<int>(&nBurnIn)->default_value(0),"Set Burn-in Period")
            ("sample,t", po::value<int>(&nSample)->default_value(1),"Sample every n generations after burn-in")
            ("output_file,f", po::value<string>(&outfileName)->default_value(string("data")),"Output File Name")
            ("seed", po::value<unsigned int>(&seed)->default_value(0), "Set PRNG seed, 0 to create random seed")
            ("landscape", po::value<string>(&bound)->default_value(string("torus")),"Set boundary conditions: torus or rectangular")
            ("transect", po::value<int>(&nTransPos)->default_value(0),"Set position of transect in X axis.")
            ("verbose", po::value<bool>(&verbose)->default_value(false),"Print data to screen")
            ("sparam", po::value<float>(&param)->default_value(0),"Extra Parameter for dispersal")
            ("fast", po::value<bool>(&f)->default_value(true),"Use fast dispersal when available")
            ;

        po::options_description hidden("Hidden Options");
        hidden.add_options()
            ("input-file", po::value<string>(&infile), "input file")
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


        if (!infile.empty())
        {
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
        }




        out << "X dimension set to " << nMaxX << ".\n"
        << "Y dimension set to " << nMaxY << ".\n"
        << "Run for " << nGenerations << " generations.\n"
        << "Burn " << nBurnIn << " generation(s).\n"
        << "Collect data every " << nSample << " Generation(s).\n"
        << "Number of Offspring set to " << nOffspring << ".\n"
        << "Mutation rate set to " << dMut<< ".\n";



        assert(fSigma>0);
        assert(nTransPos <= nMaxX);
        out << "Dispersal parameter set to " << fSigma << ".\n";

        if (seed)
            out << "User set PRNG seed to: " << seed << ".\n";
        out << "Transect position is set to: " << nTransPos << ".\n"
        << "Population landscape set to: " << bound << ".\n";

    }

    catch(exception& e)
    {
        cout<< e.what() << "\n";
        return 1;
    }

    string datafile = "IBD"+outfileName+".txt";
    string paramfile = outfileName+"_settings.txt";
    string datafile2 = "gIBD"+outfileName+".txt";
    cout << "IBD Data saved to: " << datafile << endl;
    cout << "gIBD Data saved to: " << datafile2 << endl;
    cout << "Parameters saved to: " << paramfile << endl;
    ofstream pout;
    ofstream dout;
    ofstream gout;
    pout.open(paramfile);
    dout.open(datafile);
    gout.open(datafile2);
    pout << out.str();
    cout << out.str();
	//Initialize Population
    clock_t start = clock();
	Population pop(pout, dout, gout, verbose);
	pop.initialize(nMaxX,nMaxY,nOffspring,fSigma,dMut,seed,nTransPos, nSample, dist_name, bound, param, f);
	//Run Simulation
	pop.evolve(nBurnIn, nGenerations);
	clock_t end = clock();
	float seconds = (float)(end-start)/ CLOCKS_PER_SEC;
    cout << "TIME: " << seconds << endl;



	return 0;

}
