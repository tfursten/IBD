#include "Main.h"


int main(int ac, char** av)
{
	namespace po = boost::program_options;
    using namespace std;

    static int nGenerations, nMaxX, nMaxY, nOffspring, nBurnIn, nTransPos, nSample;
    unsigned int seed;
    static double dMut;
    static float fSigma;

    enum Distribution {EXPONENTIAL, GAUSSIAN, TRIANGULAR};
    Distribution disp;
    string dist,infile;
    string outfileName;

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
            ("distribution,d", po::value<string>(&dist), "Set Dispersal Distribution")
            ("sigma,s", po::value<float>(&fSigma)->default_value(2.0), "Set dispersal parameter")
            ("burn,b", po::value<int>(&nBurnIn)->default_value(0),"Set Burn-in Period")
            ("sample", po::value<int>(&nSample)->default_value(1),"Sample every n generations after burn-in")
            ("output_file, out", po::value<string>(&outfileName)->default_value(string("data.txt")),"Output File Name")
            ("seed", po::value<unsigned int>(&seed)->default_value(0), "Set PRNG seed, 0 to create random seed")
            ("transect", po::value<int>(&nTransPos)->default_value(-1),"Set position of transect in X axis. Default: ~Center")
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
        << "Collect data every " << nSample << "Generation(s).\n"
        << "Number of Offspring set to " << nOffspring << ".\n"
        << "Mutation rate set to " << dMut<< ".\n";


        if (vm.count("Distribution"))
            {
            static std::pair<const char*, Distribution> pairs[] = {{"EXPONENTIAL", EXPONENTIAL},{"TRIANGULAR", TRIANGULAR},{"GAUSSIAN",GAUSSIAN}};
            static bool good = false;
            for (int iii=0; iii<sizeof pairs/sizeof pairs[0]; ++iii)
            {
                if (vm["Distribution"].as<string>()==pairs[iii].first)
                    {
                    out << "The dispersal distribution is " << vm["Distribution"].as<string>() << ".\n";
                    disp = pairs[iii].second;
                    good = true;
                    }
            }
            if (!good)
                {
                throw std::invalid_argument("Error: Distribution not found\n");
                }
            }
        out << "Dispersal parameter set to " << fSigma << ".\n";

        if (seed)
            out << "User set PRNG seed to: " << seed << ".\n";
        if (nTransPos == -1)
            nTransPos = floor(nMaxX/2);
        out << "Transect position is set to: " << nTransPos << ".\n";
        cout << "Data saved to: " << outfileName << endl;



    }

    catch(exception& e)
    {
        cout<< e.what() << "\n";
        return 1;
    }

    ofstream fout;
    fout.open(outfileName.c_str());
    fout << out.str();
    cout << out.str();
	//Initialize Population

	Population pop(fout);
	pop.initialize(nMaxX,nMaxY,nOffspring,fSigma,1/dMut,seed,nTransPos, nSample);
	//Run Simulation
	pop.evolve(nBurnIn, nGenerations);



	return 0;

}
