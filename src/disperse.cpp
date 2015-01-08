#include <iostream>


#include "disperse.h"

double Dispersal::set_param(std::string name, double sigma)
{
    if (name == "exponential")
        return 1.0/sigma;
    else if (name == "rayleigh")
        return sigma;
    else if (name == "triangular")
        return 2.0*sigma;
    else if (name == "normal")
        return sigma * sqrt(2.0);
    else if (name == "ring")
        return sigma;
    else if (name == "uniform")
        return sigma;
}


double Dispersal::dist_exponential(xorshift64& rand)
{
    return rand_exp(rand, param);
}

double Dispersal::dist_triangular(xorshift64& rand)
{
    //where xmin = 0 and xmax = c
    double u = rand.get_double52();
    return param*sqrt(u);
}

double Dispersal::dist_halfNormal(xorshift64& rand)
{
    return rand_abs_normal(rand, 0.0, param);
}

double Dispersal::dist_rayleigh(xorshift64& rand)
{
    //Return offset
    //double param = sigma;
    //double x = rand_normal(rand,0.0,param);
    //double y = rand_normal(rand,0.0,param);
    //return sqrt(x*x+y*y);
    //return param * sqrt(2.0 * rand_exp(rand));
    return rand_normal(rand,0.0,param);
}

double Dispersal::dist_uniform(xorshift64& rand)
{
    double u = rand.get_uint32();
    return u;
}

double Dispersal::dist_ring(xorshift64& rand)
{
    return param;
}



