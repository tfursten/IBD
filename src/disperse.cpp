#include <iostream>

#include "disperse.h"


template<class A, class B, int N>
int key_switch(A &ss, const B (&key)[N])
{
    using boost::algorithm::istarts_with;
    for(int iii=0; iii<N; ++iii)
    {
        if(istarts_with(key[iii],ss))
            return iii;
    }
    return (int) -1;

}


Dispersal::fptr Dispersal::getFunction(std::string dist_name)
{
    static const char name_keys[][16] =
        {
            "exponential", "triangular", "normal", "rayleigh"
        };

    static double (Dispersal::*dist_ops[])(xorshift64&,float) =
        {
            &Dispersal::dist_exponential,
            &Dispersal::dist_triangular,
            &Dispersal::dist_halfNormal,
            &Dispersal::dist_rayleigh
        };

    int pos = key_switch(dist_name, name_keys);
    if(pos == -1)
        std::cout << "ERROR: Invalid dispersal distribution" << std::endl;
    name = std::string(name_keys[pos]);
    return dist_ops[pos];

}

std::string Dispersal::getName()
{
    return name;
}



double Dispersal::dist_exponential(xorshift64& rand, float sigma)
{
    float param = sigma;
    return rand_exp(rand, 1.0/param);
}

double Dispersal::dist_triangular(xorshift64& rand, float sigma)
{
    //where xmin = 0 and xmax = c
    float param = 2*sigma;
    double u = rand.get_double52();
    return sqrt(param*param*u);
}

double Dispersal::dist_halfNormal(xorshift64& rand, float sigma)
{
    double param = (double)sigma * sqrt(2);
    return rand_abs_normal(rand, 0.0, param);
}

double Dispersal::dist_rayleigh(xorshift64& rand, float sigma)
{
    float param = sigma;
    double u = rand.get_double52();
    return param * sqrt(-2.0 * log(u));
}




