#ifndef DISPERSE_INCLUDED
#define DISPERSE_INCLUDED


#include <boost/algorithm/string/predicate.hpp>
#include "xorshift64.h"
#include "rexp.h"
#include "rnormal.h"


class Dispersal
{
public:
    typedef double(Dispersal::*fptr)(xorshift64&, float);
    fptr getFunction(std::string);
    std::string getName();
	double dist_exponential(xorshift64& rand, float sigma);
	double dist_triangular(xorshift64& rand, float sigma);
	double dist_halfNormal(xorshift64& rand, float sigma);
	double dist_rayleigh(xorshift64& rand, float sigma);

private:
    std::string name;
};

#endif // DISPERSE_INCLUDED
