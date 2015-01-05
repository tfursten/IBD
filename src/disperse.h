#ifndef DISPERSE_INCLUDED
#define DISPERSE_INCLUDED

#include <boost/algorithm/string/predicate.hpp>

#include "xorshift64.h"
#include "rexp.h"
#include "rnormal.h"

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

class Dispersal
{
protected:
    typedef double(Dispersal::*fptr)(xorshift64&);

public:
	template<class A>
	bool initialize(A &dist_name, double sigma) {
		static const char name_keys[][16] = {
		    "exponential", "triangular", "normal", "rayleigh", "uniform", "ring",
		};
		static const fptr dist_ops[] = {
            &Dispersal::dist_exponential,
            &Dispersal::dist_triangular,
            &Dispersal::dist_halfNormal,
            &Dispersal::dist_rayleigh,
            &Dispersal::dist_uniform,
            &Dispersal::dist_ring,
        };
        int pos = key_switch(dist_name, name_keys);
        if( pos == -1) {
        	std::cerr << "ERROR: Invalid dispersal distribution" << std::endl;
        	return false;
        }
        name = std::string(name_keys[pos]);
        op = dist_ops[pos];
        param = set_param(name,sigma);
        return true;
	}
	    
    inline std::string getName() { return name; }
    
    inline double operator()(xorshift64& rand, double sigma) {
    	return (this->*op)(rand);
    }
    double set_param(std::string name, double sigma);
	double dist_exponential(xorshift64& rand);
	double dist_triangular(xorshift64& rand);
	double dist_halfNormal(xorshift64& rand);
	double dist_rayleigh(xorshift64& rand);
    double dist_uniform(xorshift64& rand);
    double dist_ring(xorshift64& rand);
	
	Dispersal() {
		initialize("exponential",1);
	}
	
private:
    std::string name;
    fptr op;
    double param;
};

#endif // DISPERSE_INCLUDED
