
#include <iostream>
#include <vector>
#include <limits>
#include "minSolver.hpp"
#include "param.hpp"
// #include "json.hpp"

double f (const minimizer::point_type & x)
{
    // return (x[0]*x[1] + 4*x[0]*x[0]*x[0]*x[0] + x[1]*x[1] + 3*x[0]);
    return ((x[0]+1)*(x[0]+1) + (x[1]+1)*(x[1]+1));
}



minimizer::point_type df (const minimizer::point_type & x)
{
    minimizer::point_type res;
    // res[0] = (x[1] + 16*x[0]*x[0]*x[0] + 3);
    // res[1] = (x[0] + 2*x[1]);
    res[0] = 2*(x[0]+1);
    res[1] = 2*(x[1]+1);
    return res;
}

int main()// (int argc, char ** argv) 
{
    //

    // tollerance settings
    double toll = 1e-6; 

    minimizer::param p(f, df, 1000);

    minimizer::solve(p, "gradient");
    minimizer::solve(p, "exponential_decay");
    minimizer::solve(p, "inverse_decay");
    minimizer::solve(p, "armijo");

    return 0;
}