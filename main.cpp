
#include <iostream>
#include <vector>
#include <limits>
#include "minSolver.hpp"
#include "param.hpp"
#include "json.hpp"

double f (const minimizer::point_type & x)
{
    return (x[0]*x[1] + 4*x[0]*x[0]*x[0]*x[0] + x[1]*x[1] + 3*x[0]);
    // return ((x[0]+1)*(x[0]+1) + (x[1]+1)*(x[1]+1));
}



minimizer::point_type df (const minimizer::point_type & x)
{
    minimizer::point_type res(2);
    res[0] = (x[1] + 16*x[0]*x[0]*x[0] + 3);
    res[1] = (x[0] + 2*x[1]);
    // res[0] = 2*(x[0]+1);
    // res[1] = 2*(x[1]+1);
    return res;
}

int main()// (int argc, char ** argv) 
{
    //

    // // tollerance settings
    // double toll = 1e-6; 

    minimizer::param p1("gradient");
    minimizer::param p2("inverse_decay");
    minimizer::param p3("exponential_decay");
    minimizer::param p4("armijo");
    

    minimizer::solve(f, df, p1);
    minimizer::solve(f, df, p2);
    minimizer::solve(f, df, p3);
    minimizer::solve(f, df, p4);

    return 0;
}