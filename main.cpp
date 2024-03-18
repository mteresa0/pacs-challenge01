
#include <iostream>
#include <vector>
#include <limits>
#include "solvers.hpp"
#include "param.hpp"
#include "json.hpp"

double f (const minimizer::point_type & x)
{
    return (x[0]*x[1] + 4*x[0]*x[0]*x[0]*x[0] + x[1]*x[1] + 3*x[0]);
}

minimizer::point_type df (const minimizer::point_type & x)
{
    minimizer::point_type res(2);
    res[0] = (x[1] + 16*x[0]*x[0]*x[0] + 3);
    res[1] = (x[0] + 2*x[1]);
    return res;
}

int main()// (int argc, char ** argv) 
{
    // double tol = std::numeric_limits<double>::epsilon()*1000;
    double tol = 1e-6;

    // minimizer::param p1("gradient", 1000, tol, tol);
    minimizer::param p2("inverse_decay", 1000, tol, tol);
    minimizer::param p3("exponential_decay", 1000, tol, tol);
    minimizer::param p4("armijo", 1000, tol, tol);
    
    // analitic derivative
    // minimizer::solve(f, df, p1);
    minimizer::solve(f, df, p2);
    minimizer::solve(f, df, p3);
    minimizer::solve(f, df, p4);

    // finite difference
    // minimizer::solve(f, p1);
    // minimizer::solve(f, p2);
    // minimizer::solve(f, p3);
    // minimizer::solve(f, p4);

    return 0;
}