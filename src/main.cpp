
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include "solvers.hpp"
#include "param.hpp"
#include "json.hpp"

using namespace minimizer;

double f (const point_type & x)
{
    return (x[0]*x[1] + 4*x[0]*x[0]*x[0]*x[0] + x[1]*x[1] + 3*x[0]);
}

point_type df (const point_type & x)
{
    point_type res(2);
    res[0] = (x[1] + 16*x[0]*x[0]*x[0] + 3);
    res[1] = (x[0] + 2*x[1]);
    return res;
}

int main()// (int argc, char ** argv) 
{
        // double tol = std::numeric_limits<double>::epsilon()*1000;

    param p1 = read_parameters_from_json("parameters.json", "nesterov");
    param p2 = read_parameters_from_json("parameters.json", "heavy_ball");
    param p3 = read_parameters_from_json("parameters.json", "adaptive_hb");
    param p4 = read_parameters_from_json("parameters.json", "armijo");
    
    // analitic derivative
    solve(f, df, p1);
    solve(f, df, p2);
    solve(f, df, p3);
    solve(f, df, p4);

    // finite difference
    solve(f, p1);
    solve(f, p2);
    solve(f, p3);
    solve(f, p4);

    return 0;
}