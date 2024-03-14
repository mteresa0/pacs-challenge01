
#include <iostream>
#include <vector>
#include "minSolver.hpp"


double f (const mins::point_type & x)
{
    return (x[0]*x[1] + 4*x[0]*x[0]*x[0]*x[0] + x[1]*x[1] + 3*x[0]);
}

std::vector<double> df (const mins::point_type & x)
{
    std::vector<double> res(2);
    res[0] = (x[1] + 16*x[0]*x[0]*x[0] + 3);
    res[1] = (x[0] + 2*x[1]);
    return res;
}

int main()// (int argc, char ** argv) 
{
    //

    // // tollerance settings
    // double eps_r = 1e-6;
    // double eps_s = 1e-6;

    mins::minSolver solver(f,df, 1000);
    mins::point_type x_min = solver.solveByGradient();
    std::cout << "\nf(X_MIN) = " << f(x_min) << std::endl;

    std::cout << f({-0.55416, 0.27708});

    return 0;
}