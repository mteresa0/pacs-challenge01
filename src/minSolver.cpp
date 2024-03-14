#include "minSolver.hpp"
#include <cmath>

mins::point_type mins::minSolver::solveByGradient()
{
    point_type x_n, x_o;
    x_o = x0;
    unsigned long int i = 1;
    double a_k = alpha;
    x_n = new_x(x_o, a_k);
    
    bool check = true;

    while (i<k_max && check)
    {
        x_n = new_x(x_o, alpha/(1+i*mu));
        ++i;
        check = check_tol_residual(x_o, x_n) && check_tol_step(x_o, x_n);
        x_o = x_n;
    }

    std::cout << "minimum found in " << i << " itertions :" << std::endl;
    print_point(x_n);
    return x_n;
}

mins::point_type mins::minSolver::new_x(const mins::point_type & x, const double & a)
{
    // computes x_{k+1} = x_k - alpha*df(x_k)
    // where x_k and df(x_k) are R^2
    point_type new_x;
    std::vector<double> df_x = dfun(x);
    new_x.push_back(x[0] - a*df_x[0]);
    new_x.push_back(x[1] - a*df_x[1]);
    return new_x;
}

bool mins::minSolver::check_tol_step(const mins::point_type & x2, const mins::point_type & x1)
{
    return (std::abs(norm(x1)-norm(x2))>tol_step);
}

bool mins::minSolver::check_tol_residual(const mins::point_type & x2, const point_type & x1)
{
    return (std::abs(fun(x1)-fun(x2))>tol_residual);
}

double mins::norm(const mins::point_type & x)
{
    return (x[0]*x[0] + x[1]*x[1]);
}

void mins::print_point(const mins::point_type & x)
{
    std::cout << "  x = " << x[0] << std::endl;
    std::cout << "  y = " << x[1] << std::endl;
    return;
}