#include "minSolver.hpp"
#include <cmath>

namespace minimazer {

    // @param solverType choose between "gradient", "inverse_decay", "exponential_decay" and "armijo"
    // 
    point_type minSolver::solve(const std::string & solverType)
    {
        point_type x(x_0);
        solverFunPtr solver = choose_solver(solverType);
        unsigned int i = 0;
        double a_k = alpha;
        
        bool check = true;

        while (i<k_max && check)
        {
            check = (this->*solver)(x, a_k, i);
            ++i;
        }
        std::cout << "\n"<< solverType << " solver:\n";
        std::cout << "minimum found in " << i << " itertions :" << std::endl;
        print_point(x);
        std::cout << "min(f(x)) = " << fun(x) << std::endl;
        return x;
    }

    /* evaluates the Armijo condition.
    @note (fun(x0) - fun(x0 - alpha0*df(x0)) >= sigma * alpha0 * norm(df(x0))^2)*/
    bool minSolver::ArmijoCondition(const point_type & x0, const double & a0, 
    const double & sigma)
    {
        point_type x1 = new_x(x0, a0);
        double norm_df = norm2(dfun(x0));
        return (fun(x0) - fun(x1) >= sigma * a0 * norm_df*norm_df);
    }

    // computes the new point and returns true if the step and the residual tollerance are not exceeded
    // @note alpha is constant
    bool minSolver::gradient_solver(point_type & x_old, double & a0, const unsigned int & k)
    {
        point_type x_new(new_x(x_old, a0));
        bool check (check_tol_residual(x_new, x_old) && check_tol_step(x_new, x_old));
        x_old = x_new;
        return check;
    }

    // computes the new point and returns true if the step and the residual tollerance are not exceeded
    // @note alpha_k = alpha_0/(1 + k*mu)
    bool minSolver::inverse_decay_solver(point_type & x_old, double & a0, const unsigned int & k)
    {
        point_type x_new(new_x(x_old, a0/(1 + k*mu)));
        bool check (check_tol_residual(x_new, x_old) && check_tol_step(x_new, x_old));
        x_old = x_new;
        return check;
    }

    // computes the new point and returns true if the step and the residual tollerance are not exceeded
    // @note alpha_k = alpha_0 * exp(-k*mu)
    bool minSolver::exponential_decay_solver(point_type & x_old, double & a0, const unsigned int & k)
    {
        point_type x_new(new_x(x_old, a0*std::exp(-(k*mu))));
        // print_point(x_new);
        bool check (check_tol_residual(x_new, x_old) && check_tol_step(x_new, x_old));
        x_old = x_new;
        return check;
    }

    bool minSolver::armijo_solver(point_type & x_old, double & a, const unsigned int & k)
    {
        unsigned int i(0), i_max(100);
            while (!ArmijoCondition(x_old, a, sigma) && i<i_max)
            {
                ++i;
                a /= 2;
                // std::cout << a << std::endl;
            }
        point_type x_new(new_x(x_old, a));
        bool check(check_tol_residual(x_new, x_old) && check_tol_step(x_new, x_old));
        x_old = x_new;
        return check;
    }

    // @param solverType must be "gradient", "inverse_decay", "exponential_decay" or "armijo"
    // @returns pointer to function of the chose solver type
    minSolver::solverFunPtr minSolver::choose_solver(const std::string & solverType)
    {
        if (solverType=="gradient")
            return &minSolver::gradient_solver;
        else if (solverType=="exponential_decay")
            return &minSolver::exponential_decay_solver;
        else if (solverType=="inverse_decay")
            return &minSolver::inverse_decay_solver;
        else if (solverType=="armijo")
            return &minSolver::armijo_solver;
        else
        {
            std::cerr << "invalid solver name\n";
            return nullptr;
        }
    }

    point_type minSolver::new_x(const point_type &x, const double &a_k)
    {
        // computes x_{k+1} = x_k - alpha_k*df(x_k)
        // where x_k and df(x_k) are R^2
        point_type new_x;
        point_type df_x = dfun(x);
        new_x[0] = x[0] - a_k*df_x[0];
        new_x[1] = x[1] - a_k*df_x[1];
        return new_x;
    }

    bool minSolver::check_tol_step(const point_type & x1, const point_type & x0)
    {
        return (std::abs(dist(x0,x1))>tol_step);
    }

    bool minSolver::check_tol_residual(const point_type & x1, const point_type & x0)
    {
        return (std::abs(fun(x0)-fun(x1))>tol_residual);
    }

    double norm2(const point_type & x)
    {
        return (std::sqrt(x[0]*x[0] + x[1]*x[1]));
    }

    double dist(const point_type & x0, const point_type & x1)
    {
        // returns d = sqrt(|| x0 - x1 ||^2)
        point_type res;
        for (std::size_t i = 0; i<x0.size(); ++i)
        {
            res[i] = (x0[i]-x1[i]);
        }

        return (norm2(res));
    }

    void print_point(const point_type & x)
    {
        for (std::size_t i = 0; i<x.size(); ++i)
            std::cout << "  x("<< i+1 <<") = " << x[i] << std::endl;
        return;
    }

}