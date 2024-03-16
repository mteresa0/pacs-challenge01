#include "minSolver.hpp"
#include <cmath>

namespace minimizer 
{
    // @param solverType choose between "gradient", "inverse_decay", "exponential_decay" and "armijo"
    // 
    point_type solve(const param & p, const std::string & solverType)
    {
        point_type x(p.x_0);
        solverFunPtr solver = choose_solver(p, solverType);
        unsigned int i = 0;
        double a_k = p.alpha;
        
        bool check = true;

        while (i<p.k_max && check)
        {
            check = (*solver)(p, x, a_k, i);
            ++i;
        }
        std::cout << "\n"<< solverType << " solver:\n";
        std::cout << "minimum found in " << i << " itertions :" << std::endl;
        print_point(x);
        std::cout << "min(f(x)) = " << p.fun(x) << std::endl;
        return x;
    }

    /* evaluates the Armijo condition.
    @note (fun(x0) - fun(x0 - alpha0*df(x0)) >= sigma * alpha0 * norm(df(x0))^2)*/
    bool armijo_condition(const param & p, const point_type & x0, const double & a0)
    {
        point_type x1 = new_x(p, x0, a0);
        double norm_df = norm2(p.dfun(x0));
        return (p.fun(x0) - p.fun(x1) >= p.sigma * a0 * norm_df*norm_df);
    }

    // computes the new point and returns true if the step and the residual tollerance are not exceeded
    // @note alpha is constant
    bool gradient_solver(const param & p, point_type & x_old, double & a0, const unsigned int & k)
    {
        point_type x_new(new_x(p, x_old, a0));
        bool check (check_tol_residual(p, x_new, x_old) && check_tol_step(p, x_new, x_old));
        x_old = x_new;
        return check;
    }

    // computes the new point and returns true if the step and the residual tollerance are not exceeded
    // @note alpha_k = alpha_0/(1 + k*mu)
    bool inverse_decay_solver(const param & p, point_type & x_old, double & a0, const unsigned int & k)
    {
        point_type x_new(new_x(p, x_old, a0/(1 + k*p.mu)));
        bool check (check_tol_residual(p, x_new, x_old) && check_tol_step(p, x_new, x_old));
        x_old = x_new;
        return check;
    }

    // computes the new point and returns true if the step and the residual tollerance are not exceeded
    // @note alpha_k = alpha_0 * exp(-k*mu)
    bool exponential_decay_solver(const param & p, point_type & x_old, double & a0, const unsigned int & k)
    {
        point_type x_new(new_x(p, x_old, a0*std::exp(-(k*p.mu))));
        // print_point(x_new);
        bool check (check_tol_residual(p, x_new, x_old) && check_tol_step(p, x_new, x_old));
        x_old = x_new;
        return check;
    }

    bool armijo_solver(const param & p, point_type & x_old, double & a, const unsigned int & k)
    {
        unsigned int i(0), i_max(100);
            while (!armijo_condition(p, x_old, a) && i<i_max)
            {
                ++i;
                a /= 2;
                // std::cout << a << std::endl;
            }
        point_type x_new(new_x(p, x_old, a));
        bool check(check_tol_residual(p, x_new, x_old) && check_tol_step(p, x_new, x_old));
        x_old = x_new;
        return check;
    }

    // @param solverType must be "gradient", "inverse_decay", "exponential_decay" or "armijo"
    // @returns pointer to function of the chose solver type
    solverFunPtr choose_solver(const param & p, const std::string & solverType)
    {
        if (solverType=="gradient")
            return &gradient_solver;
        else if (solverType=="exponential_decay")
            return &exponential_decay_solver;
        else if (solverType=="inverse_decay")
            return &inverse_decay_solver;
        else if (solverType=="armijo")
            return &armijo_solver;
        else
        {
            std::cerr << "invalid solver name\n";
            return nullptr;
        }
    }

    point_type new_x(const param & p, const point_type &x, const double &a_k)
    {
        // computes x_{k+1} = x_k - alpha_k*df(x_k)
        // where x_k and df(x_k) are R^n
        point_type new_x;
        point_type df_x = p.dfun(x);
        for (std::size_t i = 0; i<x.size(); ++i)
            new_x[i] = x[i] - a_k*df_x[i];
        return new_x;
    }

    bool check_tol_step(const param & p, const point_type & x1, const point_type & x0)
    {
        return (std::abs(dist(x0,x1))>p.tol_step);
    }

    bool check_tol_residual(const param & p, const point_type & x1, const point_type & x0)
    {
        return (std::abs(p.fun(x0)-p.fun(x1))>p.tol_residual);
    }

    double norm2(const point_type & x)
    {
        double sum(0);
        for (auto i : x)
            sum += i*i;
        return (std::sqrt(sum));
    }

    // returns d = sqrt(|| x0 - x1 ||^2)
    double dist(const point_type & x0, const point_type & x1)
    {
        point_type res;
        for (std::size_t i = 0; i<x0.size(); ++i)
        {
            res[i] = (x0[i]-x1[i]);
        }

        return (norm2(res));
    }

    // print point
    void print_point(const point_type & x)
    {
        for (std::size_t i = 0; i<x.size(); ++i)
            std::cout << "  x("<< i <<") = " << x[i] << std::endl;
        return;
    }

}