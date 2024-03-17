#include "solvers.hpp"
#include <cmath>

namespace minimizer 
{
    // @param solverType choose between "gradient", "inverse_decay", "exponential_decay" and "armijo"
    point_type solve(const fun_type &fun, const dfun_type &dfun, const param & p)
    {
        point_type x(p.x_0);
        solverFun solver = choose_solver(p.solver_type);
        unsigned int i = 0;
        double a_k = p.alpha;
        
        bool check = true;

        while (i<p.k_max && check)
        {
            check = solver(fun, dfun, p, x, a_k, i);
            ++i;
        }
        std::cout << "\n"<< p.solver_type << " solver:\n";
        std::cout << "minimum found in " << i << " itertions :" << std::endl;
        // print_point(x);
        std::cout << "min(f(x)) = " << fun(x) << std::endl;
        return x;
    }

    point_type solve(const fun_type &fun, const param & p)
    {
        auto df = [fun](const point_type & x)
        {
            double h = std::numeric_limits<double>::epsilon() * 1000;
            minimizer::point_type df_x(x);
            for (unsigned int i = 0; i<x.size(); ++i)
            {
                minimizer::point_type x_plus(x), x_minus(x);
                x_minus[i] -= h; x_plus[i] += h;
                df_x[i] = (fun(x_plus)-fun(x_minus))/2/h;
            }
            return df_x;
        };
        
        return solve(fun, df, p);
    }

    /* evaluates the Armijo condition.
    @note (fun(x0) - fun(x0 - alpha0*df(x0)) >= sigma * alpha0 * norm(df(x0))^2)*/
    bool armijo_condition(const fun_type &fun, const dfun_type &dfun, const param & p, const point_type & x0, const double & a0)
    {
        point_type x1 = new_x(dfun, x0, a0);
        double norm_df = norm2(dfun(x0));
        return (fun(x0) - fun(x1) >= p.sigma * a0 * norm_df*norm_df);
    }

    // computes the new point and returns true if the step and the residual tollerance are not exceeded
    // @note alpha is constant
    bool gradient_solver(const fun_type &fun, const dfun_type &dfun, const param & p, point_type & x_old, double & a0, const unsigned int & k)
    {
        point_type x_new(new_x(dfun, x_old, a0));
        bool check (check_tol_residual(fun, p, x_new, x_old) && check_tol_step(p, x_new, x_old));
        x_old = x_new;
        return check;
    }

    // computes the new point and returns true if the step and the residual tollerance are not exceeded
    // @note alpha_k = alpha_0/(1 + k*mu)
    bool inverse_decay_solver(const fun_type &fun, const dfun_type &dfun, const param & p, point_type & x_old, double & a0, const unsigned int & k)
    {
        point_type x_new(new_x(dfun, x_old, a0/(1 + k*p.mu)));
        bool check (check_tol_residual(fun, p, x_new, x_old) && check_tol_step(p, x_new, x_old));
        x_old = x_new;
        return check;
    }

    // computes the new point and returns true if the step and the residual tollerance are not exceeded
    // @note alpha_k = alpha_0 * exp(-k*mu)
    bool exponential_decay_solver(const fun_type &fun, const dfun_type &dfun, const param & p, point_type & x_old, double & a0, const unsigned int & k)
    {
        point_type x_new(new_x(dfun, x_old, a0*std::exp(-(k*p.mu))));
        // print_point(x_new);
        bool check (check_tol_residual(fun, p, x_new, x_old) && check_tol_step(p, x_new, x_old));
        x_old = x_new;
        return check;
    }

    bool armijo_solver(const fun_type &fun, const dfun_type &dfun, const param & p, point_type & x_old, double & a, const unsigned int & k)
    {
        unsigned int i(0), i_max(100);
            while (!armijo_condition(fun, dfun, p, x_old, a) && i<i_max)
            {
                ++i;
                a /= 2;
            }
        point_type x_new(new_x(dfun, x_old, a));
        bool check(check_tol_residual(fun, p, x_new, x_old) && check_tol_step(p, x_new, x_old));
        x_old = x_new;
        return check;
    }

    // @param solverType must be "gradient", "inverse_decay", "exponential_decay" or "armijo"
    // @returns pointer to function of the chose solver type
    solverFun choose_solver(const std::string & solverType)
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

    point_type new_x(const dfun_type & dfun, const point_type &x, const double &a_k)
    {
        // computes x_{k+1} = x_k - alpha_k*df(x_k)
        // where x_k and df(x_k) are R^n
        point_type new_x(x);
        point_type df_x = dfun(x);
        for (std::size_t i = 0; i<x.size(); ++i)
            new_x[i] = x[i] - a_k*df_x[i];
        return new_x;
    }

    bool check_tol_step(const param & p, const point_type & x1, const point_type & x0)
    {
        return (std::abs(dist(x0,x1))>p.tol_step);
    }

    bool check_tol_residual(const fun_type & fun, const param & p, const point_type & x1, const point_type & x0)
    {
        return (std::abs(fun(x0)-fun(x1))>p.tol_residual);
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
        point_type res(2);
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
            std::cout << "     x("<< i <<") = " << x[i] << std::endl;
        return;
    }

}