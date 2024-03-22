#include "solvers.hpp"
#include <cmath>

namespace minimizer 
{
    // @param solverType choose between "gradient", "inverse_decay", "exponential_decay" and "armijo"
    point_type solve(const fun_type &fun, const dfun_type &dfun, const param & p)
    {
        std::cout<< "\nComputing minimum with " << p.solver_type << " method ...\n";
        solverFun solver = choose_solver(p.solver_type);

        point_type x_min = solver(fun, dfun, p);

        return x_min;
    }

    point_type solve(const fun_type &fun, const param & p)
    {
        auto df = [fun](const point_type & x)
        {
            double h = std::numeric_limits<double>::epsilon() * 10000;
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

    point_type fixed_step_solver(const fun_type &fun, const dfun_type &dfun, const param &p)
    {        
        point_type x_old(p.x_0);
        point_type x_new(new_x(dfun, x_old, p.alpha));
        unsigned int k = 0;
        
        bool check = check_tol(fun, p, x_new, x_old);

        while (k<p.k_max && check)
        {
            x_new = new_x(dfun, x_old, p.alpha);
            check = check_tol_residual(fun, p, x_new, x_old);
            x_old = x_new;            
            ++k;
        }

        print_results(fun, x_old, k);
        return x_old;
    }

    point_type inverse_decay_solver(const fun_type &fun, const dfun_type &dfun, const param &p)
    {
        point_type x_old(p.x_0);
        point_type x_new(new_x(dfun, x_old, p.alpha));
        unsigned int k = 0;
        double a_k = p.alpha;
        
        bool check = check_tol(fun, p, x_new, x_old);

        while (k<p.k_max && check)
        {
            x_new = new_x(dfun, x_old, a_k);
            check = check_tol(fun, p, x_new, x_old);
            a_k = p.alpha/(1 + k*p.mu);
            x_old = x_new;            
            ++k;
        }

        print_results(fun, x_old, k);
        return x_old;
    }

    point_type exponential_decay_solver(const fun_type &fun, const dfun_type &dfun, const param &p)
    {
        point_type x_old(p.x_0);
        point_type x_new(new_x(dfun, x_old, p.alpha));
        unsigned int k = 0;
        double a_k = p.alpha;
        
        bool check = check_tol(fun, p, x_new, x_old);

        while (k<p.k_max && check)
        {
            x_new = new_x(dfun, x_old, a_k);
            check = check_tol(fun, p, x_new, x_old);
            a_k = p.alpha * std::exp(-(k*p.mu));
            x_old = x_new;            
            ++k;
        }

        print_results(fun, x_old, k);
        return x_old;
    }

    point_type armijo_solver(const fun_type &fun, const dfun_type &dfun, const param &p)
    {
        point_type x_old(p.x_0);
        point_type x_new = (new_x(dfun, x_old, p.alpha));
        unsigned int k = 0;
        double a_k = p.alpha;
        
        bool check = check_tol(fun, p, x_new, x_old);
        unsigned int i_max(100), i(0);

        while (k<p.k_max && check)
        {
            i = 0;
            while (!armijo_condition(fun, dfun, p, x_old, a_k) && i<i_max)
            {
                ++i;
                a_k /= 2;
            }
            x_new = new_x(dfun, x_old, a_k);
            check = check_tol(fun, p, x_new, x_old);
            x_old = x_new;            
            ++k;
        }

        print_results(fun, x_old, k);
        return x_old;
    }
    
    point_type heavy_ball_solver(const fun_type &fun, const dfun_type &dfun, const param &p)
    {
        point_type x_older(p.x_0), x_new(p.x_0);
        point_type x_old = new_x(dfun, x_older, p.alpha);
        unsigned int k = 1;
        double a_k = p.alpha;
        
        bool check = check_tol(fun, p, x_old, x_older);
        while (k<p.k_max && check)
        {
            x_new = new_x_heavy_ball(dfun, x_old, x_older, a_k, p.eta);
            check = check_tol(fun, p, x_new, x_old);
            x_older = x_old;
            x_old = x_new;
            ++k;
        }

        print_results(fun, x_old, k);
        return x_old;
    }
    
    point_type nesterov_solver(const fun_type &fun, const dfun_type &dfun, const param &p)
    {
        point_type y(p.x_0), x_old(p.x_0);
        point_type x_new = new_x(dfun, p.x_0, p.alpha);
        unsigned int k = 1;
        double a_k = p.alpha;
        
        bool check = check_tol(fun, p, x_old, x_new);
        while (k<p.k_max && check)
        {
            for (std::size_t i = 0; i<y.size(); ++i){
                y[i] =  - p.eta*x_old[i] + (p.eta+1)*x_new[i];
            }
            x_old = x_new;
            x_new = new_x(dfun, y, a_k);
            check = check_tol(fun, p, x_new, x_old);
            ++k;
        }

        print_results(fun, x_old, k);
        return x_old;
    }

    point_type adaptive_hb_solver(const fun_type &fun, const dfun_type &dfun, const param &p)
    {
        point_type x_minus(p.x_0), x(p.x_0);
        std::size_t n = x.size();
        for (std::size_t i = 0; i<n; ++i)
            x_minus[i] += 0.1;
        point_type g_minus = dfun(x_minus);
        point_type g(n), delta_g(n), delta_x(n), nom_mu(n);
        unsigned int k = 1;
        double L, mu, alp, bet, gamma(p.alpha), norm_dx, sqrt_L, sqrt_mu;
        
        bool check = check_tol(fun, p, x_minus, x);
        while (k<p.k_max && check)
        {
            g = dfun(x);
            for (std::size_t i = 0; i<n; ++i)
            {
                delta_g[i] = g[i] - g_minus[i];
                delta_x[i] = x[i] - x_minus[i];
            }
            g_minus = g;
            x_minus = x;
            norm_dx = norm2(delta_x);
            L = gamma * norm2(delta_g) / norm_dx;
            for (std::size_t i = 0; i<n; ++i)
                nom_mu[i] = delta_g[i] - L * delta_x[i];
            mu = norm2(nom_mu) / norm_dx;
            sqrt_L  = std::sqrt(L);
            sqrt_mu = std::sqrt(mu);
            alp = 2/(sqrt_L + sqrt_mu);
            alp *= alp;
            bet = (sqrt_L - sqrt_mu)/(sqrt_L + sqrt_mu);
            bet *= bet;
            
            for (std::size_t i = 0; i<n; ++i)
                x[i] = x[i] - alp*g[i] + bet*delta_x[i];

            check = check_tol(fun, p, x, x_minus);
            ++k;
        }

        print_results(fun, x, k);
        return x;
    }

    /* evaluates the Armijo condition.
    @note (fun(x0) - fun(x0 - alpha0*df(x0)) >= sigma * alpha0 * norm(df(x0))^2)*/
    bool armijo_condition(const fun_type &fun, const dfun_type &dfun, const param & p, const point_type & x0, const double & a0)
    {
        point_type x1 = new_x(dfun, x0, a0);
        double norm_df = norm2(dfun(x0));
        return ( (fun(x0) - fun(x1)) >= p.sigma * a0 * norm_df*norm_df);
    }

    // @param solverType must be "fixed_step_solver", "inverse_decay", "exponential_decay" or "armijo"
    // @returns pointer to function of the chose solver type
    solverFun choose_solver(const std::string & solverType)
    {
        if (solverType=="fixed_step")
            return &fixed_step_solver;
        else if (solverType=="exponential_decay")
            return &exponential_decay_solver;
        else if (solverType=="inverse_decay")
            return &inverse_decay_solver;
        else if (solverType=="armijo")
            return &armijo_solver;
        else if (solverType=="heavy_ball")
            return &heavy_ball_solver;
        else if (solverType=="nesterov")
            return &nesterov_solver;
        else if (solverType=="adaptive_hb")
            return &adaptive_hb_solver;
        else
        {
            std::cerr << "invalid solver name\n";
            return nullptr;
        }
    }

    point_type new_x_heavy_ball(const dfun_type & dfun, const point_type &x_old, const point_type & x_older, const double &a_k, const double & eta)
    {
        point_type x_new = new_x(dfun, x_old, a_k);
        point_type grad_x = dfun(x_old);
        for (std::size_t i = 0; i<x_new.size(); ++i)
        {
            x_new[i] = x_old[i]*(1+eta) - a_k*grad_x[i] - eta*x_older[i];
        }
        return x_new;
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

    bool check_tol(const fun_type & fun, const param & p, const point_type & x1, const point_type & x0)
    {
        return (check_tol_residual(fun, p, x1, x0) && check_tol_step(p, x1, x0));
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

    void print_results(const fun_type &f, const point_type &x_min, const unsigned int &k)
    {
        std::cout << "Minimum found in " << k << " iterations\n";
        std::cout << "min(f(x)) = " << f(x_min) << std::endl;
        return;
    }
    

}