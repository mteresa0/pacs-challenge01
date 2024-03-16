#ifndef MINSOLVER
#define MINSOLVER

#include <iostream>
#include <functional>
#include <limits>

namespace minimazer {
        
    typedef std::array<double, 2> point_type;
    typedef std::function<double(const point_type &)> fun_type;
    typedef std::function<point_type (const point_type &)> dfun_type; 


    struct minSolver
    {
        fun_type fun;
        dfun_type dfun;
        const unsigned int k_max;
        const double tol_residual;
        const double tol_step;

        point_type x_0;
        double mu;
        double alpha;
        double sigma;

        minSolver(const fun_type & f, const dfun_type & df, 
            const unsigned int k_max_ = 100, 
            const double & tol_r = std::numeric_limits<double>::epsilon()*1000, 
            const double & tol_s = std::numeric_limits<double>::epsilon()*1000, 
            const point_type & x0_ = {0.0, 0.0}, const double & mu_= 0.5, 
            const double & a_= 0.15, const double & sigma_ = 0.5) :
            fun(f),
            dfun(df),
            k_max(k_max_),
            tol_residual(tol_r),
            tol_step(tol_s), 
            x_0(x0_),
            mu(mu_),
            alpha(a_),
            sigma(sigma_)
        {}
        
        typedef bool(minSolver::*solverFunPtr)(point_type &, double &, const unsigned int &); 

        point_type solve(const std::string &);
        bool ArmijoCondition(const point_type &, const double &, const double &);
        bool gradient_solver(point_type &, double &, const unsigned int &);
        bool inverse_decay_solver(point_type &, double &, const unsigned int &);
        bool exponential_decay_solver(point_type &, double &, const unsigned int &);
        bool armijo_solver(point_type &, double &, const unsigned int &);
        solverFunPtr choose_solver(const std::string &);
        point_type new_x(const point_type &, const double &);
        bool check_tol_step(const point_type &, const point_type &);
        bool check_tol_residual(const point_type &, const point_type &);

    };
    
    double norm2(const point_type & );
    double dist(const point_type &, const point_type &);
    void print_point(const point_type &);
}

#endif //MINSOLVER