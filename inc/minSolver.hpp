#ifndef MINSOLVER
#define MINSOLVER

#include <iostream>
#include <functional>
#include <limits>

namespace mins {

    typedef std::vector<double> point_type;
    typedef std::function<double(const point_type &)> fun_type;
    typedef std::function<point_type (const point_type &)> dfun_type;
    

    double norm(const point_type & );
    void print_point(const point_type &);

    struct minSolver
    {
        fun_type fun;
        dfun_type dfun;
        const unsigned int k_max;
        const double tol_residual;
        const double tol_step;

        point_type x0;
        double mu;
        double alpha;

        minSolver(const fun_type & f, const dfun_type & df, 
            const unsigned int k_max_ = 100, const double & tol_r = std::numeric_limits<double>::epsilon()*10, 
            const double & tol_s = std::numeric_limits<double>::epsilon()*10, const double & mu_= 0.05, 
            const double & a_= 0.11) :
            fun(f),
            dfun(df),
            k_max(k_max_),
            tol_residual(tol_r),
            tol_step(tol_s), 
            x0({0.0, 0.0}),
            mu(mu_),
            alpha(a_)
        {}

        point_type solveByGradient();
        point_type new_x(const point_type &, const double &);
        bool check_tol_step(const mins::point_type & x0, const mins::point_type & x1);
        bool check_tol_residual(const mins::point_type & x0, const mins::point_type & x1);

    };

}

#endif //MINSOLVER