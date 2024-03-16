#ifndef PARAM
#define PARAM

#include <functional>
#include <limits>


namespace minimizer{  

    typedef std::array<double, 2> point_type;
    typedef std::function<double(const point_type &)> fun_type;
    typedef std::function<point_type (const point_type &)> dfun_type; 

    struct param
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
        double eta;

        param(const fun_type & f, const dfun_type & df, 
            const unsigned int k_max_ = 100, 
            const double & tol_r = std::numeric_limits<double>::epsilon()*1000, 
            const double & tol_s = std::numeric_limits<double>::epsilon()*1000, 
            const double & mu_= 0.2, const double & a_= 0.15, 
            const double & sigma_ = 0.5, const double & eta_ = 0.9,
            const point_type & x0_ = {0.0, 0.0}) :
            fun(f),
            dfun(df),
            k_max(k_max_),
            tol_residual(tol_r),
            tol_step(tol_s), 
            x_0(x0_),
            mu(mu_),
            alpha(a_),
            sigma(sigma_),
            eta(eta_)
        {}
        
    };
}

#endif //PARAM
