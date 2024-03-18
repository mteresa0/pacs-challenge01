#ifndef PARAM
#define PARAM

#include <functional>
#include <limits>


// @todo from array to vector
namespace minimizer{  

    typedef std::vector<double> point_type;
    typedef std::function<double(const point_type &)> fun_type;
    typedef std::function<point_type (const point_type &)> dfun_type; 

    struct param
    {
        const std::string solver_type;
        const unsigned int k_max;
        const double tol_residual;
        const double tol_step;

        point_type x_0;
        double mu;
        double alpha;
        double sigma;
        double eta;

        param(
            const std::string & solver_type_ = "armijo",
            const unsigned int k_max_ = 100, 
            const double & tol_r = 1e-6, 
            const double & tol_s = 1e-6, 
            const double & mu_= 0.2, 
            const double & a_= 0.11, 
            const double & sigma_ = 0.5, 
            const double & eta_ = 0.9,
            const point_type & x0_ = {0.0, 0.0}) :
            solver_type(solver_type_),
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
