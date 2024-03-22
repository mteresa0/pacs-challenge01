#ifndef PARAM
#define PARAM

#include <fstream>
#include <functional>
#include <limits>


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
        
        const double alpha;
        const double mu;
        const double sigma;
        const double eta;

        const point_type x_0;

        param(const std::string & solver_type_ ,
            const unsigned int k_max_ = 100, 
            const double & tol_r = 1e-6, 
            const double & tol_s = 1e-6, 
            const double & a_= 0.1, 
            const double & mu_= 0.2, 
            const double & sigma_ = 0.2, 
            const double & eta_ = 0.9,
            const point_type & x0_ = {0.0, 0.0}) :
            solver_type(solver_type_),
            k_max(k_max_),
            tol_residual(tol_r),
            tol_step(tol_s),  
            alpha(a_),           
            mu(mu_),
            sigma(sigma_),
            eta(eta_),
            x_0(x0_)
        {}        

    };

    param read_parameters_from_json(const std::string & , const std::string &);
    
}

#endif //PARAM
