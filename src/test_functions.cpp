#include "test_functions.hpp"
#include "solvers.hpp"

namespace minimizer::test_functions{

    double assignment_fun(const point_type & x) {
        return (x[0]*x[1] + 4*x[0]*x[0]*x[0]*x[0] + x[1]*x[1] + 3*x[0]);
    };

    point_type assignment_grad(const point_type & x){
        point_type res(2);
        res[0] = (x[1] + 16*x[0]*x[0]*x[0] + 3);
        res[1] = (x[0] + 2*x[1]);
        return res;
    }

    double beale_fun (const point_type & x_)
    {        
        double x = x_[0]; double y = x_[1];
        return (solvers::power(1.5 - x + x*y, 2) + solvers::power(2.25 - x + x*y*y, 2) + solvers::power(2.625 - x + x*y*y*y, 2));
    }

    point_type beale_grad(const point_type & x_)
    {
        point_type res(2);
        double x = x_[0]; double y = x_[1];
        res[0] = 2*(1.5 - x + x*y)*(-1 + y) + 2*(2.25 - x + x*y*y)*(-1 + y*y) + 2*(2.625 - x + x*y*y*y)*(-1 + y*y*y);
        res[1] = 2*(1.5 - x + x*y)*(x) + 2*(2.25 - x + x*y*y)*(2*x*y) + 2*(2.625 - x + x*y*y*y)*(3*x*y*y);
        return res;
    }

    const std::map<std::string, std::pair<fun_type, dfun_type>> functions = 
    {
        {"assignment", {assignment_fun,     assignment_grad}}, 
        {"beale",      {beale_fun,          beale_grad}}, 

    };

    std::pair<fun_type, dfun_type> get_functions(const std::string & fun_name, const bool & use_analitic_grad)
    {
        auto it = functions.find(fun_name);

        if (use_analitic_grad) {
            if (it!=functions.end())
            {
                return it->second;
            } else {
                std::cerr << "Function name in configuration file is invalid.\n";
                return {nullptr, nullptr};
            }
        } else {
            if (it!=functions.end())
            {
                auto fun = it->second.first;
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
                return {fun, df};
            } else {
                std::cerr << "Function name in configuration file is invalid.\n";
                return {nullptr, nullptr};
            }
        }
    }

}

