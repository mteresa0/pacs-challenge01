#ifndef TEST_FUNCTION
#define TEST_FUNCTION

#include "param.hpp"
#include <map>
#include <string>

namespace minimizer::test_functions{

    double assignment_fun(const point_type & );
    point_type assignment_grad(const point_type &);

    double beale_fun(const point_type & );
    point_type beale_grad(const point_type &);

    extern const std::map<std::string, std::pair<fun_type, dfun_type>> functions;
    
    std::pair<fun_type, dfun_type> get_functions(const std::string &, const bool &);
    // dfun_type get_gradient(const std::string &);
}
#endif 