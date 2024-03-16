#ifndef MINSOLVER
#define MINSOLVER

#include <iostream>
#include <functional>
#include "param.hpp"

namespace minimizer 
{        
    typedef bool(*solverFunPtr)(const param &, point_type &, double &, const unsigned int &); 

    point_type solve(const param &, const std::string &);
    bool armijo_condition(const param &, const point_type &, const double &);
    bool gradient_solver(const param &, point_type &, double &, const unsigned int &);
    bool inverse_decay_solver(const param &, point_type &, double &, const unsigned int &);
    bool exponential_decay_solver(const param &, point_type &, double &, const unsigned int &);
    bool armijo_solver(const param &, point_type &, double &, const unsigned int &);
    solverFunPtr choose_solver(const param &, const std::string &);

    point_type new_x(const param & p, const point_type &, const double &);
    bool check_tol_step(const param & p, const point_type &, const point_type &);
    bool check_tol_residual(const param & p, const point_type &, const point_type &);
    
    double norm2(const point_type & );
    double dist(const point_type &, const point_type &);
    void print_point(const point_type &);
}

#endif //MINSOLVER