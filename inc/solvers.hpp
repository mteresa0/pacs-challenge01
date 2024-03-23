#ifndef MINSOLVERS
#define MINSOLVERS

#include <iostream>
#include <functional>
#include "param.hpp"

namespace minimizer::solvers
{        
    typedef std::function<point_type(const point_type &, const fun_type &, const grad_type &, const param &)> solverFun;

    point_type solve(const point_type &, const std::pair<fun_type,grad_type> &, const param &);
    point_type solve(const point_type &, const fun_type &, const param &);

    // solvers
    point_type fixed_step_solver        (const point_type &, const fun_type &, const grad_type &, const param &);
    point_type inverse_decay_solver     (const point_type &, const fun_type &, const grad_type &, const param &);
    point_type exponential_decay_solver (const point_type &, const fun_type &, const grad_type &, const param &); 
    point_type armijo_solver            (const point_type &, const fun_type &, const grad_type &, const param &);
    point_type heavy_ball_solver        (const point_type &, const fun_type &, const grad_type &, const param &);
    point_type nesterov_solver          (const point_type &, const fun_type &, const grad_type &, const param &);
    point_type adaptive_hb_solver       (const point_type &, const fun_type &, const grad_type &, const param &);
    point_type adam_solver              (const point_type &, const fun_type &, const grad_type &, const param &);
    
    bool armijo_condition(const fun_type &, const grad_type &,const param &, const point_type &, const double &);
    solverFun choose_solver(const std::string &);

    point_type new_x_heavy_ball(const grad_type & , const point_type &x_old, const point_type &x_older, const double &, const double & );
    point_type new_x(const grad_type &, const point_type &, const double &);
    bool check_tol_step(const param & p, const point_type &, const point_type &);
    bool check_tol_residual(const fun_type &, const param & p, const point_type &, const point_type &);
    bool check_tol(const fun_type &, const param & p, const point_type &, const point_type &);
    
    double norm2(const point_type & );
    double dist(const point_type &, const point_type &);
    double power(const double &, const unsigned int);
    void print_point(const point_type &);
    void print_results(const fun_type &, const point_type &, const unsigned int &);
    
}

#endif