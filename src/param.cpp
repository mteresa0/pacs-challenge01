
#include <fstream>
#include <iostream>
#include "param.hpp"
#include "json.hpp"

namespace minimizer{

    // armijo_param::armijo_param()
    // {
    //     std::ifstream ifile("parameters.json");
    //     using json=nlohmann::json;
    //     json jp;
    //     ifile >> jp;
    //     ifile.close();

    //     double alpha = jp["solvers"]["armijo"]["alpha"].get<double>();
    //     sigma = jp["solvers"]["armijo"]["sigma"].get<double>();
        
    // }

    param read_parameters_from_json(const std::string & filename, const std::string & solver_name)
    {
        std::ifstream ifile(filename);
        using json=nlohmann::json;
        json jp;
        ifile >> jp;
        ifile.close();

        auto tol_step = jp["step_tolerance"].get<double>();
        auto tol_res = jp["residual_tolerance"].get<double>();
        auto k_max = jp["max_iterations"].get<unsigned int>();

        auto nan = std::numeric_limits<double>::quiet_NaN();

        if (solver_name=="armijo")
        {
            auto alpha = jp["solvers"][solver_name]["alpha"].get<double>();
            auto sigma = jp["solvers"][solver_name]["sigma"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, alpha, nan, sigma, nan);
        } else if (solver_name=="inverse_decay")
        {
            auto alpha = jp["solvers"][solver_name]["alpha"].get<double>();
            auto mu = jp["solvers"][solver_name]["mu"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, alpha, mu, nan, nan);
        
        } else if (solver_name=="exponential_decay")
        {
            auto alpha = jp["solvers"][solver_name]["alpha"].get<double>();
            auto mu = jp["solvers"][solver_name]["mu"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, alpha, mu, nan, nan);      
        } else if (solver_name=="fixed_step")
        {
            auto alpha = jp["solvers"][solver_name]["alpha"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, alpha, nan, nan, nan);       
        } else if (solver_name=="heavy_ball")
        {
            auto eta = jp["solvers"][solver_name]["memory"].get<double>();
            auto alpha = jp["solvers"][solver_name]["alpha"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, alpha, nan, nan, eta);
        }else if (solver_name=="nesterov")
        {
            auto eta = jp["solvers"][solver_name]["memory"].get<double>();
            auto alpha = jp["solvers"][solver_name]["alpha"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, alpha, nan, nan, eta);
        }else if (solver_name=="adaptive_hb")
        {
            auto alpha = jp["solvers"][solver_name]["alpha"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, alpha, nan, nan, nan);
        }

        return param("");
    }



}
