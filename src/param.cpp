
#include <fstream>
#include <iostream>
#include "param.hpp"
#include "json.hpp"

namespace minimizer{

    // reads solvers parameters from the parameters.json file
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

        // i really do not like it, it would be better if i could use the inheritance property of classes
        if (solver_name=="armijo")
        {
            auto alpha = jp["solvers"][solver_name]["alpha"].get<double>();
            auto sigma = jp["solvers"][solver_name]["sigma"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, alpha, nan, sigma, nan);
        } 
        else if (solver_name=="inverse_decay")
        {
            auto alpha = jp["solvers"][solver_name]["alpha"].get<double>();
            auto mu = jp["solvers"][solver_name]["mu"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, alpha, mu, nan, nan);
        
        } 
        else if (solver_name=="exponential_decay")
        {
            auto alpha = jp["solvers"][solver_name]["alpha"].get<double>();
            auto mu = jp["solvers"][solver_name]["mu"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, alpha, mu, nan, nan);      
        } 
        else if (solver_name=="fixed_step")
        {
            auto alpha = jp["solvers"][solver_name]["alpha"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, alpha, nan, nan, nan);       
        } 
        else if (solver_name=="heavy_ball")
        {
            auto alpha = jp["solvers"][solver_name]["alpha"].get<double>();
            auto eta   = jp["solvers"][solver_name]["memory"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, alpha, nan, nan, eta);
        }
        else if (solver_name=="nesterov")
        {
            auto alpha = jp["solvers"][solver_name]["alpha"].get<double>();
            auto eta   = jp["solvers"][solver_name]["memory"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, alpha, nan, nan, eta);
        }
        else if (solver_name=="adaptive_hb")
        {
            auto gamma = jp["solvers"][solver_name]["gamma"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, gamma, nan, nan, nan);
        }
        else if (solver_name=="adam")
        {
            auto alpha = jp["solvers"][solver_name]["alpha"].get<double>();
            auto beta1 = jp["solvers"][solver_name]["beta1"].get<double>();
            auto beta2 = jp["solvers"][solver_name]["beta2"].get<double>();
            auto eps   = jp["solvers"][solver_name]["eps"].get<double>();
            return param(solver_name, k_max, tol_res, tol_step, alpha, beta1, beta2, eps);
        }

        return param("");
    }



}
