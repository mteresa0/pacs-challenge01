
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include "solvers.hpp"
#include "param.hpp"
#include "json.hpp"
#include "test_functions.hpp"

using namespace minimizer;

int main()
{
    std::ifstream ifile("configuration.json");
    using json=nlohmann::json;
    json config;
    ifile >> config;
    ifile.close();

    bool use_analitic_gradient = config["use_analitic_gradient"].get<bool>();
    auto func_name = config["function"].get<std::string>();

    auto funcs  = test_functions::get_functions(func_name, use_analitic_gradient);

    for (const auto & solver_name : config["solvers"].get<std::vector<std::string>>())
    {
        param p = read_parameters_from_json("parameters.json", solver_name);
        solvers::solve(funcs, p);
    }

    return 0;
}