# CHALLENGE 01 - PACS
First challenge of the course Advaced Programming for Scientific Computing, 23/24, Politecnico di Milano.

## Installation
To clone the repository run:

```bash
git clone git@github.com:mteresa0/pacs-challenge01.git
```

## Usage
To compile the main executable and to run the program:
```bash
make run
```
To just run the program:
```bash
./main
```

## Configuration files
To change the settings of the minimizer please change `configuration.json` and  `parameters.json` files accorgingly to the following instructions.

In `configuration.json` you can set multiple types of solvers, with which solve your minimization problem, the function you want to solve, the method the gradient is computed and the initial guess. An example:
```json
{
    "solvers": [
        "armijo",
        "nesterov"
    ],
    "function": "rosenbrock",
    "use_analitic_gradient": true,
    "initial_guess": [0, 0]
}
```
The solvers name and the function name you can use are listed and explained in the subsections below.


All the parameters for each solver are set in `parameters.json`.
Also, in `parameters.json` the tollerances are set as follows:
```json
    "step_tolerance": 1e-6,
    "residual_tolerance": 1e-6,
    "max_iterations": 1000
```

### Solvers
Here are listed all the implemented solver names you can add in the `configuration.json` file:
- `fixed_step` -
    This is the simplest solver: it solves the minimization problem through the gradient method with a fixed step, $\alpha_0$. 
    $$ x_{k+1} = x_{k} - \alpha_0 \nabla f(x_k) $$

    The parameter of `fixed_step` is set in `parameters.json`:
    ```json
    "solvers":{
        "fixed_step": {
            "alpha": 0.005
        }
    }
    ```
    where `alpha` is $\alpha_0$, the fixed step.

- `inverse_decay` -
    This solver finds the minimum using the gradient method with a variable step, $\alpha_k$:
    <!-- $$ x_{k+1} = x_{k} - \alpha_k \nabla f(x_k) $$ -->
    $$ \alpha_k = \frac{\alpha_0}{1 + k\mu} $$

    The parameters of `inverse_decay` are set in `parameters.json`:
    ```json
    "solvers":{
        "inverse_decay": {
            "alpha": 0.005,
            "mu": 0.1
        }
    }
    ```
    where `alpha` is $\alpha_0$, the initial step size, and `mu` is $\mu$.

- `exponential_decay` - 
    This solver finds the minimum using the gradient method with a variable step, $\alpha_k$:
        <!-- $$ x_{k+1} = x_{k} - \alpha_k \nabla f(x_k) $$ -->
        $$ \alpha_k = \alpha_0 e^{- k\mu} $$
    The parameterss of `exponential_decay` are set in `parameters.json`:
    ```json
    "solvers":{
        "exponential_decay": {
            "alpha": 0.005,
            "mu": 0.2
        }
    }
    ```
    where `alpha` is $\alpha_0$, the initial step size, and `mu` is $\mu$.

- `armijo` -
    This solver utilizes the gradient method with a step size that satisfies the Armijo rule. Given a step size, $\alpha$, and a parameter $\sigma \in (0, 0.5)$, if 
    $$ f(x_{k}) - f(x_k - \alpha \nabla f({x_k})) \ge \sigma \alpha ||\nabla f(x_k)||^2 $$
    is not satisfied the step size is halved, $\alpha = \frac{\alpha}{2}$.

    The parameters of `armijo` are set in `parameters.json`:
    ```json
    "solvers":{
        "armijo": {
            "alpha": 0.05,
            "sigma": 0.25
        }
    }
    ```
    where `alpha` is $\alpha$, the initial step size, and `sigma` is $\sigma$, the Armijo parameter.

- `heavy_ball` - 
    This solver utilizes the heavy-ball method: it solves the minimization problem by taking the following update rule:
    $$ x_{k+1} = x_{k} - \alpha \nabla f(x_k) + \eta (x_k - x_{k-1})$$ 
    where $\alpha$ is the step size (in this implementation the step size is fixed) and $\eta$ is the memory parameter.

    The parameters of `heavy_ball` are set in `parameters.json`:
    ```json
    "solvers":{
        "heavy_ball": {
            "alpha": 0.001,
            "memory": 0.875
        }
    }
    ```
    where `alpha` is $\alpha$, the fixed step size, and `memory` is $\eta$, the memory parameter.

- `nesterov` -
    This solver uses a similar strategy to the heavy-ball method. It updates the coordinate by taking the following rule:
    $$ y = x_k + \eta (x_k - x_{k-1}) $$
    $$ x_{k+1} = y - \alpha \nabla f(y)$$ 

    where `alpha` is $\alpha$, the fixed step size, and `memory` is $\eta$, the memory parameter.

    The parameters of `nesterov` are set in `parameters.json`:
    ```json
    "solvers":{
        "nesterov": {
            "alpha": 0.001,
            "memory": 0.9
        }
    }
    ```
    where `alpha` is $\alpha$, the fixed step size, and `memory` is $\eta$, the memory parameter.

- `adaptive_hb` - 
    This solver utilizes an adaptive Heavy-Ball method that estimates the Polyak’s optimal hyper-parameters at each iteration. The theory and implementation follwed of this strategy are outlined in this [paper](https://link.springer.com/content/pdf/10.1007/s10994-022-06215-7.pdf).

    The parameter of `adaptive_hb` is set in `parameters.json`:
    ```json
    "solvers":{
        "adaptive_hb": {
            "gamma": 1.8
        }
    }
    ```
    where `gamma` is the only hyper-parameter to set.


- `adam` -
    This solver utilizes the Adaptive Moment Estimation (Adam) method. It utilizes estimates of the first and second moments of the gradients to dynamically adjust the learning rate for each parameter during the optimization process. The theory and the implementation followed here are outlined in this [paper](https://arxiv.org/pdf/1412.6980.pdf).

    The parameters of this solver are set in `parameters.json`:
    ```json
    "solvers":{
        "adam": {
            "alpha": 0.1, 
            "beta1": 0.9,
            "beta2": 0.999,
            "eps":   1e-8
        }
    }
    ```
    where `alpha`, `beta1`, `beta2` and `eps` are respectly $\alpha$, $\beta_1$, $\beta_2$ and $\epsilon$ in the paper. The values have been set as suggested.

### Test Functions

Several test functions have been implemented in order to test the aforementioned solvers. The implemented functions are listed as follows:

- `default` - this is the function from the challenge assignment.
    $$ f(x,y) = xy + 4x^4 + y^2+ 3x$$
- `rosenbrock` - the Rosenbrock function is defined as follows:
    $$f(x,y)=(a-x)^{2}+b(y-x^{2})^{2}$$
    where $a = 1$ and $b = 100$.
    The global minimum is in $\hat{x} = (1, 1)$ and $f(\hat{x}) = 0$.
- `beale` - the Beale function is defined as follows:
    $$ f(x,y)=\left(1.5-x+xy\right)^{2}+\left(2.25-x+xy^{2}\right)^{2}+ ( 2.625 − x + x y^3 )^2$$
    The global minimum is in $\hat{x} = (3, 0.5)$ and $f(\hat{x}) = 0$.

- `rastrigin` - the Rastrigin function is defined as follows:
$$ f(x) = An + \sum_{i=1}^n [(x_1-1)^2 - A \cos(2\pi (x_i-1))]$$
where A = 10 and n is the dimension of the vector $x$. This is the decentered version, therefore the global minima is $\hat{x} = (1, 1)$ and $f(\hat{x}) = 0$.

Note that all the solvers aforementioned could not find the global minumum for the Rastrigin function due to its large number of local minima.


