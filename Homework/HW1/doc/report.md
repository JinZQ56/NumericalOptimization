## Homework_1: Linear-search Steepest Gradient Descent

#### 1. Workflow

* implement a class named Rosenbrock  which can instantiate Rosenbrock function with any dimension and calculate the value and gradient of such function
* implement a function which optimize the given Rosenbrock function with Linear-search Steepest Gradient Descent

#### 2. Result

* Visualization of a 2D Rosenbrock function

![Rosenbrock2D](/run/user/1000/doc/fe6633fa/Rosenbrock2D.png) 

* Optimization result of a 2D RB function and 10D RB function

  ![Result10D](/run/user/1000/doc/1b21334a/Result10D.png)

#### 3. Analysis

* Armijo condition can help the optimization algorithm a lot. Without using of Armijo condition, it is hard to converge to the reult[1, 1, ..., 1].