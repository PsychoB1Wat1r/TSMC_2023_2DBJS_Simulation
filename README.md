## FTROA networked filtering of 2DBJS

This repository contains the code implementation of our TSMC_2023 paper [here](https://ieeexplore.ieee.org/abstract/document/10287368).

If you find this repository useful, please cite our paper.
```
  @article{zhou2023fault,
  title={Fault-Tolerant Reduced-Order Asynchronous Networked Filtering of {2-D Bernoulli} Jump Systems},
  author={Zhou, Jianping and Ma, Xiaofeng and Yan, Zhilian and Ahn, Choon Ki},
  journal={IEEE Transactions on Systems, Man, and Cybernetics: Systems},
  volume={54},
  number={2},
  pages={891--902},
  year={2024},
  publisher={IEEE}}
```
## Explanation:
We use **MATLAB** as the programming language and employ the **Mosek solver** along with the **YALMIP toolbox** as development tools to perform numerical solution and simulation for the proposed **FTROA** method.

## Mosek Solver:
Mosek is a high-performance mathematical optimization solver specifically designed for convex optimization problems. For official documentation, please refer to this [link](https://www.mosek.com/documentation/).

## YALMIP Toolbox:
YALMIP (Yet Another LMI Parser) is an optimization modeling toolbox for MATLAB that facilitates the formulation and solution of various optimization problems. Developed by [**Johan Löfberg**](https://scholar.google.com/citations?user=No-9sDUAAAAJ&hl=en), YALMIP is primarily used for convex optimization but can also handle certain non-convex problems.

It is important to note that YALMIP itself is **not a solver**; rather, it serves as a modeling interface that translates user-defined optimization problems into a standard form and then calls external solvers such as **Mosek**, **Gurobi**, **SDPT3**, **SeDuMi**, and others for solution. For official documentation, please refer to this [link](https://yalmip.github.io/).
