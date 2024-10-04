*****
# StoDARS algorithm for stochastic derivative-free optimization
*****
![GitHub](https://img.shields.io/github/license/poptus/randfo) [![miss_hit](https://github.com/POptUS/RanDFO/actions/workflows/miss_hit.yml/badge.svg)](https://github.com/POptUS/RanDFO/actions/workflows/miss_hit.yml)

StoDARS is a direct-search algorithm for stochastic optimization that uses random subspaces. 
It is developed for the optimization of stochastically noisy objective functions in large dimensions, that is, whose values can be computed through a blackbox corrupted with random noise. 
It is a direct-search algorithm that takes as input a vector of decision variables, and achieves scalability using pollsets constructed in low-dimensional affine subspaces.
Users may also be interested in [STARS](https://github.com/POptUS/RanDFO/STARS), a trust-region-based algorithm for similar problems.

## Prerequisites

* StoDARS is implemented using Matlab R2021a.
* The 40 problems managed by `stodars_yatsop.m` are from the [YATSOp](https://github.com/POptUS/YATSOp) repository, which is automatically loaded as a submodule.
* The 53 problems managed by `stodars_benchmark.m` are from the [BenDFO](https://github.com/POptUS/BenDFO) repository, which is loaded automatically as a submodule.

## Installation of StoDARS

* No installation is required. Users simply need to clone the StoDARS repository making use of the following command from a terminal:

```
git clone --recursive https://github.com/POptUS/RanDFO
```
* The algorithm can then be run using either stodars_applications.m, stodars_yatsop.m or stodars_benchmark.m located in the StoDARS folder (more details below).

## Getting started

In the StoDARS folder, users have three options, using either stodars_applications.m, stodars_yatsop.m or stodars_benchmark.m to run StoDARS.

* `stodars_applications.m` runs StoDARS on unconstrained problems, or problems with bound constraints. It aims to show users how to provide problems to the algorithm. Other very detailed information is provided in the file.
* `stodars_yatsop.m` runs StoDARS in an automated way on the 40 problems (from the [YATSOp](https://github.com/POptUS/YATSOp) repository) considered in the numerical section of the [StoDARS paper](), for various subspace dimensions, various types of noise (additive, multiplicative, Gaussian, uniform, etc.), and various noise levels via their standard deviations. It generates solutions files, stats files and history files in a 'stodars_outputs' folder, which can be used to generate data profiles, performance profiles, trajectory plots, convergence graphs, etc.
* Users are referred to the numerical section of the [STODARS paper](https://arxiv.org/abs/2403.13320) for more details on the use of this script. Other very detailed information is provided in the stodars_yatsop.m file.
* `stodars_benchmark.m` runs StoDARS the same way as `stodars_yatsop.m`, but on 53 problems from the [BenDFO](https://github.com/POptUS/BenDFO) repository.

## Citing STODARS

If you use StoDARS, please cite the following [paper](https://arxiv.org/abs/2403.13320):


```
@techreport{DzahiniWildDS24,
      title = {Direct Search for Stochastic Optimization in Random Subspaces with Zeroth-, First-, and Second-Order Convergence and Expected Complexity}, 
      author = {Kwassi Joseph Dzahini and Stefan M. Wild},
      year = {2024},
      url = {https://arxiv.org/abs/2403.13320},
      institution = {ArXiv},
       number      = {2403.13320},
}
```
