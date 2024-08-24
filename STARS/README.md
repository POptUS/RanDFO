*****
# STARS algorithm for stochastic derivative-free optimization
*****
![GitHub](https://img.shields.io/github/license/poptus/randfo) [![miss_hit](https://github.com/POptUS/RanDFO/actions/workflows/miss_hit.yml/badge.svg)](https://github.com/POptUS/RanDFO/actions/workflows/miss_hit.yml)

[STARS](https://doi.org/10.1137/22M1524072) is an extension of [STORM](https://doi.org/10.1007/s10107-017-1141-8) to random subspaces, and a stochastic variant of a simplified version of [RSDFO](https://doi.org/10.1007/s10107-022-01836-1). It is developed for the optimization of stochastically noisy objective functions in large dimensions, that is, whose values can be computed through a blackbox corrupted with random noise. It is a trust-region algorithm which takes as input a vector of decision variables, and which achieves scalability using random models constructed in low-dimensional affine subspaces.

## Prerequisites

* [STARS](https://doi.org/10.1137/22M1524072) is implemented using Matlab R2021a.
* The 40 problems managed by `stars_yatsop.m` are from the [YATSOp](https://github.com/POptUS/YATSOp) repository, which is automatically loaded as a submodule.
* The 53 problems managed by `stars_benchmark.m` are from the [BenDFO](https://github.com/POptUS/BenDFO) repository, which is loaded automatically as a submodule.

## Installation of STARS

* No installation is required. Users simply need to clone the STARS repository making use of the following command from a terminal:

```
git clone --recursive https://github.com/POptUS/RanDFO
```
* The algorithm can then be run using either stars_applications.m, stars_yatsop.m or stars_benchmark.m located in the STARS folder (more details below).

## Getting started

In the STARS folder, users have three options, using either stars_applications.m, stars_yatsop.m or stars_benchmark.m to run [STARS](https://doi.org/10.1137/22M1524072).

* `stars_applications.m` runs [STARS](https://doi.org/10.1137/22M1524072) on unconstrained problems, or problems with bound constraints. It aims to show users how to provide problems to the algorithm. Other very detailed information is provided in the file.
* `stars_yatsop.m` runs [STARS](https://doi.org/10.1137/22M1524072) in an automated way on the 40 problems (from the [YATSOp](https://github.com/POptUS/YATSOp) repository) considered in the numerical section of the [STARS paper](https://doi.org/10.1137/22M1524072), for various subspace dimensions, various types of noise (additive, multiplicative, Gaussian, uniform, etc.), and various noise levels via their standard deviations. It generates solutions files, stats files and history files in a 'stars_outputs' folder, which can be used to generate data profiles, performance profiles, trajectory plots, convergence graphs, etc. Users are referred to the numerical section of the [STARS paper](https://doi.org/10.1137/22M1524072) for more details on the use of this script. Other very detailed information is provided in the stars_yatsop.m file.
* `stars_benchmark.m` runs STARS the same way as `stars_yatsop.m`, but on 53 problems from the [BenDFO](https://github.com/POptUS/BenDFO) repository.

## Citing STARS

If you use [STARS](https://doi.org/10.1137/22M1524072), please cite the following [paper](https://doi.org/10.1137/22M1524072):


```
@article{DzaWildSub2024,
	Author    = {Kwassi Joseph Dzahini and Stefan M. Wild},
	Title     = {Stochastic Trust-Region Algorithm in Random Subspaces with Convergence and Expected Complexity Analyses},
	Journal   = {SIAM Journal on Optimization},
	Year      = {2024},
	volume    = {34},
	number    = {3},
	pages     = {2671--2699},
	doi       = {10.1137/22M1524072}
}
```
[![DOI](https://doi.org/10.1137/22M1524072)](https://doi.org/10.1137/22M1524072)
