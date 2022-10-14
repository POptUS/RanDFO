*****
# STARS algorithm for stochastic derivative-free optimization
*****
[STARS](https://arxiv.org/abs/2207.06452) is an extension of [STORM](https://doi.org/10.1007/s10107-017-1141-8) to random subspaces, and a stochastic variant of a simplified version of [RSDFO](https://doi.org/10.1007/s10107-022-01836-1). It is developed for the optimization of stochastically noisy objective functions in large dimensions, that is, whose values can be computed through a blackbox corrupted with random noise. It is a trust-region algorithm which takes as input a vector of decision variables, and which achieves scalability using random models constructed in low-dimensional affine subspaces.

## Prerequisites

* [STARS](https://arxiv.org/abs/2207.06452) is implemented using Matlab R2021a.
* The 40 problems managed by stars_yatsop.m are from the [YATSOp](https://github.com/POptUS/YATSOp) repository, or available upon request via email: poptus@mcs.anl.gov (Mathematics and Computer Science division of Argonne National Laboratory), in case the latter link is unavailable (that is, the repository is not yet public). But when [YATSOp](https://github.com/POptUS/YATSOp) is public, it is automatically loaded as submodule.
* The 53 problems managed by stars_benchmark.m are from the [BenDFO](https://github.com/POptUS/BenDFO) repository which is loaded automatically as a submodule.

## Installation of STARS

* No installation is required. Users simply need to clone the STARS repository making use of the following command from a terminal:

```
git clone --recursive https://github.com/POptUS/RanDFO
```
* The algorithm can then be run using either stars_applications.m, stars_yatsop.m or stars_benchmark.m located in the STARS folder (more details below).

## Getting started

In the STARS folder, users have three options, using either stars_applications.m, stars_yatsop.m or stars_benchmark.m to run [STARS](https://arxiv.org/abs/2207.06452).

* stars_applications.m runs [STARS](https://arxiv.org/abs/2207.06452) on unconstrained problems, or problems with bound constraints. It aims to show users how to provide problems to the algorithm. Other very detailed information is provided in the file.
* stars_yatsop.m runs [STARS](https://arxiv.org/abs/2207.06452) in an automated way on the 40 problems (from the [YATSOp](https://github.com/POptUS/YATSOp) repository) considered in the numerical section of the [STARS paper](https://arxiv.org/abs/2207.06452), for various subspace dimensions, various types of noise (additive, multiplicative, Gaussian, uniform, etc.), and various noise levels via their standard deviations. It generates solutions files, stats files and history files in a 'stars_outputs' folder, which can be used to generate data profiles, performance profiles, trajectory plots, convergence graphs, etc. Users are referred to the numerical section of the [STARS paper](https://arxiv.org/abs/2207.06452) for more details on the use of this script. Other very detailed information is provided in the stars_yatsop.m file.
* stars_benchmark.m runs STARS the same way as stars_yatsop.m, but on 53 problems from the [BenDFO](https://github.com/POptUS/BenDFO) repository.

## Citing STARS

If you use [STARS](https://arxiv.org/abs/2207.06452), please cite the following [paper](https://arxiv.org/abs/2207.06452):


```
@article{DzaWildSub2022,
	Author      = {K. J. Dzahini and S. M. Wild},
	Title       = {Stochastic trust-region algorithm in random subspaces with convergence and expected complexity analyses},
	Journal = {ArXiv},
	Year        = {2022},
	ArxivUrl    = {https://arxiv.org/abs/2207.06452},
	Url         = {https://arxiv.org/abs/2207.06452},
}
```
[![DOI](https://arxiv.org/abs/2207.06452)](https://arxiv.org/abs/2207.06452)
