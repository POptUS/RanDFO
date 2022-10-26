# RanDFO -- A Collection of Randomized Algorithms for Derivative-Free Optimization
![GitHub](https://img.shields.io/github/license/poptus/randfo) 

Currently available are:

* [STARS](STARS) a randomized algorithm for large-scale derivative-free stochastic optimization  
Implemented in Matlab/Octave and described in [Arxiv:2207.06452](https://arxiv.org/abs/2207.06452) 


## Submodules

RanDFO currently includes dependencies via git submodules. Currently, the following submodules are employed:

* [STARS/problems/](STARS/problems/)BenDFO - A set of test problems
* [STARS/problems/](STARS/problems/)YATSOp - Yet another set of test problems

As a consequence, when cloning the RanDFO repository, the submodules can be retrieved automatically via
- `git clone --recursive` (in place of the usual `git clone`)

If you have already cloned the repository, the following modified `git` commands can be used:
- `git submodule update --init --recursive` (to obtain all submodules and any submodules those submodules have)
- `git pull --recurse-submodules=yes` (in place of the usual `git pull`)
- `git submodule update --init` (additionally, after `git pull`)
- `git submodule update --init STARS/problems/BenDFO`
  (variant of the previous item, in case you want to only get the BenDFO submodule)
  
## Contributing to RanDFO

Contributions are welcome in a variety of forms; please see [CONTRIBUTING](CONTRIBUTING.rst).

## License 

All code included in RanDFO is open source, with the particular form of license contained in the top-level 
subdirectories.  If such a subdirectory does not contain a LICENSE file, then it is automatically licensed 
as described in the otherwise encompassing RanDFO [LICENSE](/LICENSE).  


## Resources

To seek support or report issues, e-mail:

 * ``poptus@mcs.anl.gov``
