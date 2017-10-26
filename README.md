# GMMAT (Generalized Linear Mixed Model Association Test)

GMMAT is an R package for performing genetic association tests for
outcomes with distribution in the exponential family (e.g. binary
outcomes) based on the generalized linear mixed model. It can be used to
analyze genetic data from individuals with population structure and
relatedness. GMMAT fits a generalized linear mixed model under the null
hypothesis of no genetic association, and then performs a score test for
each individual genetic variant. See user manual
[here](https://content.sph.harvard.edu/xlin/dat/GMMAT_user_manual_v0.7.pdf).
A light version of GMMAT which does not depend on the C++ library boost
is available [here](https://github.com/lin-lab/GMMAT_lite) (it cannot
take .gz and .bz2 compressed genotype files).

## References

+ Breslow NE and Clayton DG. (1993) Approximate Inference in Generalized
  Linear Mixed Models. *Journal of the American Statistical Association*
  88: 9-25.
  doi:[10.1080/01621459.1993.10594284](http://dx.doi.org/10.1080/01621459.1993.10594284).
+ Chen H, Wang C, Conomos MP, Stilp AM, Li Z, Sofer T, Szpiro AA, Chen
  W, Brehm JM, Celedon JC, Redline S, Papanicolaou GJ, Thornton TA,
  Laurie CC, Rice K and Lin X. (2016) Control for Population Structure
  and Relatedness for Binary Traits in Genetic Association Studies Using
  Logistic Mixed Models. *The American Journal of Human Genetics*,
  98(4): 653-666.
  doi:[10.1016/j.ajhg.2016.02.012](https://doi.org/10.1016/j.ajhg.2016.02.012).

## License

This software is licensed under GPL-3.
