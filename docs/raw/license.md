# Licensing and Citation Guide {#license}

The GRHayL library contains methods and code contributions from several
sources, and these should be properly cited according to usage.
Additionally, the various codes have different licenses, which we document here.

## Licensing

The Noble routines are adapted from the HARM code, which is under the
GNU General Public License.

The Palenzuela and Newman routines are adapted from the [GRMHD_con2prim](https://bitbucket.org/dsiegel/grmhd_con2prim/src/master/) code, which is under the Creative
Commons Attribution-ShareAlike 4.0 International License.

## Citations

The following chart gives a quick reference guide, followed by a more detailed
description of the citations.

 | Component                             | Additional Citation                                                    |
 |:-------------------------------------:|:----------------------------------------------------------------------:|
 | @grhayl                               | TBD and IllinoisGRMHD \cite Etienne_IGM                                |
 | Tabulated EOS                         | Werneck _et al_. \cite Werneck_IGM                                     |
 | Noble Con2Prim routines               | HARM \cite HARM and its solver \cite Noble_2006                        |
 | Palenzuela & Newman Con2Prim routines | GRMHD_con2prim \cite Siegel_2018_code and the method \cite Siegel_2018 |
 | PPM method                            | method paper \cite Siegel_2018_code                                    |

Any work using GRHayL should cite its release paper [TBD](). Additionally,
all components of the library contain elements adapted from IllinoisGRMHD
and should cite its [release paper](https://iopscience.iop.org/article/10.1088/0264-9381/32/17/175009)
([arxiv](https://arxiv.org/abs/1501.07276)) as well.

### EOS citations:
The tabulated EOS is adapted from the Werneck _et al_. fork of IllinoisGRMHD,
and its [release paper](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.107.044037)
([arxiv](https://arxiv.org/abs/2208.14487)) should be cited.

### Con2Prim citations:
The Noble routines are adapted from the HARM code, and the code papers
[here](https://iopscience.iop.org/article/10.1086/374594)
([arxiv](https://arxiv.org/abs/astro-ph/0301509)) and
[here](https://iopscience.iop.org/article/10.1086/500349)
([arxiv](https://arxiv.org/abs/astro-ph/0512420)) should be cited.

The Palenzuela and Newman routines are adapted from the
[GRMHD_con2prim](https://bitbucket.org/dsiegel/grmhd_con2prim/src/master/) code, and
the [code](https://zenodo.org/records/1213306) and associated
[method paper](https://iopscience.iop.org/article/10.3847/1538-4357/aabcc5)
([arxiv](https://arxiv.org/abs/1712.07538)) should be cited.

### Reconstruction citations:
The PPM method is adapted from the work by Del Zanna et al, and the
[release paper](https://www.aanda.org/articles/aa/abs/2003/11/aa3107/aa3107.html)
([arxiv](https://arxiv.org/abs/astro-ph/0210618)) should be cited.

The WENO-z method is adapted from the Phoebus code developed by
Los Alamos National Lab.
