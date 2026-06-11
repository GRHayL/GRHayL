# Licensing and Citation Guide {#license}

The GRHayL library contains methods and code contributions from several
sources, and these should be properly cited according to usage.
Additionally, the various codes have different licenses, which we document here.

## Licensing

The Noble routines are adapted from the HARM code, which is under the
GNU General Public License.

The Palenzuela and Newman routines are adapted from the
[GRMHD_con2prim](https://bitbucket.org/dsiegel/grmhd_con2prim/src/master/)
code, which is under the Creative Commons Attribution-ShareAlike 4.0
International License. These routines also use a modified version of the toms748
routine from the Boost library, which is subject to the Boost Software License,
Version 1.0.

## Citations

The following chart gives a quick reference guide, followed by a more detailed
description of the citations.

 | Component                             | Additional Citation                                                    |
 |:-------------------------------------:|:----------------------------------------------------------------------:|
 | @grhayl                               | GRHayL release paper and IllinoisGRMHD \cite Etienne_IGM               |
 | Tabulated EOS                         | Werneck _et al_. \cite Werneck_IGM                                     |
 | Noble Con2Prim routines               | HARM \cite HARM and its solver \cite Noble_2006                        |
 | Palenzuela & Newman Con2Prim routines | GRMHD_con2prim \cite Siegel_2018_GRMHD_con2prim and the method \cite Siegel_2018_recovery_schemes |
 | PPM method                            | corrected method/reference papers \cite Marti_1996 and \cite DelZanna_2003 |

Any work using GRHayL should cite its release paper
[arXiv:2512.15846](https://arxiv.org/abs/2512.15846). Additionally,
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
For the PPM method, cite the corrected method/reference papers listed in the
table above.
The original Colella and Woodward PPM reference \cite Colella_1984 contains
known formula errors; later papers and implementation references correct these,
so the table above points PPM users to \cite Marti_1996 and \cite DelZanna_2003.

The WENO-z method is adapted from the Phoebus code developed by
Los Alamos National Lab.
