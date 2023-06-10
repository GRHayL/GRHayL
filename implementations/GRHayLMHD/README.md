# GRHayLMHD

---

GRHayLMHD is a thorn for the Cactus Framework designed
to work in the Einstein Toolkit. It serves as the
replacement for the IllinoisGRMHD thorn using the GRHayL
library.

## Transition From IllinoisGRMHD

Since GRHayLMHD is the future of IllinoisGRMHD, it is
important to provide a smooth transition. As part of this,
we have endeavored to minimize the changes needed for
old parfiles to function with the new code. However,
the cleanup and simplification of the code, along with
new additions like tabulated equation of state, necessitate
some changes.

### Condensing Thorns

Two thorns in WVUThorns only work with IllinoisGRMHD, and
IllinoisGRMHD can't work without them. These are
ID\_converter\_ILGRMHD and Convert\_to\_HydroBase. These
thorns also only consist of single functions, so there is
no reason for them to be separate. GRHayLMHD also includes
these functions.

### ActiveThorn Changes

Before, a simulation with IllinoisGRMHD includes

ID\_converter\_ILGRMHD
Convert\_to\_HydroBase
IllinoisGRMHD

and ususally includes either

Seed\_Magnetic\_Fields

or

Seed\_Magnetic\_Fields\_BNS

The Seed... thorns do not depend on IllinoisGRMHD, so they
can be used with GRHayLMHD. They have been updated to
schedule their functions purely based on HydroBase instead
of the old IllinoisGRMHD-based scheduling. The other 3
must be removed and replaced with the thorns

GRHayLib
GRHayLMHD

If users do not want to evolve magnetic fields, please
consider using the GRHayLHD thorn which explicitly assumes
A\_i = 0. This code will naturally be slightly faster since
there are less grid functions and evolution equations.

### Parameter Changes

Most parameters don't change, so parfiles will simply
need to change the thorn name. However, we give a
detailed list of parameters from the three old thorns and
their equivalent GRHayLMHD parameters.

First, those which have only changed thornnames are as follows:

ID\_converter\_ILGRMHD
random\_seed -> GRHayLMHD
random\_pert -> GRHayLMHD

Convert\_to\_HydroBase
Convert\_to\_HydroBase\_every -> GRHayLMHD

IllinoisGRMHD
rho\_b\_atm   -> GRHayLib
rho\_b\_max   -> GRHayLib
Psi6threshold -> GRHayLib
neos          -> GRHayLib
gamma\_th     -> GRHayLib
damp\_lorenz  -> GRHayLib
verbose       -> GRHayLMHD
update\_Tmunu -> GRHayLMHD
Matter\_BC    -> GRHayLMHD
EM\_BC        -> GRHayLMHD
Symmetry      -> GRHayLMHD
Sym\_Bz       -> GRHayLMHD

The remaining parameters are either deprecated or renamed.
First, the renamed parameters are

ID\_converter\_ILGRMHD
Gamma\_Initial -> GRHayLib::Gamma\_ppoly\_in[0]

IllinoisGRMHD
GAMMA\_SPEED\_LIMIT -> GRHayLib::W\_max
K\_poly             -> GRHayLib::k\_ppoly0

The following parameters are deprecated and should be
removed:

ID\_converter\_ILGRMHD
pure\_hydro\_run -> Use the GRHayLHD thorn
K\_Initial       -> Since the thorns are combined, this and
                    K\_poly are duplicates

IllinoisGRMHD
tau\_atm -> tau\_atm is now automatically computed by GRHayLib;
            to set it manually, simply set the value of
            grhayl\_eos->tau\_atm
            after the GRHayLib initialization function has run.

conserv\_to\_prims\_debug -> this section of the con2prim
            routine is removed; if users feel this is particularly
            useful, it could be added in a future update

### Initial Data and EOS\_Omni

Many thorns use EOS\_Omni, but the GRHayL library has its own
EOS. To use these ID thorns (like the Meudon thorns), simply
set up EOS\_Omni as usual and make sure that GRHayLib is set
to have the same EOS. This behavior is mechanically identical
to IllinoisGRMHD, which also handled the EOS internally instead
of relying on EOS\_Omni.

## Transition Summary
1) Activate the thorns GRHayLib and GRHayLMHD and deactivate
   IllinoisGRMHD thorns

2) Change parameters:
   ID\_converter\_ILGRMHD::Gamma\_Initial -> GRHayLib::Gamma\_ppoly\_in[0]
   ID\_converter\_ILGRMHD::random\_seed   -> GRHayLMHD::random\_seed
   ID\_converter\_ILGRMHD::random\_pert   -> GRHayLMHD::random\_pert
   
   Convert\_to\_HydroBase::Convert\_to\_HydroBase\_every -> GRHayLMHD::Convert\_to\_HydroBase\_every
   
   IllinoisGRMHD::GAMMA\_SPEED\_LIMIT -> GRHayLib::W\_max
   IllinoisGRMHD::K\_poly             -> GRHayLib::k\_ppoly0
   IllinoisGRMHD::rho\_b\_atm         -> GRHayLib::rho\_b\_atm
   IllinoisGRMHD::rho\_b\_max         -> GRHayLib::rho\_b\_max
   IllinoisGRMHD::Psi6threshold       -> GRHayLib::Psi6threshold
   IllinoisGRMHD::neos                -> GRHayLib::neos
   IllinoisGRMHD::gamma\_th           -> GRHayLib::gamma\_th
   IllinoisGRMHD::damp\_lorenz        -> GRHayLib::damp\_lorenz
   IllinoisGRMHD::verbose             -> GRHayLMHD::verbose
   IllinoisGRMHD::update\_Tmunu       -> GRHayLMHD::update\_Tmunu
   IllinoisGRMHD::Matter\_BC          -> GRHayLMHD::Matter\_BC
   IllinoisGRMHD::EM\_BC              -> GRHayLMHD::EM\_BC
   IllinoisGRMHD::Symmetry            -> GRHayLMHD::Symmetry
   IllinoisGRMHD::Sym\_Bz             -> GRHayLMHD::Sym\_Bz

3) Removed deprecated parameters:
   ID\_converter\_ILGRMHD::K\_Initial
   ID\_converter\_ILGRMHD::pure\_hydro\_run
   IllinoisGRMHD::tau\_atm
   IllinoisGRMHD::conserv\_to\_prims\_debug

4) If doing pure hydro simulations, change GRHayLMHD to GRHayLHD
   and remove EM\_BC and Sym\_Bz parameters
   
