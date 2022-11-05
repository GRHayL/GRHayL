# Con2Prim_beta
This is a beta version of the infrastructure-agnostic conservative-to-primitive code
my collaborators and I are developing. It is designed to simply take primitives and
conservative structs, perform c2p routines, and return the structs with the reconstructed
primitive data. This code is being developed with the long-term goal of converting
IllinoisGRMHD into a set of libraries which are extensible and easy to plug-n-play
with other codes. This code currently includes pieces of code which will eventually
migrate to other modules, such as an EOS library.

This code does *not* depend on the Einstein Toolkit or Cactus framework, though this
repository is currently being tested within that infrastructure. All functions except
the included con2prim_test.cc have no references to Cactus parameters or header files.

