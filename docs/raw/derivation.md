/*! @page derivation GRMHD Derivation

Re-derived and written by Samuel Cupp

\section toc Table of Contents
1. @ref intro
2. @ref tmunu
   1. @ref tmunu_term_1
   2. @ref tmunu_term_2
   3. @ref tmunu_final
3. @ref cons_eq
   1. @ref dens_eq
   2. @ref mom_eq
   3. @ref energy_eq
   3. @ref mag_eq
   3. @ref eq_summary
4. @ref appendix
   1. @ref levi_civita
   2. @ref curv_id

---

\section intro Introduction

When considering the evolution of magnetohydrodynamics in a general
relativistic setting, there are generally two sets of variables that are used.
These are the primitive variables

\f[
\mathbf{P} = \left( \rho_b, P, v^i, B^i \right)
\f]

and the conservative variables

\f[
\mathbf{C} = \left( \rho_*, \tilde{\tau}, \tilde{S}_i, \tilde{B}^i \right).
\f]

The primitive variables are the more 'natural' choice, and most equations are
defined in terms of these variables, like expressions for the stress-energy
tensor. However, evolution of \f$ \mathbf{P} \f$ directly is unstable, and so a
more reliable choice of evolution variables is required. This leads to the
conservative variables. Evolution of \f$ \mathbf{C} \f$ also has the benefit of
automatically ensuring that the constraints are satisfied, but these variables
come at the cost of requiring complicated reconstruction of the primitive
variables. GRHayL provides functionality to facilitate the evolution of
\f$ \mathbf{C} \f$, and so I derive these core evolution equations here. The
derivation in this document is based on the work by \cite Duez_IGM, and I have
simply re-derived their results. The only deviation from the original paper is
a rescaling we apply to \f$ B^i \f$, which is also described below.

---

\section tmunu Deriving the stress-energy tensor

In order to find evolution equations for the conserved quantities for the
system, we start by finding a convenient form of the electromagnetic
energy-momentum tensor \f$ T_{EM}^{\mu\nu} \f$, which is defined as

\f[
T_{EM}^{\mu\nu} = \frac{1}{4\pi}\left( F^{\mu\lambda}F^{\nu}_{\hphantom{\nu}\lambda}
                                   - \frac{1}{4}g^{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}
                                \right)
\f]

where \f$ F^{\mu\nu} \f$ is the Faraday tensor. The magnetic field is given by

\f[
\mathcal{B}^\mu = \frac{1}{2} \epsilon^{\mu\nu\beta\alpha}n_\nu F_{\alpha\beta}
                = n_\nu F^{*\nu\mu},
\f]

where

\f[
F^{*\mu\nu} = \frac{1}{2} \epsilon^{\mu\nu\alpha\beta}F_{\alpha\beta}
\f]

is the dual of the Faraday tensor. **However, in GRHayL we instead use the
rescaled magnetic field**

\f[
B^\mu = \frac{1}{\sqrt{4\pi}}\mathcal{B}^\mu.
\f]

By rescaling \f$ B^\mu \f$ by a factor of \f$ \sqrt{4\pi} \f$, we remove this
factor from all subsequent hydrodynamic evolution equations. As the induction
equations relate \f$ B^i \f$, \f$ E^i \f$, and \f$ A_i \f$, rescaling all
quantities results in identical evolution equations for the magnetic sector
regardless of this choice. This effectively removes many unneeded mathematical
operations that would occur on every point, at every time step, for nearly
every hydrodynamic function. For the remainder of this documentation--and
everywhere else in the GRHayL documentation--magnetic quantities are assumed to
be rescaled by this factor.

Now, we want to examine this with respect to an arbitrary observer with
normalized four-velocity \f$ \xi^\mu \f$. If we consider the electromagnetic
fields in the rest frame of this observer,

\f[
\xi_\mu E^\mu = 0 = \xi_\mu B^\mu
\f]

since \f$ E \f$ and \f$ B \f$ are purely spatial. The Faraday tensor can then
be decomposed in terms of \f$ E \f$ and \f$ B \f$:

\f[
F^{\mu\nu} = \sqrt{4\pi} \left( \xi^\mu E^\nu
                              - \xi^\nu E^\mu
                              + \xi_\gamma \epsilon^{\gamma\mu\nu\delta} B_{\delta}
                         \right).
\f]

where \f$ \epsilon^{\gamma\mu\nu\delta} \f$ is the Levi-Civita tensor. As many
elements of the subsequent derivation revolves around ~~abusing~~ taking
advantage of various identities and symmetries of this tensor, I refer readers
to the [Appendix](@ref levi_civita) for more details about these identities.

---

\subsection tmunu_term_1 Solving for F^2

Now, we simply want to find an expression for of \f$ T_{EM}^{\mu\nu} \f$ in terms
of \f$ \xi \f$, \f$ E \f$, and \f$ B \f$. Let's first consider the relatively
simpler term

\f[
\begin{aligned}
\frac{1}{\sqrt{4\pi}} F_{\alpha\beta}F^{\alpha\beta}
&= \left(\xi_\alpha E_\beta
       - \xi_\beta E_\alpha
       + \xi_\gamma \epsilon^{\gamma\hphantom{\alpha\beta}\delta}_{\hphantom{\gamma}\alpha\beta} B_{\delta}
   \right)
   \left(\xi^\alpha E^\beta
       - \xi^\beta E^\alpha
       + \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} B_{\mu}
   \right) \\

&= \xi_\alpha \xi^\alpha E_\beta E^\beta
 - \xi_\alpha E^\alpha E_\beta \xi^\beta
 - E_\alpha \xi^\alpha \xi_\beta E^\beta
 + \xi_\beta \xi^\beta E_\alpha E^\alpha
 + \xi_\alpha E_\beta \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} B_{\mu}
 - \xi_\beta E_\alpha \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} B_{\mu} \\
&\hspace{5mm}
 + \xi^\alpha E^\beta \xi_\gamma \epsilon^{\gamma\hphantom{\alpha\beta}\delta}_{\hphantom{\gamma}\alpha\beta} B_{\delta}
 - \xi^\beta E^\alpha \xi_\gamma \epsilon^{\gamma\hphantom{\alpha\beta}\delta}_{\hphantom{\gamma}\alpha\beta} B_{\delta}
 + \xi_\sigma \xi_\gamma \epsilon^{\gamma\hphantom{\alpha\beta}\delta}_{\hphantom{\gamma}\alpha\beta} \epsilon^{\sigma\alpha\beta\mu} B_{\delta} B_{\mu} \\ 

&= - 2E_\beta E^\beta
   + 4\xi_\alpha \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} E_\beta B_{\mu}
   + \xi_\sigma \xi^\gamma \epsilon_{\gamma\alpha\beta\delta} \epsilon^{\sigma\alpha\beta\mu} B^{\delta} B_{\mu}
\end{aligned}
\f]

where we have used \f$ \xi^\mu \xi_\mu=-1 \f$ and the permutation rules of the
Levi-Civita tensor. By consideration of the properties of the Levi-Civita
tensor, we can further simplify this expression:

\f[
\begin{aligned}
\frac{1}{\sqrt{4\pi}} F_{\alpha\beta}F^{\alpha\beta}
&= - 2E_\beta E^\beta
   + 4\xi_\alpha \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} E_\beta B_{\mu}
   + 2\xi_\sigma \xi^\gamma (\delta^\sigma_\delta \delta^\mu_\gamma
                           - \delta^\sigma_\gamma \delta^\mu_\delta) B^{\delta} B_{\mu} \\

&= - 2E_\beta E^\beta
   + 4\xi_\alpha \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} E_\beta B_{\mu}
   + 2\xi_\sigma \xi^\mu B^{\sigma} B_{\mu}
   - 2\xi_\sigma \xi^\sigma B^{\mu} B_{\mu} \\

&= - 2E_\beta E^\beta
   + 4\xi_\alpha \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} E_\beta B_{\mu}
   + 2B^{\mu} B_{\mu}
\end{aligned}
\f]

Now, since we are in the rest frame of the observer, the spatial components of
the observer's velocity are zero. Then,

\f[
\begin{align}
\xi_\alpha \xi_\sigma \epsilon^{\sigma\alpha\beta\mu}
&= \xi_0 \xi_\sigma \epsilon^{\sigma 0\beta\mu} \\
&= \xi_0 \xi_0 \epsilon^{0 0\beta\mu} \\
&= 0
\end{align}
\f]

since any components of the Levi-Civita tensor with repeated indices are zero.
Therefore,

\f[
\frac{1}{\sqrt{4\pi}} F_{\alpha\beta}F^{\alpha\beta} = -2E^\beta E_\beta + 2B^{\mu} B_{\mu}
\f]

---

\subsection tmunu_term_2 Solving for F^{ai} F^b_i

In this section, we tackle the term \f$ F^{\mu\lambda}F^{\nu}_{\hphantom{\nu}\lambda} \f$
in the energy-momentum tensor.

\f[
\begin{aligned}
\frac{1}{\sqrt{4\pi}} F^{\mu\lambda}F^{\nu}_{\hphantom{\nu}\lambda}
&= \left(\xi^\mu E^\lambda
       - \xi^\lambda E^\mu
       + \xi_\gamma \epsilon^{\gamma\mu\lambda\delta} B_{\delta}
   \right)
   \left(\xi^\nu E_\lambda
       - \xi_\lambda E^\nu
       + \xi_\alpha \epsilon^{\alpha\nu\hphantom{\lambda}\beta}_{\hphantom{\alpha\nu}\lambda} B_{\beta}
   \right) \\

&= \xi^\mu E^\lambda \xi^\nu E_\lambda
 - \xi^\mu E^\lambda \xi_\lambda E^\nu
 - \xi^\lambda E^\mu \xi^\nu E_\lambda
 + \xi^\lambda E^\mu \xi_\lambda E^\nu
 + \xi^\mu E_\lambda \xi_\alpha \epsilon^{\alpha\nu\lambda\beta} B_{\beta} \\
&\hspace{5mm}
 - \xi_\lambda E^\mu \xi_\alpha \epsilon^{\alpha\nu\lambda\beta} B_{\beta}
 + \xi_\gamma \epsilon^{\gamma\mu\lambda\delta} B_{\delta} \xi^\nu E_\lambda
 - \xi_\gamma \epsilon^{\gamma\mu\lambda\delta} B_{\delta} \xi_\lambda E^\nu
 + \xi_\gamma \epsilon^{\gamma\mu\lambda\delta} B_{\delta} \xi_\alpha \epsilon^{\alpha\nu\hphantom{\lambda}\beta}_{\hphantom{\alpha\nu}\lambda} B_{\beta} \\

&= \xi^\mu \xi^\nu E^\lambda E_\lambda
 - E^\mu E^\nu
 + \xi_\alpha E_\lambda B_{\beta} \left(\xi^\mu \epsilon^{\alpha\nu\lambda\beta}
                                      + \xi^\nu \epsilon^{\alpha\mu\lambda\beta}
                                  \right) \\
&\hspace{5mm}
 - \xi_\lambda \xi_\alpha B_{\beta} (\epsilon^{\alpha\nu\lambda\beta} E^\mu
 + \epsilon^{\alpha\mu\lambda\beta} E^\nu)
 + \xi_\alpha \xi_\gamma \epsilon^{\gamma\mu\lambda\delta} \epsilon^{\alpha\nu\hphantom{\lambda}\beta}_{\hphantom{\alpha\nu}\lambda} B_{\delta} B_{\beta}
\end{aligned}
\f]

where we have again used the relations \f$ \xi^\mu \xi_\mu=-1 \f$ and
\f$ \xi_\mu E^\mu = 0 = \xi_\mu B^\mu \f$. We've also taken the liberty to
re-label dummy indices to simplify the resulting expression. To simplify the
final term, we must (unfortunately) once again rely on the generalized
Kronecker delta. This time we use the second (far more complicated) expression
in the [Appendix](@ref levi_civita):

\f[
\begin{aligned}
\xi_\alpha \xi_\gamma \epsilon^{\gamma\mu\lambda\delta} \epsilon^{\alpha\nu\hphantom{\lambda}\beta}_{\hphantom{\alpha\nu}\lambda} B_{\delta} B_{\beta}
&= -g^{\nu\phi} \xi^\alpha \xi_\gamma B_{\delta} B^{\beta} \left(\delta^\gamma_\alpha \left(\delta^\mu_\phi \delta^\delta_\beta
                                                                                          - \delta^\mu_\beta \delta^\delta_\phi
                                                                                      \right)
                                                                - \delta^\gamma_\phi \left(\delta^\mu_\alpha \delta^\delta_\beta
                                                                                         - \delta^\mu_\beta \delta^\delta_\alpha
                                                                                     \right)
                                                                + \delta^\gamma_\beta \left(\delta^\mu_\alpha \delta^\delta_\phi
                                                                                          - \delta^\mu_\phi \delta^\delta_\alpha
                                                                                      \right)
                                                           \right) \\

&= - g^{\nu\phi} \xi^\alpha \xi_\alpha B_{\delta} B^{\beta} \left(\delta^\mu_\phi \delta^\delta_\beta
                                                                - \delta^\mu_\beta \delta^\delta_\phi
                                                            \right)
   + g^{\nu\gamma} \xi^\alpha \xi_\gamma B_{\delta} B^{\beta} \left(\delta^\mu_\alpha \delta^\delta_\beta
                                                                  - \delta^\mu_\beta \delta^\delta_\alpha
                                                              \right)
   - g^{\nu\phi} \xi^\alpha \xi_\beta B_{\delta} B^{\beta} \left(\delta^\mu_\alpha \delta^\delta_\phi
                                                               - \delta^\mu_\phi \delta^\delta_\alpha
                                                           \right) \\

&= - g^{\mu\nu} \xi^\alpha \xi_\alpha B_{\beta} B^{\beta}
   + g^{\nu\delta} \xi^\alpha \xi_\alpha B_{\delta} B^{\mu}
   + g^{\nu\gamma} \xi^\mu \xi_\gamma B_{\beta} B^{\beta}
   - g^{\nu\gamma} \xi^\alpha \xi_\gamma B_{\alpha} B^{\mu}
   - g^{\nu\delta} \xi^\mu \xi_\beta B_{\delta} B^{\beta}
   + g^{\mu\nu} \xi^\alpha \xi_\beta B_{\alpha} B^{\beta} \\

&= g^{\mu\nu} B_{\beta} B^{\beta}
 - B^{\mu} B^{\nu}
 + \xi^\mu \xi^\nu B_{\beta} B^{\beta} \\
\end{aligned}
\f]

which leaves us with

\f[
\begin{aligned}
\frac{1}{\sqrt{4\pi}} F^{\mu\lambda}F^{\nu}_{\hphantom{\nu}\lambda}
&= \xi^\mu \xi^\nu E^\lambda E_\lambda
 - E^\mu E^\nu
 + \xi_\alpha E_\lambda B_{\beta} \left(\xi^\mu \epsilon^{\alpha\nu\lambda\beta}
                                      + \xi^\nu \epsilon^{\alpha\mu\lambda\beta}
                                  \right) \\
&\hspace{5mm}
 - \xi_\lambda \xi_\alpha B_{\beta} \left(\epsilon^{\alpha\nu\lambda\beta} E^\mu
                                        + \epsilon^{\alpha\mu\lambda\beta} E^\nu
                                    \right)
 + g^{\mu\nu} B_{\beta} B^{\beta}
 - B^{\mu} B^{\nu}
 + \xi^\mu \xi^\nu B_{\beta} B^{\beta} \\

&= \xi^\mu \xi^\nu E^\lambda E_\lambda
 + \left(g^{\mu\nu}
 + \xi^\mu \xi^\nu\right) B_{\beta} B^{\beta}
 - E^\mu E^\nu
 - B^{\mu} B^{\nu}
 + \xi_\alpha E_\lambda B_{\beta} \left(\xi^\mu \epsilon^{\alpha\nu\lambda\beta}
                                      + \xi^\nu \epsilon^{\alpha\mu\lambda\beta}
                                  \right)
\end{aligned}
\f]

where we have again used the fact that
\f$ \xi_\alpha \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} = 0 \f$.

---

\subsection tmunu_final Finishing Touches

Now we just need to put all these pieces together. The final expression for
the electromagnetic stress-energy tensor for an arbitrary observer is

\f[
\begin{aligned}
T_{EM}^{\mu\nu}|_{\xi\rightarrow n}
&= \xi^\mu \xi^\nu E^\lambda E_\lambda
 + \left(g^{\mu\nu} + \xi^\mu \xi^\nu\right) B_{\beta} B^{\beta}
 - E^\mu E^\nu
 - B^{\mu} B^{\nu} \\
&\hspace{5mm}
 + \xi_\alpha E_\lambda B_{\beta} \left(\xi^\mu \epsilon^{\alpha\nu\lambda\beta}
                                      + \xi^\nu \epsilon^{\alpha\mu\lambda\beta}
                                  \right)
 - \frac{1}{4}g^{\mu\nu}\left(-2E^\beta E_\beta + 2B^{\mu} B_{\mu} \right) \\

&= \frac{1}{2}\left(g^{\mu\nu}
 + 2\xi^\mu \xi^\nu\right)\left(E^\alpha E_\alpha + B^{\alpha} B_{\alpha}\right)
 - \left(E^\mu E^\nu + B^{\mu} B^{\nu}\right)
 + \xi_\alpha E_\lambda B_{\beta} \left(\xi^\mu \epsilon^{\alpha\nu\lambda\beta}
                                      + \xi^\nu \epsilon^{\alpha\mu\lambda\beta}
                                  \right)
\end{aligned}
\f]

As a reminder, all the factors of \f$ \frac{1}{4\pi} \f$ have been absorbed by
\f$ E \f$ and \f$ B \f$ due to our choice of rescaling. Now, we can consider
two special cases. First, let's consider the normal observer, in which case
the normalized 4-velocity becomes the unit normal vector \f$ \xi\rightarrow n \f$

\f[
\begin{aligned}
n^\nu &= \left( \frac{1}{\alpha}, -\frac{\beta^i}{\alpha} \right) \\
n_\nu &= \left( -\alpha, 0 \right).
\end{aligned}
\f]

If we also define

\f[
\xi_\alpha \epsilon^{\alpha\nu\lambda\beta} \equiv \epsilon^{\nu\lambda\beta}
\f]

then the stress-energy tensor becomes

\f[
T_{EM}^{\mu\nu}|_{\xi\rightarrow n} = \frac{1}{2}\left(g^{\mu\nu} + 2n^\mu n^\nu\right)
                                      \left(E^\alpha E_\alpha + B_{\alpha} B^{\alpha}\right)
                                    - \left(E^\mu E^\nu + B^{\mu} B^{\nu}\right)
                                    + E_\alpha B_{\beta} n^{(\mu} \epsilon^{\nu)\alpha\beta}
\f]

We can also consider the observer co-moving with the fluid. We can call this
co-moving velocity \f$ u_\mu \f$. In the case of ideal conditions (that is to
say, a fluid with perfect conductivity), we have (from Ohm's law) that

\f[
u_\mu F^{\mu\nu} = 0
\f]

which implies that

\f[
E^\mu = \frac{1}{\sqrt{4\pi}} u_\nu F^{\mu\nu} = 0
\f]

Since \f$ E^\mu = 0 \f$, the stress-energy tensor simplifies significantly
in this frame:

\f[
T_{EM}^{\mu\nu}|_{\xi\rightarrow u} = \frac{1}{2}\left(g^{\mu\nu}
                                                     + 2u^\mu u^\nu
                                                 \right) B_{\beta} B^{\beta}
                                    - B^{\mu} B^{\nu}
\f]

For clarity, we henceforth refer to the magnetic field in the co-moving frame
as \f$ B^\mu_{(u)} \f$. The more common expression for the stress-energy tensor
defines a variable \f$ b^\mu \f$. However, this variable is identical to our
rescaled magnetic field in the co-moving frame:

\f[
b^\mu = B^\mu_{(u)}
\f]

In terms of \f$ b \f$, the electromagnetic energy-momentum tensor is given by

\f[
T_{EM}^{\mu\nu}|_{\xi\rightarrow u} = b^2\left(\frac{1}{2}g^{\mu\nu} + u^\mu u^\nu\right)
                                    - b^{\mu} b^{\nu}
\f]

---

\section cons_eq Evolution Equations for the Conservative Variables

Before discussing the evolution variables, recall the 3+1 decomposition which
is used for GRMHD evolution. The metric is given by the line element

\f[
ds^2 = -\alpha^2 dt^2
     + \gamma_{i j}\left(dx^i + \beta^j dt\right)\left(dx^j + \beta^i dt\right)
\f]

where \f$ \alpha \f$, \f$ \beta^i \f$, and \f$ \gamma_{i j} \f$ are the lapse,
shift, and spatial 3-metric, respectively. These quantities are used in the
definitions of the conservative variables which we will soon introduce.

In GRMHD simulations, there are two sets of variables which are used during
evolution. The first set, called primitive variables, are

\f[
\mathbf{P} = \left( \rho_b, P, v^i, B^i \right)
\f]

where \f$ \rho_b \f$ is the (baryonic) fluid rest-mass density, \f$ P \f$ is
the pressure, \f$ v^i = u^i/u^0 \f$ is the rescaled three-velocity, and
\f$ B^i \f$ is the magnetic field. An important sidenote here is that
different GRMHD codes define different velocities. As an example, the widely
used `GRHydro` thorn in the Einstein Toolkit uses the Valencia formalism defined by

\f[
\tilde{v}^i = \frac{u^i}{W} + \frac{\beta^i}{\alpha}
\f]

where \f$ W=\sqrt{1-\tilde{v}^i \tilde{v}_i} \f$ is the Lorentz factor. The
quantity \f$ v^i \f$ used by `IllinoisGRMHD` is **not** the Valencia
three-velocity, so those interested in seeing more about the Valencia
formulation of the evolution equations can look to GRHydro\cite GRHydro and
its references.

The second set, the conservative variables, are

\f[
\mathbf{C} = \left( \tilde{D}, \tilde{\tau}, \tilde{S}_i, \tilde{B}^i \right)
\f]

where

\f[
\mathbf{C} = 
\begin{bmatrix}
  \tilde{D} \\
  \tilde{\tau} \\
  \tilde{S}_{i} \\
  \tilde{B}^{i}
\end{bmatrix}
\equiv
\sqrt{\gamma}
\begin{bmatrix}
  W \rho_b \\
  \tau \\
  S_{i} \\
  B^{i}
\end{bmatrix}
= \sqrt{\gamma}
\begin{bmatrix}
  W \rho_b \\
  \alpha^{2} T^{00} - W \rho_b \\
  \alpha T^{0}_{i} \\
  B^{i}
\end{bmatrix}
\f]

The variable \f$ \tilde{D} \equiv \sqrt{\gamma} W \rho_b \f$ is also referred
to as \f$ \rho_* \f$, and both notations are used in the literature. All of
our previous derivation, of course, only considered the electromagnetic
contribution. We do not derive the hydrodynamic part of the stress-energy
tensor and instead merely state that in the ideal MHD limit the full
stress-energy tensor takes the form

\f[
T^{\mu\nu} = \rho_b h u^\mu u^\nu + P g^{\mu\nu} + T_{EM}^{\mu\nu}
\f]

where \f$ h = 1 + \epsilon + P/\rho_b \f$ is the enthalpy and \f$ \epsilon \f$
is the specific internal energy. Our next step is to derive the actual
evolution equations for our evolution variables \f$ \mathbf{C} \f$. To this
end, we will use the identity

\f[
\Gamma^\nu_{\gamma\nu} = \frac{1}{\sqrt{-g}} \partial_\gamma\sqrt{-g}
\f]

to simplify the equations. We will be frequently using this to condense two
terms into one in the following way:

\f[
\begin{aligned}
\partial_\nu T_\mu^\nu + \Gamma^\nu_{\gamma\nu} T^\gamma_\mu
&= \partial_\nu T_\mu^\nu + \frac{1}{\sqrt{-g}} T^\gamma_\mu \partial_\gamma\sqrt{-g} \\
&= \partial_\nu T_\mu^\nu + \frac{1}{\sqrt{-g}} T^\gamma_\mu \partial_\gamma\sqrt{-g} \\

&= \frac{1}{\sqrt{-g}} \left(\sqrt{-g} \partial_\nu T_\mu^\nu + T^\gamma_\mu \partial_\gamma\sqrt{-g}\right) \\

&= \frac{1}{\sqrt{-g}} \partial_\gamma \left(\sqrt{-g}T^\gamma_\mu\right)
\end{aligned}
\f]

---

\subsection dens_eq Evolution Equation for the Fluid Density

The first evolution equation comes from the baryon number conservation equation

\f[
\begin{aligned}
0 &= \nabla_\nu\left( \rho_b u^\nu \right) \\

&= \partial_\nu\left( \rho_b u^\nu \right)
 + \Gamma^\nu_{\gamma\nu} \rho_b u^\gamma \\

&= \partial_\nu\left( \rho_b u^\nu \right)
 + \rho_b u^\gamma \frac{1}{\sqrt{-g}} \partial_\gamma\sqrt{-g} \\

&= \frac{1}{\sqrt{-g}}\partial_\nu \left(\sqrt{-g}\rho_b u^\nu \right) \\

&= \partial_t \left( \alpha\sqrt{\gamma}\rho_b u^0 \right)
 + \partial_i \left( \alpha\sqrt{\gamma}\rho_0 u^i \right) \\

&= \partial_t \left( \sqrt{\gamma} W \rho_b \right)
 + \partial_i \left( \sqrt{\gamma} u^i \frac{W}{u^0} \rho_0 \right) \\

&= \partial_t \tilde{D} + \partial_i \left( \frac{\tilde{D}}{u^0} u^i \right) \\

&= \partial_t \tilde{D} + \partial_i \left( \tilde{D} v^i \right)
\end{aligned}
\f]

---

\subsection mom_eq Evolution Equation for Momentum Density Variable

The next two sections derive their equations from the energy-momentum
conservation equation

\f[
\begin{aligned}
\nabla_\nu T_\mu^\nu &= 0 \\

&= \partial_\nu T_\mu^\nu
 + \Gamma^\nu_{\gamma\nu} T^\gamma_\mu
 - \Gamma^\gamma_{\mu\nu} T^\nu_\gamma
\end{aligned}
\f]

Rearranging and using the identity for \f$ \Gamma^\nu_{\gamma\nu} \f$, we get

\f[
\frac{1}{\sqrt{-g}}\partial_\nu \left(\sqrt{-g} T_\mu^\nu \right)
      = \Gamma^\gamma_{\mu\nu} T^\nu_\gamma
\f]

For the momentum density variable \f$ \tilde{S}_i \f$, we consider the indices
\f$ \mu\in\{1,2,3\} \f$. Replacing \f$ \mu \f$ with \f$ i \f$ to find the
evolution equation,

\f[
\begin{aligned}
\frac{1}{\sqrt{-g}}\partial_\nu \left(\sqrt{-g} T_i^\nu \right)
&= \Gamma^\gamma_{i\nu} T^\nu_\gamma \\

\frac{1}{\sqrt{-g}}\partial_t \left(\sqrt{-g} T_i^0 \right)
 + \frac{1}{\sqrt{-g}}\partial_j \left(\sqrt{-g} T_i^j \right)
&= \frac{1}{2}T^\nu_\gamma g^{\gamma\beta}\left( g_{\beta\nu, i}
 + g_{i\beta, \nu} - g_{i\nu, \beta} \right) \\

\partial_t \tilde{S}_i + \partial_j \left(\alpha\sqrt{\gamma} T_i^j \right)
&= \frac{\alpha\sqrt{\gamma}}{2} T^{\beta\nu}g_{\beta\nu, i}
 + \frac{\alpha\sqrt{\gamma}}{2} T^{\beta\nu}\left( g_{i\beta, \nu}
 - g_{i\nu, \beta} \right) \\

\partial_t \tilde{S}_i + \partial_j \left(\alpha\sqrt{\gamma} T_i^j \right)
&= \frac{\alpha\sqrt{\gamma}}{2} T^{\beta\nu}g_{\beta\nu,i}
\end{aligned}
\f]

---

\subsection energy_eq Evolution Equation for Energy Density Variable

The energy evolution equation comes from the same conservation equation as the
previous section, but with \f$ \mu=0 \f$. However, we choose to start with
\f$ \mu \f$ raised instead of lowered. Then,

\f[
\begin{aligned}
0 &= \nabla_\nu T^{\mu\nu} \\

&= \partial_\nu T^{\mu\nu}
 + \Gamma^\mu_{\sigma\nu} T^{\sigma\nu}
 + \Gamma^\nu_{\sigma\nu} T^{\mu\sigma} \\

&= \frac{1}{\sqrt{-g}}\partial_\nu \left(\sqrt{-g} T^{\mu\nu} \right)
 + \Gamma^\mu_{\sigma\nu} T^{\sigma\nu}
\end{aligned}
\f]

Setting \f$ \mu=0 \f$,

\f[
\begin{aligned}
\frac{1}{\sqrt{-g}}\partial_\nu \left(\sqrt{-g} T^{0\nu} \right)
&= -\Gamma^0_{\sigma\nu} T^{\sigma\nu} \\

\partial_t \left(\alpha\sqrt{\gamma} T^{00} \right) + \partial_i \left(\alpha\sqrt{\gamma} T^{0i} \right)
&= -\alpha\sqrt{\gamma}\Gamma^0_{\sigma\nu} T^{\sigma\nu} \\

&= -\frac{\alpha\sqrt{\gamma}}{2} T^{\sigma\nu} g^{0\beta}\left(g_{\beta\sigma,\nu}
                                                              + g_{\beta\nu,\sigma}
                                                              - g_{\sigma\nu,\beta}
                                                          \right) \\

&= -\frac{\alpha\sqrt{\gamma}}{2} T^{\sigma\nu} g^{0\beta}\left(2 g_{\beta\sigma,\nu}
                                                              - g_{\sigma\nu,\beta}
                                                          \right)
\end{aligned}
\f]

where the final step follows from the fact that \f$ T^{\sigma\nu} \f$ is
symmetric. Now, we will consider the right-hand side for various components of
\f$ T^{\sigma\nu} \f$. For these derivations, we will use derivations relating
to the extrinsic curvature in the [Appendix](@ref curv_id). For \f$ T^{00} \f$,

\f[
\begin{aligned}
-\frac{\alpha\sqrt{\gamma}}{2} T^{00} g^{0\beta}\left(2 g_{\beta 0,0}
                                                    - g_{00,\beta} \right)
&= - \frac{\alpha\sqrt{\gamma}}{2} T^{00} \left[2 \left(g^{00} g_{00,0}
                                                      + g^{0i}g_{0i,0}\right)
                                              - g^{00}g_{00,0}
                                              - g^{0i}g_{00,i}\right] \\

&= -\frac{\alpha\sqrt{\gamma}}{2} T^{00} \left(g^{00} g_{00,0}
                                             + 2g^{0i}g_{0i,0}
                                             - g^{0i}g_{00,i} \right) \\

&= -\frac{\alpha\sqrt{\gamma}}{2} T^{00} \left[-\alpha^{-2} \partial_t \left(\beta^2
                                                                           - \alpha^2
                                                                       \right)
                                               + \frac{2\beta^i}{\alpha^2} \partial_t \beta_i
                                               - \frac{\beta^i}{\alpha^2} \partial_i \left(\beta^2
                                                                                         - \alpha^2
                                                                                     \right)
                                         \right] \\

&= -\frac{\sqrt{\gamma}}{2\alpha} T^{00} \left[-\partial_t \left(\beta^2 - \alpha^2\right)
                                             + 2\beta^i \partial_t \beta_i
                                             - \beta^i \partial_i \left(\beta^2
                                                                      - \alpha^2\right)
                                         \right] \\

&= -\frac{\sqrt{\gamma}}{\alpha} T^{00} \left(\alpha\partial_t \alpha
                                            + \beta^i \alpha\partial_i \alpha
                                            - \beta^i \beta^j \partial_i \beta_j
                                         \right) \\

&= \sqrt{\gamma} T^{00} \left(\beta^i \beta^j K_{ij}
                            - \partial_t \alpha
                            - \beta^i \partial_i \alpha
                        \right)
\end{aligned}
\f]

where we have used the identity in the [Appendix](@ref curv_id). Next, we look
at the mixed term \f$ T^{0i} + T^{i0} \f$:

\f[
\begin{aligned}
-\frac{\alpha\sqrt{\gamma}}{2} T^{0i} g^{0\beta}\left(2 g_{\beta 0,i}
                                                    - g_{0i,\beta}
                                                    + 2 g_{\beta i,0}
                                                    - g_{i 0,\beta} \right)
&= -\alpha\sqrt{\gamma} T^{0i} \left[g^{00}\left( g_{00,i} + g_{0i,0} \right)
                                   + g^{0j}\left( g_{0j,i} + g_{ij,0} \right)
                                   - \left(g^{00} g_{0i,0} + g^{0j} g_{0i,j}\right)
                               \right] \\

&= -\alpha\sqrt{\gamma} T^{0i} \left[g^{00} g_{00,i}
                                   + g^{0j}\left( g_{0j,i} + g_{ij,0} - g_{0i,j} \right)
                               \right] \\

&= -\alpha\sqrt{\gamma} T^{0i} \left[-\alpha^{-2} \partial_i \left( \beta^2
                                                                  - \alpha^2 \right)
                                     + \frac{\beta^j}{\alpha^2} \left(\partial_i \beta_j
                                                                    + \partial_t \gamma_{ij}
                                                                    - \partial_j \beta_i
                                                                \right)
                               \right] \\

&= -\frac{\sqrt{\gamma}}{\alpha} T^{0i} \left[2\alpha\partial_i \alpha
                                            + \beta^j \left(\partial_t \gamma_{ij}
                                                          - \partial_i \beta_j
                                                          - \partial_j \beta_i
                                                      \right)
                                        \right] \\

&= -\frac{\sqrt{\gamma}}{\alpha} T^{0i} \left[2\alpha\partial_i \alpha
                                            + \beta^j \left(-2\alpha K_{ij}
                                                           - 2\Gamma^k_{ij}\beta_k
                                                      \right)
                                        \right] \\

&= 2\sqrt{\gamma} T^{0i} \left( \beta^j K_{ij} - \partial_i \alpha \right)
 + \frac{2\sqrt{\gamma}}{\alpha} T^{0i} \beta^j \Gamma^k_{ij}\beta_k \\

&= 2\sqrt{\gamma} T^{0i} \left( \beta^j K_{ij} - \partial_i \alpha \right)
\end{aligned}
\f]

Finally, for \f$ T^{ij} \f$ we have

\f[
\begin{aligned}
-\frac{\alpha\sqrt{\gamma}}{2} T^{ij} g^{0\beta}\left(2 g_{\beta i,j}
                                                    - g_{ij,\beta} \right)
&= -\frac{\alpha\sqrt{\gamma}}{2} T^{ij} \left[2\left( g^{00} g_{0i,j}
                                                     + g^{0k} g_{ki,j} \right)
                                             - g^{0\beta} g_{ij,\beta}
                                         \right] \\

&= -\frac{\alpha\sqrt{\gamma}}{2} T^{ij} \left[2\left(g^{00} \partial_j \beta_i
                                                    + g^{0k} \partial_j \gamma_{ki}
                                                \right)
                                             - g^{0\beta} \partial_\beta \gamma_{ij}
                                         \right] \\

&= -\frac{\alpha\sqrt{\gamma}}{2} T^{ij} \left[g^{00}\left(2\partial_j \beta_i
                                                         - \partial_t \gamma_{ij} \right)
                                             + g^{0k}\left(2\partial_j \gamma_{ki}
                                                         - \partial_k \gamma_{ij} \right)
                                         \right] \\

&= -\frac{\alpha\sqrt{\gamma}}{2} T^{ij} \left[-\alpha^{-2}\left(2\partial_j \beta_i
                                                               - \partial_t \gamma_{ij}
                                                           \right)
                                             + \frac{\beta^k}{\alpha^2}\left(2\partial_j \gamma_{ki}
                                                                           - \partial_k \gamma_{ij}
                                                                       \right)
                                         \right] \\
\end{aligned}
\f]

We can see that, thanks to the symmetry of \f$ T^{ij} \f$,

\f[
\begin{aligned}
T^{ij}\beta^k\left(2\partial_j \gamma_{ki} - \partial_k \gamma_{ij}\right)
&= T^{ij}\beta_n g^{kn}\left(\partial_j \gamma_{ki}
                           + \partial_i \gamma_{kj}
                           - \partial_k \gamma_{ij}\right) \\

&= 2T^{ij}\Gamma^k_{ij}\beta_k
\end{aligned}
\f]

Using this and using the definition of \f$ K_{ij} \f$ to replace
\f$ \partial_t \gamma_{ij} \f$,

\f[
\begin{aligned}
-\frac{\alpha\sqrt{\gamma}}{2} T^{ij} g^{0\beta}\left(2 g_{\beta i,j}
                                                    - g_{ij,\beta} \right)
&= \frac{\sqrt{\gamma}}{2\alpha} T^{ij} \left(2\partial_j \beta_i
                                            + 2\alpha K_{ij}
                                            - \partial_i\beta_j
                                            - \partial_j\beta_i
                                            + 2\Gamma^k_{ij}\beta_k
                                            - 2\Gamma^k_{ij}\beta_k \right) \\

&= \frac{\sqrt{\gamma}}{2\alpha} T^{ij} \left(2\alpha K_{ij}
                                            + \partial_j \beta_i
                                            - \partial_i\beta_j \right) \\

&= \sqrt{\gamma} T^{ij} K_{ij} \\
\end{aligned}
\f]

where the final step again follows from the symmetry of \f$ T^{ij} \f$. Putting
all this together and multiplying through by \f$ \alpha \f$,

\f[
\begin{aligned}
\alpha\partial_t \left(\alpha\sqrt{\gamma} T^{00} \right)
  + \alpha\partial_i \left(\alpha\sqrt{\gamma} T^{0i} \right)
&= \alpha\sqrt{\gamma} T^{00} \left(\beta^i \beta^j K_{ij}
                                  - \partial_t \alpha
                                  - \beta^i \partial_i \alpha \right) \\
&\hspace{5mm} + 2\alpha\sqrt{\gamma} T^{0i} \left(\beta^j K_{ij}
                                                - \partial_i \alpha \right)
              + \alpha\sqrt{\gamma} T^{ij} K_{ij} \\

\alpha\sqrt{\gamma} T^{00}\partial_t \alpha
  + 2\alpha\sqrt{\gamma} T^{0i}\partial_i \alpha
  + \alpha\partial_t \left(\alpha\sqrt{\gamma} T^{00} \right)
  + \alpha\partial_i \left(\alpha\sqrt{\gamma} T^{0i} \right)
&= \alpha\sqrt{\gamma} T^{00} \left(\beta^i \beta^j K_{ij}
                                  - \beta^i \partial_i \alpha \right)
 + 2\alpha\sqrt{\gamma} T^{0i}\beta^j K_{ij} \\
&\hspace{5mm} - \alpha\sqrt{\gamma} T^{0i} \partial_i \alpha
              + \alpha\sqrt{\gamma} T^{ij} K_{ij} \\

\partial_t \left(\alpha^2\sqrt{\gamma} T^{00} \right)
  + \partial_i \left(\alpha^2\sqrt{\gamma} T^{0i} \right)
&= \alpha\sqrt{\gamma} \left[\left(T^{00} \beta^i \beta^j
                                 + 2T^{0i}\beta^j
                                 + T^{ij} \right) K_{ij}
                           - \left(T^{00} \beta^i
                                 + T^{0i} \right)\partial_i \alpha
                       \right]
\end{aligned}
\f]

We have finally arrived at the evolution equation. However, we still need to
get it in terms of the evolution variable
\f$ \tilde{\tau} = \alpha^2 \sqrt{\gamma} T^{00} - \tilde{D} \f$. To do so, we
can add the fluid density evolution equation(\f$ =0 \f$) to the left-hand side. We also
define the source term

\f[
s = \alpha\sqrt{\gamma} \left[\left(T^{00} \beta^i \beta^j
                                  + 2T^{0i}\beta^j
                                  + T^{ij} \right) K_{ij}
                            - \left(T^{00} \beta^i
                                  + T^{0i} \right)\partial_i \alpha
                        \right]
\f]

Then, the energy evolution equation becomes

\f[
\begin{aligned}
\partial_t \left(\alpha^2\sqrt{\gamma} T^{00} \right)
 + \partial_i \left(\alpha^2\sqrt{\gamma} T^{0i} \right)
 - \partial_t \tilde{D} - \partial_i \left(\tilde{D} v^i \right) &= s \\

\partial_t \tilde{\tau}
 + \partial_i \left(\alpha^2\sqrt{\gamma} T^{0i}
 - \tilde{D} v^i \right) &= s
\end{aligned}
\f]

---

\subsection mag_eq Evolution Equation for the Magnetic Field

For evolution of magnetic fields, we elect to evolve the magnetic vector
potential \f$ A_i \f$ instead of directly evolving \f$ B^i \f$. This ensures
that the constraints are always satisfied (i.e. numerical errors will not
introduce magnetic monopoles). While unimportant for this derivation, using a
vector potential which is on a staggered grid with respect to the hydrodynamic
variables is vital for preserving the numerical accuracy on mesh boundaries in
a mesh refinement scheme (fixed or adaptive). \f$ A_i \f$ evolution can still
violate the constraint equations around mesh boundaries if \f$ A_i \f$ is
not edge-centered.

To find a convenient start for the derivation, I find the dual of Maxwell's
equation. Recall that the dual of the Faraday tensor is

\f[
F^{*\mu\nu} = \frac{1}{2} \epsilon^{\mu\nu\alpha\beta}F_{\alpha\beta}
\f]

Also, since the covariant derivative of the metric is 0,

\f[
\epsilon^{\alpha\lambda\mu\nu} F_{\mu\nu;\lambda}
 = \nabla_\lambda \left(\epsilon^{\alpha\lambda\mu\nu} F_{\mu\nu}\right)
 = 2\nabla_\lambda F^{*\alpha\lambda}
\f]

Therefore, Maxwell's equation becomes

\f[
\begin{aligned}
0 &= \epsilon^{\alpha\lambda\mu\nu} F_{[\mu\nu;\lambda]} \\

&= \epsilon^{\alpha\lambda\mu\nu} \left(F_{\mu\nu;\lambda}
                                      - F_{\lambda\nu;\mu}
                                      + F_{\nu\lambda;\mu}
                                      - F_{\mu\lambda;\nu}
                                      + F_{\lambda\mu;\nu}
                                      - F_{\nu\mu;\lambda} \right) \\

&= \epsilon^{\alpha\lambda\mu\nu}F_{\mu\nu;\lambda}
 + \epsilon^{\alpha\lambda\mu\nu} \nabla_\mu \left(F_{\nu\lambda}
                                                 - F_{\lambda\nu} \right)
 + \epsilon^{\alpha\lambda\mu\nu} \nabla_\nu \left(F_{\lambda\mu}
                                                 - F_{\mu\lambda}\right)
 - \epsilon^{\alpha\lambda\mu\nu}\nabla_\lambda F_{\nu\mu} \\

&= 2\nabla_\lambda F^{*\alpha\lambda}
 + 2\epsilon^{\alpha\lambda\mu\nu} \nabla_\mu F_{\nu\lambda}
 + 2\epsilon^{\alpha\lambda\mu\nu} \nabla_\nu F_{\lambda\mu}
 - \epsilon^{\alpha\lambda\mu\nu}\nabla_\lambda F_{\nu\mu}
\end{aligned}
\f]

By exploiting the symmetries of the Levi-Civita tensor,

\f[
\begin{aligned}
0 &= 2\nabla_\lambda F^{*\alpha\lambda}
   + 2\epsilon^{\alpha\mu\nu\lambda} \nabla_\mu F_{\nu\lambda}
   + 2\epsilon^{\alpha\nu\lambda\mu} \nabla_\nu F_{\lambda\mu}
   + \epsilon^{\alpha\lambda\nu\mu}\nabla_\lambda F_{\nu\mu} \\

&= 2\nabla_\lambda F^{*\alpha\lambda}
 + 4\nabla_\mu F^{*\alpha\mu}
 + 4\nabla_\nu F^{*\alpha\nu}
 + 2\nabla_\lambda F^{*\alpha\lambda} \\

&= 12\nabla_\lambda F^{*\alpha\lambda} \\

&= \nabla_\lambda F^{*\alpha\lambda} \\

&= \partial_\lambda F^{*\alpha\lambda}
 + \Gamma^\alpha_{\beta\lambda} F^{*\beta\lambda}
 + \Gamma^\lambda_{\beta\lambda} F^{*\alpha\beta} \\

&= \frac{1}{\sqrt{-g}}\partial_\lambda \left(\sqrt{-g} F^{*\alpha\lambda}\right)
\end{aligned}
\f]

where the first \f$ \Gamma \f$ term disappears because it is summing over the
multiplication of symmetric and anti-symmetric objects. We can now find the
magnetic equations from this. Taking the time component simply gives us the
no-monopole constraint. To see this, we need to first consider the individual
components of the dual. Recall that we defined \f$ B \f$ to be

\f[
B^\mu = \frac{1}{2\sqrt{4\pi}} \epsilon^{\mu\nu\beta\alpha}n_\nu F_{\alpha\beta}
 = \frac{1}{\sqrt{4\pi}} n_\nu F^{*\nu\mu}
\f]

Since \f$ B \f$ for the normal observer is purely spatial (\f$ B^\mu n_\mu=0 \f$),
this implies that \f$ F^{*00}=0 \f$. Then, the time component of the magnetic
equations is

\f[
\begin{aligned}
0 &= \frac{1}{\sqrt{-g}}\partial_\lambda \left(\sqrt{-g} F^{*0\lambda}\right) \\
&= \partial_i \left(\alpha\sqrt{\gamma} F^{*0i}\right) \\
&= \partial_i \left( \alpha\sqrt{\gamma} \frac{B^i}{\alpha}\right) \\
&= \partial_i \left(\tilde{B}^i \right)
\end{aligned}
\f]

Before examining the spatial components of Maxwell's equations, we will do
some preliminary work to make our lives easier later. In the co-moving frame,

\f[
\begin{aligned}
\frac{1}{\sqrt{4\pi}} F^{\mu\nu}
&= u^\mu E^\nu - u^\nu E^\mu + u_\gamma \epsilon^{\gamma\mu\nu\delta} B_{\delta} \\
&= u_\gamma \epsilon^{\gamma\mu\nu\delta} B_{\delta}
\end{aligned}
\f]

Then, the dual is

\f[
\begin{aligned}
\frac{1}{\sqrt{4\pi}} F^{*\mu\nu}
&= \frac{1}{2\sqrt{4\pi}} \epsilon^{\mu\nu\alpha\beta} F_{\alpha\beta} \\

&= \frac{1}{2} \epsilon^{\mu\nu\alpha\beta} u^\gamma
               \epsilon_{\gamma\alpha\beta\delta} B^{\delta} \\

&= \frac{1}{2} \epsilon^{\alpha\beta\mu\nu}\epsilon_{\alpha\beta\gamma\delta}
                                                                u^\gamma B^{\delta} \\

&= \left(\delta^\mu_\delta \delta^\nu_\gamma
       - \delta^\mu_\gamma \delta^\nu_\delta \right) u^\gamma B^{\delta} \\

&= u^\nu B^{\mu} - u^\mu B^{\nu}
\end{aligned}
\f]

where we have again used the Levi-Civita identity from the
[Appendix](@ref levi_civita). Next, we need to find a relationship between
the magnetic field of the normal observer and the co-moving observer. For
clarity, let the co-moving magnetic field be \f$ B^\mu_{(u)} \f$ and the normal
magnetic field be \f$ B^\mu \f$. Then, we can define a projection operator

\f[
P^{\mu\nu} = g^{\mu\nu} + u^\mu u^\nu
\f]

Naturally, the projection of the co-moving magnetic field should simply project
back into the same field:

\f[
\begin{aligned}
P^\mu_\nu B^\nu_{(u)} &= \left( \delta^\mu_\nu + u^\mu u_\nu \right)B^\nu_{(u)} \\
&= B^\mu_{(u)}
\end{aligned}
\f]

where we have used the orthogonality relation \f$ u_\nu B^\nu_{(u)}=0 \f$.
Projecting the normal observer's magnetic field,

\f[
\begin{aligned}
P^\mu_\nu B^\nu
&= \frac{1}{\sqrt{4\pi}} P^\mu_\nu n_\alpha F^{*\alpha\nu} \\

&= P^\mu_\nu n_\alpha \left(u^\nu B^{\alpha}_{(u)}
                          - u^\alpha B^{\nu}_{(u)} \right) \\

&= n_\alpha \left(\delta^\mu_\nu
                + u^\mu u_\nu \right) \left(u^\nu B^{\alpha}_{(u)}
                                          - u^\alpha B^{\nu}_{(u)} \right) \\

&= n_\alpha \left(u^\mu B^{\alpha}_{(u)}
                + u^\mu u_\nu u^\nu B^{\alpha}_{(u)}
                - u^\alpha B^{\mu}_{(u)}
                - u^\mu u_\nu u^\alpha B^{\nu}_{(u)} \right) \\

&= n_\alpha u^\mu B^{\alpha}_{(u)} \left(1 + u_\nu u^\nu \right)
 - n_\alpha u^\alpha \left( B^{\mu}_{(u)} + u^\mu u_\nu B^{\nu}_{(u)} \right) \\

&= - \alpha u^0 B^{\mu}_{(u)} \\

\Rightarrow B^{\mu}_{(u)} &= \frac{B^\mu + u^\mu u_\nu B^\nu}{\alpha u^0}
\end{aligned}
\f]

Since we will only need the spatial component for our purposes,

\f[
\begin{aligned}
B^{i}_{(u)} &= \frac{B^i + u^i u_\nu B^\nu}{\alpha u^0} \\
&= \frac{B^i + u^i u_j B^j}{\alpha u^0} \\
&= \frac{B^i}{\alpha u^0} + \frac{v^i u_j B^j}{\alpha}
\end{aligned}
\f]

where \f$ v^i \f$ is the conservative variable \f$ u^i/u^0 \f$. Finally, the
spatial components of the dual of Maxwell's equations gives the induction
equation for the magnetic field:

\f[
\begin{aligned}
0 &= \frac{1}{\sqrt{-g}}\partial_\nu \left(\sqrt{-g} F^{*i \nu}\right) \\
&= \partial_t \left(\alpha\sqrt{\gamma} F^{*i 0}\right)
 + \partial_j \left(\alpha\sqrt{\gamma} F^{*i j}\right) \\

&= \partial_t \left(\sqrt{\gamma} B^i \right)
 + \partial_j \left(\alpha\sqrt{\gamma} \left[u^j B^i_{(u)} - u^i B^j_{(u)}
              \right]\right) \\

&= \partial_t \tilde{B}^i
 + \partial_j \left(\alpha\sqrt{\gamma} \left[u^j \left(\frac{B^i}{\alpha u^0}
                                                      + \frac{v^i u_k B^k}{\alpha}
                                                  \right)
                                            - u^i \left(\frac{B^j}{\alpha u^0}
                                                      + \frac{v^j u_k B^k}{\alpha}
                                                  \right)
                                         \right]
              \right) \\

&= \partial_t \tilde{B}^i
 + \partial_j \left(v^j \tilde{B}^i
                  - v^i \tilde{B}^j
                  + \frac{\sqrt{\gamma}}{u^0}\left[u^j u^i u_k B^k
                                                 - u^i u^j u_k B^k \right]
              \right) \\

&= \partial_t \tilde{B}^i
 + \partial_j \left(v^j \tilde{B}^i - v^i \tilde{B}^j \right)
\end{aligned}
\f]

---

\subsection eq_summary Summary of the Conservative Variable Evolution Equations

In the previous sections, we have derived the evolution equations for the
conservative variables \f$ \mathbf{C} \f$. To summarize, these are

\f[
\begin{aligned}
\partial_t \tilde{D} + \partial_i \left(\tilde{D} v^i \right) &= 0 \\
\partial_t \tilde{\tau} + \partial_i \left(\alpha^2\sqrt{\gamma} T^{0i}
                                         - \tilde{D} v^i \right) &= s \\
\partial_t \tilde{S}_i + \partial_j \left(\alpha\sqrt{\gamma} T_i^j \right)
    &= \frac{\alpha\sqrt{\gamma}}{2} T^{\beta\nu}g_{\beta\nu,i} \\
\partial_t \tilde{B}^i
 + \partial_j \left(v^j \tilde{B}^i
 - v^i \tilde{B}^j \right) &= 0
\end{aligned}
\f]

where

\f[
s = \alpha\sqrt{\gamma} \left[\left(T^{00} \beta^i \beta^j
                                  + 2T^{0i}\beta^j
                                  + T^{ij}
                              \right) K_{ij}
                            - \left(T^{00} \beta^i
                                  + T^{0i}
                              \right)\partial_i \alpha
                        \right]
\f]

---

\section appendix Appendices

\subsection levi_civita Levi-Civita Contractions

In order to simplify the various combinations of \f$ F^{\mu\nu} \f$,
contractions of the Levi-Civita tensor are required. First, summing over the
first indices of two such tensors yields

\f[
\begin{aligned}
\epsilon^{\alpha\beta\sigma\mu} \epsilon_{\alpha\beta\gamma\delta}
&= -\delta^{\beta\sigma\mu}_{\hphantom{\beta\sigma\mu}\beta\gamma\delta} \\

&= - \delta^\beta_\beta \delta^{\sigma\mu}_{\gamma\delta}
   + \delta^\beta_\gamma \delta^{\sigma\mu}_{\beta\delta}
   - \delta^\beta_\delta \delta^{\sigma\mu}_{\beta\gamma} \\

&= - 4\left(\delta^\sigma_\gamma \delta^\mu_\delta
          - \delta^\sigma_\delta \delta^\mu_\gamma\right)
   + \delta^\beta_\gamma \left(\delta^\sigma_\beta \delta^\mu_\delta
                             - \delta^\sigma_\delta \delta^\mu_\beta\right)
   - \delta^\beta_\delta \left(\delta^\sigma_\beta \delta^\mu_\gamma
                             - \delta^\sigma_\gamma \delta^\mu_\beta\right) \\

&= - 4\left(\delta^\sigma_\gamma \delta^\mu_\delta
          - \delta^\sigma_\delta \delta^\mu_\gamma\right)
   + 2\left(\delta^\sigma_\gamma \delta^\mu_\delta
          - \delta^\sigma_\delta \delta^\mu_\gamma\right) \\

&= 2\left(\delta^\sigma_\delta \delta^\mu_\gamma
        - \delta^\sigma_\gamma \delta^\mu_\delta\right)
\end{aligned}
\f]

Second, summing over the third index (with all other components raised) yields

\f[
\begin{aligned}
\epsilon^{\gamma\mu\lambda\delta}
\epsilon^{\alpha\nu\hphantom{\lambda}\beta}_{\hphantom{\alpha\nu}\lambda}
&= g^{\alpha\tau}g^{\nu\phi}g^{\beta\theta}
   \epsilon^{\gamma\mu\lambda\delta} \epsilon_{\tau\phi\lambda\theta} \\

&= -g^{\alpha\tau} g^{\nu\phi} g^{\beta\theta}
   \delta^{\gamma\mu\delta}_{\tau\phi\theta} \\

&= -g^{\alpha\tau} g^{\nu\phi} g^{\beta\theta}
   \left(\delta^\gamma_\tau\delta^{\mu\delta}_{\phi\theta}
       - \delta^\gamma_\phi \delta^{\mu\delta}_{\tau\theta}
       + \delta^\gamma_\theta \delta^{\mu\delta}_{\tau\phi} \right) \\

&= -g^{\alpha\tau}g^{\nu\phi}g^{\beta\theta}
   \left(\delta^\gamma_\tau \left(\delta^\mu_\phi \delta^\delta_\theta
                                - \delta^\mu_\theta \delta^\delta_\phi \right)
       - \delta^\gamma_\phi \left(\delta^\mu_\tau \delta^\delta_\theta
                                - \delta^\mu_\theta \delta^\delta_\tau \right)
       + \delta^\gamma_\theta \left(\delta^\mu_\tau \delta^\delta_\phi
                                  - \delta^\mu_\phi \delta^\delta_\tau \right)
    \right) \\
\end{aligned}
\f]

---

\subsection curv_id Extrinsic Curvature Identities

The source term of the energy equation involves the extrinsic curvature.
Defined in terms of the 3+1 formalism quantities, the extrinsic curvature is

\f[
\begin{aligned}
K_{ij}
&= -\frac{1}{2\alpha}\left(\partial_t \gamma_{ij}
                         - \nabla_i\beta_j
                         - \nabla_j\beta_i \right) \\

&= -\frac{1}{2\alpha}\left(\partial_t \gamma_{ij}
                         - \partial_i\beta_j
                         - \partial_j\beta_i
                         + 2\Gamma^k_{ij}\beta_k \right)
\end{aligned}
\f]

As for how the extrinsic curvature enters into the energy equations, it
involves several relations between it and the other quantities which we must
show. First, consider the contraction of the Christoffel symbol with
\f$ \beta^i \beta^j \f$:

\f[
\begin{aligned}
\beta^j \Gamma^k_{ij}\beta_k
&= \beta^j \beta^k \left( \gamma_{k i,j} + \gamma_{k j,i} - \gamma_{i j,k} \right) \\

&= \beta^j \beta^k \gamma_{k j,i} \\

&= \alpha^2 \left( \gamma^{jk} - g^{jk} \right) \gamma_{k j,i} \\

&\propto \gamma^{jk}\partial_i \gamma_{k j} - g_{jk} \partial_i \gamma^{k j} \\

&\propto \gamma^{jk}\partial_i \gamma_{k j} - \gamma_{jk} \partial_i \gamma^{k j} \\

&\propto \gamma^{jk} \partial_i \gamma_{k j} - \gamma^{jk} \partial_i \gamma_{k j} \\

&= 0
\end{aligned}
\f]

The change \f$ g_{jk}\rightarrow\gamma_{jk} \f$ is allowed because the 3-metric
is identical to the 4-metric with only spatial indices. We can use this to find
a relationship between the partial derivative of \f$ \beta_j \f$ and the
extrinsic curvature:

\f[
\begin{aligned}
\beta^i \beta^j \partial_i \beta_j
&= \beta^i \beta^j \left(2\alpha K_{ij}
                       + \partial_t \gamma_{ij}
                       - \partial_j\beta_i
                       + 2\Gamma^k_{ij}\beta_k \right) \\

&= \beta^i \beta^j \left(2\alpha K_{ij}
                       + \partial_t \gamma_{ij}
                       - \partial_i\beta_j
                       + 2\Gamma^k_{ij}\beta_k \right) \\

&= \beta^i \beta^j \left(2\alpha K_{ij}
                       + \partial_t \gamma_{ij} \right)
 - \beta^i \beta^j\partial_i\beta_j \\

&= \beta^i \beta^j \left(\alpha K_{ij} + \frac{1}{2}\partial_t \gamma_{ij} \right) \\

&= \beta^i \beta^j \alpha K_{ij}
\end{aligned}
\f]

where in the second step we use the symmetry of \f$ \beta^i \beta^j \f$. The
final step uses the same trick as with \f$ \beta^j \Gamma^k_{ij}\beta_k \f$,
just with the derivative being \f$ \partial_t\f$ instead of \f$\partial_i \f$.
