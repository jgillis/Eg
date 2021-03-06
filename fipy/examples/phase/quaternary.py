#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "quaternary.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 5/15/06 {3:11:15 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

r""" 
The same procedure used to construct the two-component phase field
diffusion problem in ``examples.phase.binary`` can be used to build up a
system of multiple components. Once again, we'll focus on 1D.

.. raw:: latex

   \IndexClass{Grid1D}

..

    >>> nx = 400
    >>> dx = 0.01
    >>> L = nx * dx
    >>> from fipy.meshes.grid1D import Grid1D
    >>> mesh = Grid1D(dx = dx, nx = nx)

.. raw:: latex

   We consider a free energy density \( f(\phi, C_0,\ldots,C_N, T) \) 
   that is a function of phase \( \phi \)
   \IndexClass{CellVariable}

..

    >>> from fipy.variables.cellVariable import CellVariable
    >>> phase = CellVariable(mesh=mesh, name='phase', value=1, hasOld=1)

.. raw:: latex

   interstitial components \( C_0 \ldots C_M \)
   
..

    >>> interstitials = [
    ...     CellVariable(mesh=mesh, name='C0', hasOld=1)
    ... ]

.. raw:: latex

   substitutional components \( C_{M+1} \ldots C_{N-1} \)
   
..

    >>> substitutionals = [
    ...     CellVariable(mesh=mesh, name='C1', hasOld=1),
    ...     CellVariable(mesh=mesh, name='C2', hasOld=1),
    ... ]

.. raw:: latex

   a ``solvent'' \( C_N \) that is constrained by the concentrations of the
   other substitutional species, such that \( C_N = 1 - \sum_{j=M}^{N-1}
   C_j \),

..

    >>> solvent = 1
    >>> for Cj in substitutionals:
    ...     solvent -= Cj
    >>> solvent.name = 'CN'


.. raw:: latex

   and temperature \( T \)
   
..

    >>> T = 1000

.. raw:: latex

   The free energy density of such a system can be written as
   \begin{align*}
       f(\phi, C_0, \ldots, C_N, T) 
       &= \sum_{j=0}^N C_j \left[ \mu^\circ_j(\phi, T) + R T \ln \frac{C_j}{\rho} \right]
   \end{align*}
   where
   
..

    >>> R = 8.314 # J / (mol K)

.. raw:: latex
    
   is the gas constant. As in the binary case, 
   \[
       \mu^\circ_j(\phi, T) = p(\phi) \mu_j^{\circ S}(T)
       + \left(1 - p(\phi)\right) \mu_j^{\circ L}(T) + \frac{W_j}{2} g(\phi)
   \]
   is constructed with the free energies of the pure components in each
   phase, given the ``tilting'' function
   
..

    >>> def p(phi):
    ...     return phi**3 * (6 * phi**2 - 15 * phi + 10)

and the "double well" function

    >>> def g(phi):
    ...     return (phi * (1 - phi))**2

.. raw:: latex

   We consider a very simplified model that has partial molar volumes \(
   \bar{V}_0 = \cdots = \bar{V}_{M} = 0 \) for the ``interstitials'' and \(
   \bar{V}_{M+1} = \cdots = \bar{V}_{N} = 1 \) for the ``substitutionals''.
   This approximation has been used in a number of models where density
   effects are ignored, including the treatment of electrons in
   electrodeposition processes \cite{ElPhFI,ElPhFII}. Under these
   constraints
   \begin{align*}
       \frac{\partial f}{\partial \phi}
       &= \sum_{j=0}^N C_j \frac{\partial f_j}{\partial \phi} 
       \nonumber \\
       &= \sum_{j=0}^N C_j \left[
           \mu_j^{\circ SL}(T) p'(\phi) + \frac{W_j}{2} g'(\phi)
       \right]
       \\
       \frac{\partial f}{\partial C_j}
       &= \left[\mu^\circ_j(\phi, T) + R T \ln \frac{C_j}{\rho} \right] 
       \nonumber \\
       &= \mu_j(\phi, C_j , T)
       \qquad\text{for \( j = 0\ldots M \)}
       \\
     \intertext{and}
       \frac{\partial f}{\partial C_j}
       &= \left[\mu^\circ_j(\phi, T) + R T \ln \frac{C_j}{\rho} \right] 
       - \left[\mu^\circ_N(\phi, T) + R T \ln \frac{C_N}{\rho} \right] 
       \nonumber \\
       &= \left[\mu_j(\phi, C_j, T) - \mu_N(\phi, C_N, T) \right]
       \qquad\text{for \( j = M+1\ldots N-1 \)}
   \end{align*}
   where \( \mu_j^{\circ SL}(T) \equiv \mu_j^{\circ S}(T) - \mu_j^{\circ
   L}(T) \) and where \( \mu_j \) is the classical chemical potential of
   component \( j \) for the binary species and \( \rho = 1 +
   \sum_{j=0}^{M} C_j \) is the total molar density.
   
..

    >>> rho = 1.
    >>> for Cj in interstitials:
    ...     rho += Cj
   
.. raw:: latex

   \( p'(\phi) \) and \( g'(\phi) \) are the partial derivatives of of \( p
   \) and \( g \) with respect to \( \phi \)

..    
    
    >>> def pPrime(phi):
    ...     return 30. * g(phi)

    >>> def gPrime(phi):
    ...     return 2. * phi * (1 - phi) * (1 - 2 * phi)

We "cook" the standard potentials to give the desired solid and liquid
concentrations, with a solid phase rich in interstitials and the solvent
and a liquid phase rich in the two substitutional species.
    
    >>> interstitials[0].S = 0.3
    >>> interstitials[0].L = 0.4
    >>> substitutionals[0].S = 0.4
    >>> substitutionals[0].L = 0.3
    >>> substitutionals[1].S = 0.2
    >>> substitutionals[1].L = 0.1
    >>> solvent.S = 1.
    >>> solvent.L = 1.
    >>> for Cj in substitutionals:
    ...     solvent.S -= Cj.S
    ...     solvent.L -= Cj.L

    >>> rhoS = rhoL = 1.
    >>> for Cj in interstitials:
    ...     rhoS += Cj.S
    ...     rhoL += Cj.L

.. raw:: latex

   \IndexModule{numerix}
   \IndexFunction{log}

..

    >>> from fipy.tools.numerix import log
    >>> for Cj in interstitials + substitutionals + [solvent]:
    ...     Cj.standardPotential = R * T * (log(Cj.L/rhoL) - log(Cj.S/rhoS))

    >>> for Cj in interstitials:
    ...     Cj.diffusivity = 1.
    ...     Cj.barrier = 0. 

    >>> for Cj in substitutionals:
    ...     Cj.diffusivity = 1.
    ...     Cj.barrier = R * T

    >>> solvent.barrier = R * T

-----

.. raw:: latex 

   We create the phase equation 
   \[
       \frac{1}{M_\phi}\frac{\partial \phi}{\partial t}
       = \kappa_\phi \nabla^2 \phi
       - \sum_{j=0}^N C_j \left[
                  \mu_j^{\circ SL}(T) p'(\phi) + \frac{W_j}{2} g'(\phi)
              \right]
   \]

with a semi-implicit source just as in ``examples.phase.simple.input`` and
``examples.phase.binary``

    >>> enthalpy = 0.
    >>> barrier = 0.
    >>> for Cj in interstitials + substitutionals + [solvent]:
    ...     enthalpy += Cj * Cj.standardPotential
    ...     barrier += Cj * Cj.barrier

    >>> mPhi = -((1 - 2 * phase) * barrier + 30 * phase * (1 - phase) * enthalpy)
    >>> dmPhidPhi = 2 * barrier - 30 * (1 - 2 * phase) * enthalpy
    >>> S1 = dmPhidPhi * phase * (1 - phase) + mPhi * (1 - 2 * phase)
    >>> S0 = mPhi * phase * (1 - phase) - S1 * phase * (S1 < 0)

.. raw:: latex

   \IndexClass{TransientTerm}
   \IndexClass{ImplicitDiffusionTerm}
   \IndexClass{ImplicitSourceTerm}

..
    
    >>> from fipy.terms.transientTerm import TransientTerm
    >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
    >>> from fipy.terms.implicitSourceTerm import ImplicitSourceTerm

    >>> phase.mobility = 1.
    >>> phase.gradientEnergy = 25
    >>> phase.equation = TransientTerm(coeff=1/phase.mobility) \
    ...   == ImplicitDiffusionTerm(coeff=phase.gradientEnergy) \
    ...      + S0 + ImplicitSourceTerm(coeff = S1 * (S1 < 0))

-----

We could construct the diffusion equations one-by-one, in the manner of
``examples.phase.binary``, but it is better to take advantage of the full
scripting power of the Python language, where we can easily loop over
components or even make "factory" functions if we desire. For the
interstitial diffusion equations, we arrange in canonical form as before:
    
.. raw:: latex

   \begin{align*}
       \underbrace{
           \frac{\partial C_j}{\partial t}
           \vphantom{\left\{
               \overbrace{
                   \left[
                       \mu_j^{\circ SL} \nabla p(\phi)
                   \right] 
               }^{\text{phase transformation}}
           \right\}}
       }_{\text{transient}}
       &= \underbrace{
           D_j\nabla^2 C_j 
           \vphantom{\left\{
               \overbrace{
                   \left[
                       \mu_j^{\circ SL} \nabla p(\phi)
                   \right] 
               }^{\text{phase transformation}}
           \right\}}
       }_{\text{diffusion}} \\
       & \qquad + \underbrace{
           D_j\nabla\cdot 
           \frac{C_j}{1 + \sum_{\substack{k=0\\ k \neq j}}^{M} C_k}
           \left\{
               \overbrace{
                   \frac{\rho}{R T}
                   \left[
                       \mu_j^{\circ SL} \nabla p(\phi)
                       + \frac{W_j}{2} \nabla g(\phi) 
                   \right] 
               }^{\text{phase transformation}}
               -
               \overbrace{
                   \sum_{\substack{i=0\\ i \neq j}}^{M} \nabla C_i
               }^{\text{counter diffusion}}
           \right\}
       }_{\text{convection}}
   \end{align*}
   \IndexClass{PowerLawConvectionTerm}

..

    >>> from fipy.terms.powerLawConvectionTerm import PowerLawConvectionTerm

    >>> for Cj in interstitials:
    ...     phaseTransformation = (rho.getHarmonicFaceValue() / (R * T)) \
    ...       * (Cj.standardPotential * p(phase).getFaceGrad() 
    ...          + 0.5 * Cj.barrier * g(phase).getFaceGrad())
    ...                            
    ...     CkSum = CellVariable(mesh=mesh, value=0.)
    ...     for Ck in [Ck for Ck in interstitials if Ck is not Cj]:
    ...         CkSum += Ck
    ...         
    ...     counterDiffusion = CkSum.getFaceGrad()
    ...     
    ...     convectionCoeff = counterDiffusion + phaseTransformation
    ...     convectionCoeff *= (Cj.diffusivity
    ...                         / (1. + CkSum.getHarmonicFaceValue()))
    ...                         
    ...     diffusionTerm = ImplicitDiffusionTerm(coeff=Cj.diffusivity)
    ...     convectionTerm = PowerLawConvectionTerm(coeff=convectionCoeff, 
    ...                                             diffusionTerm=diffusionTerm)
    ...                                             
    ...     Cj.equation = TransientTerm() == diffusionTerm + convectionTerm

-----

The canonical form of the substitutional diffusion equations is

.. raw:: latex

   \begin{align*}
      \underbrace{
           \frac{\partial C_j}{\partial t}
      }_{\text{transient}}
       &= \underbrace{
           D_{j}\nabla^2 C_j
           \vphantom{\frac{\partial C_j}{\partial t}}
       }_{\text{diffusion}} \\
       & \qquad + \underbrace{
           D_{j}\nabla\cdot 
           \frac{C_j}{1 - \sum_{\substack{k=2\\ k \neq j}}^{n-1} C_k}
           \left\{
              \overbrace{ 
                   \frac{C_N}{R T}
                   \left[
                       \left(\mu_j^{\circ SL} - \mu_N^{\circ SL}\right) \nabla p(\phi)
                       + \frac{W_j - W_N}{2} \nabla g(\phi) 
                   \right] 
                   \vphantom{\sum_{\substack{i=M+1\\ i \neq j}}^{N} \nabla C_i}
              }^{\text{phase transformation}}
              +
              \overbrace{
                  \sum_{\substack{i=M+1\\ i \neq j}}^{N} \nabla C_i
              }^{\text{counter diffusion}}
           \right\}
       }_{\text{convection}}
   \end{align*}

..
    
    >>> for Cj in substitutionals:
    ...     phaseTransformation = (solvent.getHarmonicFaceValue() / (R * T)) \
    ...       * ((Cj.standardPotential - solvent.standardPotential) * p(phase).getFaceGrad() 
    ...          + 0.5 * (Cj.barrier - solvent.barrier) * g(phase).getFaceGrad())
    ...                            
    ...     CkSum = CellVariable(mesh=mesh, value=0.)
    ...     for Ck in [Ck for Ck in substitutionals if Ck is not Cj]:
    ...         CkSum += Ck
    ...         
    ...     counterDiffusion = CkSum.getFaceGrad()
    ...     
    ...     convectionCoeff = counterDiffusion + phaseTransformation
    ...     convectionCoeff *= (Cj.diffusivity
    ...                         / (1. - CkSum.getHarmonicFaceValue()))
    ...                         
    ...     diffusionTerm = ImplicitDiffusionTerm(coeff=Cj.diffusivity)
    ...     convectionTerm = PowerLawConvectionTerm(coeff=convectionCoeff, 
    ...                                             diffusionTerm=diffusionTerm)
    ...                                             
    ...     Cj.equation = TransientTerm() == diffusionTerm + convectionTerm

-----

We start with a sharp phase boundary

.. raw:: latex

   $$ \xi =
   \begin{cases}
       1& \text{for $x \le L/2$,} \\
       0& \text{for $x > L/2$,}
   \end{cases} $$

..

    >>> x = mesh.getCellCenters()[...,0]
    >>> phase.setValue(1.)
    >>> phase.setValue(0., where=x > L / 2)

and with uniform concentration fields, initially equal to the average of
the solidus and liquidus concentrations

    >>> for Cj in interstitials + substitutionals:
    ...     Cj.setValue((Cj.S + Cj.L) / 2.)

If we're running interactively, we create a viewer
    
.. raw:: latex

   \IndexModule{viewers}

..

    >>> if __name__ == '__main__':
    ...     from fipy import viewers
    ...     viewer = viewers.make(vars = [phase] \
    ...                                  + interstitials + substitutionals \
    ...                                  + [solvent],
    ...                           limits = {'datamin': 0, 'datamax': 1})
    ...     viewer.plot()

and again iterate to equilibrium

.. raw:: latex

   \IndexClass{LinearLUSolver}

..

    >>> from fipy.solvers.linearLUSolver import LinearLUSolver
    >>> solver = LinearLUSolver(tolerance=1e-3)

    >>> dt = 10000
    >>> for i in range(5):
    ...     for field in [phase] + substitutionals + interstitials:
    ...         field.updateOld()
    ...     phase.equation.solve(var = phase, dt = dt)
    ...     for field in substitutionals + interstitials:
    ...         field.equation.solve(var = field, 
    ...                              dt = dt,
    ...                              solver = solver)
    ...     if __name__ == '__main__':    
    ...         viewer.plot()

.. image:: examples/phase/quaternary.pdf
   :scale: 50
   :align: center

..

We can confirm that the far-field phases have remained separated

.. raw:: latex

   \IndexModule{numerix}
   \IndexFunction{take}
   \IndexFunction{allclose}

..

    >>> from fipy.tools.numerix import take, allclose
    >>> ends = take(phase, (0,-1))
    >>> allclose(ends, (1.0, 0.0), rtol = 1e-5, atol = 1e-5)
    1
    
and that the concentration fields have appropriately segregated into 
their equilibrium values in each phase

    >>> equilibrium = True
    >>> for Cj in interstitials + substitutionals:
    ...     ends = take(Cj, (0,-1))
    ...     equilibrium &= allclose(ends, (Cj.S, Cj.L), rtol = 3e-3, atol = 3e-3)
    >>> print equilibrium
    1
"""
__docformat__ = 'restructuredtext'

if __name__ == "__main__": 
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    raw_input("finished")

## if __name__ == '__main__':
##     ## from fipy.tools.profiler.profiler import Profiler
##     ## from fipy.tools.profiler.profiler import calibrate_profiler
## 
##     # fudge = calibrate_profiler(10000)
##     # profile = Profiler('profile', fudge=fudge)
## 
##     import fipy.tests.doctestPlus
##     exec(fipy.tests.doctestPlus.getScript())
## 
##     # profile.stop()
## 	    
##     raw_input("finished")
## 
