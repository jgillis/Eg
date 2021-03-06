#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "inputSimpleTrenchSystem.py"
 #                                    created: 8/26/04 {10:29:10 AM} 
 #                                last update: 5/18/06 {8:41:56 PM} { 1:23:41 PM}
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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
This input file

.. raw:: latex

    \label{inputSimpleTrench} is a demonstration of the use of \FiPy{}
    for modeling electrodeposition using the CEAC mechanism. The
    material properties and experimental parameters used are roughly
    those that have been previously
    published~\cite{NIST:damascene:2003}.

To run this example from the base fipy directory type::
    
    $ examples/levelSet/electroChem/inputSimpleTrenchSystem.py

at the command line. The results of the simulation will be displayed
and the word `finished` in the terminal at the end of the
simulation. In order to alter the number of timesteps, the python
function that encapsulates the system of equations must first be
imported (at the python command line),

.. raw:: latex

   \IndexFunction{runSimpleTrenchSystem}

..

    >>> from examples.levelSet.electroChem.inputSimpleTrenchSystem import runSimpleTrenchSystem

and then the function can be run with a different number of time steps
with the `numberOfSteps` argument as follows,

    >>> runSimpleTrenchSystem(numberOfSteps=5, displayViewers=False)
    1

Change the `displayViewers` argument to `True` if you wish to see the

.. raw:: latex

    results displayed on the
    screen. Example~\ref{inputWriteScriptHowTo} gives explanation for
    writing new scripts or modifying existing scripts that are
    encapsulated by functions.

Any argument parameter can be changed. For example if the initial
catalyst coverage is not 0, then it can be reset,

    >>> runSimpleTrenchSystem(catalystCoverage=0.1, displayViewers=False)
    0

The following image shows a schematic of a trench geometry along with
the governing equations for modeling electrodeposition with the CEAC
mechanism. All of the given equations are implemented in the
`runSimpleTrenchSystem` function. As stated above, all the parameters
in the equations can be changed with function arguments.

.. image:: examples/levelSet/electroChem/schematicOfEquations.pdf
   :scale: 60
   :align: center
   :alt: schematic of superfill equations

The following table shows the symbols used in the governing equations
and their corresponding arguments to the `runSimpleTrenchSystem`
function. The boundary layer depth is intentionally small in this
example in order not to complicate the mesh. Further examples will
simulate more realistic boundary layer depths but will also have more
complex meshes requiring the `gmsh` software.

.. raw:: latex

    \IndexSoftware{gmsh}
    %
    \begin{tabular}{|rllr@{.}ll|}
    \hline
    Symbol                & Description                       & Keyword Argument                      & \multicolumn{2}{l}{Value} & Unit                               \\
    \hline
    \multicolumn{6}{|c|}{Deposition Rate Parameters}                                                                                                                   \\
    \hline
    $v$                   & deposition rate                   &                                       & \multicolumn{2}{l}{}      & m s$^{-1}$                         \\
    $i$                   & current density                   &                                       & \multicolumn{2}{l}{}      & A m$^{-2}$                         \\
    $\Omega$              & molar volume                      & \verb+molarVolume+                    & 7&1$\times$10$^{-6}$      & m$^3$ mol$^{-1}$                   \\
    $n$                   & ion charge                        & \verb+charge+                         & \multicolumn{2}{c}{2}     &                                    \\
    $F$                   & Faraday's constant                & \verb+faradaysConstant+               & 9&6$\times$10$^{-4}$      & C mol$^{-1}$                       \\
    $i_0$                 & exchange current density          &                                       & \multicolumn{2}{l}{}      & A m$^{-2}$                         \\
    $\alpha$              & transfer coefficient              & \verb+transferCoefficient+            & 0&5                       &                                    \\
    $\eta$                & overpotential                     & \verb+overpotential+                  & -0&3                      & V                                  \\
    $R$                   & gas constant                      & \verb+gasConstant+                    & 8&314                     & J K$^{-1}$ mol$^{-1}$               \\
    $T$                   & temperature                       & \verb+temperature+                    & 298&0                     & K                                  \\
    $b_0$                 & current density for $\theta^0$    & \verb+currentDensity0+                & 0&26                      & A m$^{-2}$                         \\
    $b_1$                 & current density for $\theta$      & \verb+currentDensity1+                & 45&0                      & A m$^{-2}$                         \\
    \hline
    \multicolumn{6}{|c|}{Metal Ion Parameters}                                                                                                                         \\
    \hline
    $c_m$                 & metal ion concentration           & \verb+metalConcentration+             & 250&0                     & mol m$^{-3}$                       \\
    $c_m^{\infty}$        & far field metal ion concentration & \verb+metalConcentration+             & 250&0                     & mol m$^{-3}$                       \\
    $D_m$                 & metal ion diffusion coefficient   & \verb+metalDiffusion+                 & 5&6$\times$10$^{-10}$     & m$^2$ s$^{-1}$                     \\
    \hline
    \multicolumn{6}{|c|}{Catalyst Parameters}                                                                                                                          \\
    \hline
    $\theta$              & catalyst surfactant concentration & \verb+catalystCoverage+               & 0&0                       &                                    \\
    $c_{\theta}$          & bulk catalyst concentration       & \verb+catalystConcentration+          & 5&0$\times$10$^{-3}$      & mol m$^{-3}$                       \\
    $c_{\theta}^{\infty}$ & far field catalyst concentration  & \verb+catalystConcentration+          & 5&0$\times$10$^{-3}$      & mol m$^{-3}$                       \\
    $D_{\theta}$          & catalyst diffusion coefficient    & \verb+catalystDiffusion+              & 1&0$\times$10$^{-9}$      & m$^2$ s$^{-1}$                     \\
    $\Gamma$              & catalyst site density             & \verb+siteDensity+                    & 9&8$\times$10$^{-6}$      & mol m$^{-2}$                       \\
    $k$                   & rate constant                     &                                       & \multicolumn{2}{l}{}      & m$^3$ mol$^{-1}$ s$^{-1}$          \\
    $k_0$                 & rate constant for $\eta^0$        & \verb+rateConstant0+                  & 1&76                      & m$^3$ mol$^{-1}$ s$^{-1}$          \\
    $k_3$                 & rate constant for $\eta^3$        & \verb+rateConstant3+                  & -245&0$\times$10$^{-6}$   & m$^3$ mol$^{-1}$ s$^{-1}$ V$^{-3}$ \\
    \hline
    \multicolumn{6}{|c|}{Geometry Parameters}                                                                                                                          \\
    \hline
    $D$                   & trench depth                      & \verb+trenchDepth+                    & 0&5$\times$10$^{-6}$      & m                                  \\
    $D / W$               & trench aspect ratio               & \verb+aspectRatio+                    & 2&0                       &                                    \\
    $S$                   & trench spacing                    & \verb+trenchSpacing+                  & 0&6$\times$10$^{-6}$      & m                                  \\
    $\delta$              & boundary layer depth              & \verb+boundaryLayerDepth+             & 0&3$\times$10$^{-6}$      & m                                  \\
    \hline
    \multicolumn{6}{|c|}{Simulation Control Parameters}                                                                                                                \\
    \hline
                          & computational cell size           & \verb+cellSize+                       & 0&1$\times$10$^{-7}$      & m                                  \\
                          & number of time steps              & \verb+numberOfSteps+                  & \multicolumn{2}{c}{5}     &                                    \\
                          & whether to display the viewers    & \verb+displayViewers+                 & \multicolumn{2}{c}{\texttt{True}} &                           \\
    \hline
    \end{tabular}

If the MayaVi plotting
software is

.. raw:: latex

    installed (see Chapter~\ref{chap:Installation}) then a plot should
    appear that is updated every 20 time steps and will eventually

resemble the image below.

.. image:: examples/levelSet/electroChem/inputSimpleTrenchSystem.pdf
   :scale: 60
   :align: center
   :alt: resulting image

"""
__docformat__ = 'restructuredtext'

def runSimpleTrenchSystem(faradaysConstant=9.6e4,
                          gasConstant=8.314,
                          transferCoefficient=0.5,
                          rateConstant0=1.76,
                          rateConstant3=-245e-6,
                          catalystDiffusion=1e-9,
                          siteDensity=9.8e-6,
                          molarVolume=7.1e-6,
                          charge=2,
                          metalDiffusion=5.6e-10,
                          temperature=298.,
                          overpotential=-0.3,
                          metalConcentration=250.,
                          catalystConcentration=5e-3,
                          catalystCoverage=0.,
                          currentDensity0=0.26,
                          currentDensity1=45.,
                          cellSize=0.1e-7,
                          trenchDepth=0.5e-6,
                          aspectRatio=2.,
                          trenchSpacing=0.6e-6,
                          boundaryLayerDepth=0.3e-6,
                          numberOfSteps=5,
                          displayViewers=True):

    cflNumber = 0.2
    numberOfCellsInNarrowBand = 10
    cellsBelowTrench = 10
    
    yCells = cellsBelowTrench \
             + int((trenchDepth + boundaryLayerDepth) / cellSize)

    xCells = int(trenchSpacing / 2 / cellSize)

    from fipy.meshes.grid2D import Grid2D
    mesh = Grid2D(dx = cellSize,
                  dy = cellSize,
                  nx = xCells,
                  ny = yCells)

    narrowBandWidth = numberOfCellsInNarrowBand * cellSize
    from fipy.models.levelSet.distanceFunction.distanceVariable import \
         DistanceVariable        

    distanceVar = DistanceVariable(
        name = 'distance variable',
        mesh = mesh,
        value = -1,
        narrowBandWidth = narrowBandWidth,
        hasOld = 1)

    bottomHeight = cellsBelowTrench * cellSize
    trenchHeight = bottomHeight + trenchDepth
    trenchWidth = trenchDepth / aspectRatio
    sideWidth = (trenchSpacing - trenchWidth) / 2

    x, y = mesh.getCellCenters()[...,0], mesh.getCellCenters()[...,1]
    distanceVar.setValue(1, where=(y > trenchHeight) | ((y > bottomHeight) & (x < xCells * cellSize - sideWidth)))

    distanceVar.calcDistanceFunction(narrowBandWidth = 1e10)

    from fipy.models.levelSet.surfactant.surfactantVariable import \
         SurfactantVariable
    
    catalystVar = SurfactantVariable(
        name = "catalyst variable",
        value = catalystCoverage,
        distanceVar = distanceVar)
    
    from fipy.variables.cellVariable import CellVariable

    bulkCatalystVar = CellVariable(
        name = 'bulk catalyst variable',
        mesh = mesh,
        value = catalystConcentration)

    metalVar = CellVariable(
        name = 'metal variable',
        mesh = mesh,
        value = metalConcentration)
    
    expoConstant = -transferCoefficient * faradaysConstant \
                   / (gasConstant * temperature)
    
    tmp = currentDensity1 * catalystVar.getInterfaceVar()

    exchangeCurrentDensity = currentDensity0 + tmp

    import fipy.tools.numerix as numerix
    expo = numerix.exp(expoConstant * overpotential)
    currentDensity = expo * exchangeCurrentDensity * metalVar \
                     / metalConcentration

    depositionRateVariable = currentDensity * molarVolume \
                             / (charge * faradaysConstant)

    extensionVelocityVariable = CellVariable(
        name = 'extension velocity',
        mesh = mesh,
        value = depositionRateVariable)   

    from fipy.models.levelSet.surfactant.adsorbingSurfactantEquation \
                import AdsorbingSurfactantEquation

    surfactantEquation = AdsorbingSurfactantEquation(
        surfactantVar = catalystVar,
        distanceVar = distanceVar,
        bulkVar = bulkCatalystVar,
        rateConstant = rateConstant0 + rateConstant3 * overpotential**3)

    from fipy.models.levelSet.advection.higherOrderAdvectionEquation \
                   import buildHigherOrderAdvectionEquation

    advectionEquation = buildHigherOrderAdvectionEquation(
        advectionCoeff = extensionVelocityVariable)

    from fipy.boundaryConditions.fixedValue import FixedValue
    from fipy.models.levelSet.electroChem.metalIonDiffusionEquation \
                         import buildMetalIonDiffusionEquation

    metalEquation = buildMetalIonDiffusionEquation(
        ionVar = metalVar,
        distanceVar = distanceVar,
        depositionRate = depositionRateVariable,
        diffusionCoeff = metalDiffusion,
        metalIonMolarVolume = molarVolume,
    )

    metalEquationBCs = FixedValue(mesh.getFacesTop(), metalConcentration)

    from fipy.models.levelSet.surfactant.surfactantBulkDiffusionEquation \
                    import buildSurfactantBulkDiffusionEquation

    bulkCatalystEquation = buildSurfactantBulkDiffusionEquation(
        bulkVar = bulkCatalystVar,
        distanceVar = distanceVar,
        surfactantVar = catalystVar,
        diffusionCoeff = catalystDiffusion,
        rateConstant = rateConstant0 * siteDensity
    )

    catalystBCs = FixedValue(mesh.getFacesTop(), catalystConcentration)

    if displayViewers:
        try:
            from fipy.viewers.mayaviViewer.mayaviSurfactantViewer import MayaviSurfactantViewer
            viewers = (MayaviSurfactantViewer(distanceVar, catalystVar.getInterfaceVar(), zoomFactor = 1e6, limits = { 'datamax' : 0.5, 'datamin' : 0.0 }, smooth = 1, title = 'catalyst coverage'),)
        except:
            from fipy.viewers import make
            viewers = (
                make(distanceVar, limits = { 'datamin' :-1e-9 , 'datamax' : 1e-9 }),
                make(catalystVar.getInterfaceVar()))
    else:
        viewers = ()

    levelSetUpdateFrequency = int(0.8 * narrowBandWidth \
                                  / (cellSize * cflNumber * 2))

    for step in range(numberOfSteps):

        if step % 5 == 0:
            for viewer in viewers:
                viewer.plot()

        if step % levelSetUpdateFrequency == 0:
            distanceVar.calcDistanceFunction()
            
        extensionVelocityVariable.setValue(depositionRateVariable())

        distanceVar.updateOld()
        catalystVar.updateOld()
        metalVar.updateOld()
        bulkCatalystVar.updateOld()

        distanceVar.extendVariable(extensionVelocityVariable)
        dt = cflNumber * cellSize / numerix.max(extensionVelocityVariable)

        advectionEquation.solve(distanceVar, dt = dt) 
        surfactantEquation.solve(catalystVar, dt = dt)
        metalEquation.solve(metalVar, dt = dt, 
                            boundaryConditions = metalEquationBCs)
        bulkCatalystEquation.solve(bulkCatalystVar, dt = dt,
                                   boundaryConditions = catalystBCs)

    try:
        import os
        import examples.levelSet.electroChem
        filepath = os.path.join(examples.levelSet.electroChem.__path__[0], 'test.gz')
        
        from fipy.tools import dump
        from fipy.tools import numerix
        print catalystVar.allclose(numerix.array(dump.read(filepath)), rtol = 1e-4)
    except:
        return 0

if __name__ == '__main__':
    runSimpleTrenchSystem(numberOfSteps = 800, cellSize = 0.05e-7)
    raw_input("finished")
