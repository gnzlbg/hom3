HOM3
====

A laboratory for High-Order Multiphysics Multiscale Methods

## Right now

- Hierarchical Cartesian Grids
- Low-order finite volume methods
- Solver coupling (e.g. compressible fluid with heat conducting solid)
- VTK output

## Current solvers

- Low-order finite volume methods for:
  - the heat equation,
  - the Euler equations,
  - the Navier-Stokes equations.

## Solver coupling

Solving the shock-tube problem using two finite-volume solvers on two different domains: an exterior hollow cube and an interior cube.

<img src="docs/images/euler_euler_coupling.gif" align="center" width=75% />

## Short-term priorities

- boundary and volume coupling between different solvers
- parallelization of the grid generatior
- parallel I/O

## Acknowledgements

The financial support of the research cluster [Fuel production with renewable
raw materials](http://www.brenaro.rwth-aachen.de/pages/pageHomeDE.html)
(BrenaRo) at RWTH Aachen University is gratefully acknowledged.
