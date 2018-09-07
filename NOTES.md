# NBodyIPFitting Development Notes

## List of variations to explore

* Laplace Regulariser
   * Integrate w.r.t. r or x
   * Apply Laplace w.r.t. r, x, or R?
* How to incorporate environment
   * next logical step: simple sums of pair potentials for the environment representation
   * then bond-angles?
* "Correct" weights for forces and virials
* Choice of space transform
   * speed vs accuracy, polynmomials vs exponentials
   * possibly optimise the meta-parameters e.g. the a in exp( - a k r )
* high-pressure data and match with 2B for core (forgot the name)
* different choices of the angle variable, angle, cos, or just ip (current)
* BL vs BA
* Is there some way to use orthogonal invariant polynomials?
