# Cellular-Potts

The basic idea here is that we have a 2D grid of collagen fibers, on top of which we put a number of cells.

The program keeps track of this 2D lattice, with the value at each lattice point representing its state:
  air
  cell #1
  cell #1 perimeter
  cell #2
  cell #2 perimeter
  etc.
  
We loop through time, and at each time point try to cleverly vary the state of a certain lattice point, keeping cells intact.

We define an energy of the system based on cell shape and collagen overlap, and accept or reject variations with via Metropolis.

In the end, we spit out a bunch of text files that can be read into Mathematica to make movies of what happened.
