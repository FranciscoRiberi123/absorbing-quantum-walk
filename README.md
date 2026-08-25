# Continuous-time quantum walks with a Lindblad boundary sink

This repository contains the Wolfram Language / Mathematica code accompanying the paper

**"Boundary-localized non-Hermitian modes in an absorbing quantum walk"**  
Francisco Riberi

The code reproduces the numerical figures in the manuscript and includes
the numerical convergence checks discussed in the referee response.

## Contents

- `absorbing_QW_Fig2_survival_first_passage_with_convergence.m`  
  Generates the survival-probability and first-passage-density curves.
  Includes checks of probability conservation, finite-chain convergence,
  population near the artificial boundary, and the ballistic-front estimate.

- `absorbing_QW_Wigner_with_convergence.m`  
  Generates the physical and conditionally normalized doubled-lattice
  Wigner-function figures, including the boundary-localized pole contribution.
  Includes finite-chain, boundary-population, ballistic-front, and
  momentum-grid convergence checks.

## Numerical method

The semi-infinite lattice is approximated by a finite chain, and the state
is propagated directly under the finite-dimensional non-Hermitian effective
Hamiltonian. 

The scripts contain explicit convergence tests for the finite-chain cutoff.
For the Wigner-function figures, the momentum grid is also tested for
resolution convergence.

## Requirements

The calculations were performed using Wolfram Mathematica 14.3.

No external Mathematica packages are required.

## Usage

Each script is self-contained and can be evaluated from a fresh Mathematica
kernel. The convergence checks are included after the figure-generation
sections and print their results directly when evaluated.

## Paper

https://arxiv.org/abs/2605.08056

## License

This code is distributed under the MIT License.
