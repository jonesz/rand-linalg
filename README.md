# Random linear algebra package for Futhark

Routines for randomized numerical linear algebra.

## Installation

```
$ futhark pkg add github.com/jonesz/rand-linalg
$ futhark pkg sync
```

## References

Significant code is based off of Tropp's Winter 2020
 [*Randomized Algorithms for Matrix Computations*](https://tropp.caltech.edu/notes/Tro20-Randomized-Algorithms-LN.pdf).

| Algorithm | Source |
| --------- | ------ |
| Hutch++ | [Hutch++: Optimal Stochastic Trace Estimation (2021)](https://arxiv.org/pdf/2010.09649)
| TSQR    | [Communication-avoiding parallel and sequential QR factorizations](https://bebop.cs.berkeley.edu/pubs/mhoemmen2008-tsqr-tech-report.pdf)
