# Random linear algebra package for Futhark

Routines for randomized numerical linear algebra; sketches, randomized rangefinder, rsvd.

## Installation

```
$ futhark pkg add github.com/jonesz/rand-linalg
$ futhark pkg sync
```

## Examples

#### Compute the SVD of an `[m][m]` matrix (without randomness).


```futhark
import "lib/github.com/jonesz/rand-linalg/svd"
module OSJS = mk_one_sided_jacobi_slow f32

let A: [3][3]f32 = [[1, 2, -1], [3, -4, -5], [6, 7, 100]]
let (U, S, V_T) = OSJS.svd A

-- ([[-0.20220032, 0.9793124, -7.906095e-3],
--   [0.978168, 0.20155467, -5.0621767e-2],
--   [4.7980998e-2, 1.7969193e-2, 0.99868655]],
--  [[4.998353, 0.0, 0.0],
--   [0.0, 2.1288712, 0.0],
--   [0.0, 0.0, 100.555885]],
--  [[0.6042373, -0.79650337, 2.1899257e-2],
--   [0.79469055, 0.60040605, -8.932731e-2],
--   [5.800106e-2, 7.137804e-2, 0.9957616]])

import "lib/github.com/diku-dk/linalg/linalg"
module LA = mk_linalg f32

-- See that A = U * S * V_T
let A_reconstructed = L.matmul U S |> flip (LA.matmul) V_T

--  [[1.0, 1.9999996, -0.99999994],
--  [2.9999998, -3.9999993, -5.0],
--  [5.9999995, 7.0000005, 100.000015]]

```

A matrix that is non-square can be handled by computing the QR/LQ decomposition,
then computing the SVD of the square R/L.

#### Compute the RSVD of a `[m][n]` matrix where `m >> n`

```futhark
import "lib/github.com/jonesz/rand-linalg/svd"
module RSVD = mk_rsvd_default f32

-- See the tests in `test/test_svd.fut` for a matrix that isn't structured.
let A = tabulate_2d 1000 20 (\a b -> a + b |> f32.i64)
let seed = RSVD.dist.engine.construct 1337i64
let (U, S, V_T) = SVD.rsvd seed A 15

import "lib/github.com/diku-dk/linalg/linalg"
module L = mk_linalg f32

let A_reconstructed = L.matmul U S |> flip (L.matmul) V_T
```

The target rank should be `l <= n`.

## References

Significant code is based off of Tropp's Winter 2020 [*Randomized Algorithms for Matrix Computations*](https://tropp.caltech.edu/notes/Tro20-Randomized-Algorithms-LN.pdf).
