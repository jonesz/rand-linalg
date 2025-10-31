# Random linear algebra package for Futhark

Routines for randomized numerical linear algebra; sketches, randomized rangefinder, rsvd.

## Installation

```
$ futhark pkg add github.com/jonesz/rand-linalg
$ futhark pkg sync
```

## Examples

#### Compute the SVD of a `[m][m]` matrix (without randomness).

```futhark
import "lib/github.com/jonesz/rand-linalg/svd"
module SVD = mk_one_sided_jacobi_serial f32

let A = tabulate_2d 10 10 (\a b -> a + b |> f32.i64)
let (U, S, V_T) = SVD.svd A

-- See that A = U * S * V_T

import "lib/github.com/diku-dk/linalg/linalg"
module L = mk_linalg f32

let A_reconstructed = L.matmul U S |> flip (L.matmul) V_T

```

A matrix that is non-square can be handled by computing the QR/LQ decomposition,
then computing the SVD of the square {R, L}.

#### Compute the RSVD of a `[m][n]` matrix where `m >> n`

```futhark
import "lib/github.com/jonesz/cbrng-fut/cbrng"
import "lib/github.com/jonesz/rand-linalg/sketch"
import "lib/github.com/jonesz/rand-linalg/rangefinder"
import "lib/github.com/jonesz/rand-linalg/svd"

-- Choose the sketch matrix.
module SK = mk_gaussian_embedding f32 u32 squares32
-- Choose the randomized rangefinder.
module RF = mk_rangefinder f32 SK.dense.right
module SVD = mk_rsvd f32 RF

let A = tabulate_2d 1000 20 (\a b -> a + b |> f32.i64)
let seed = SVD.dist.engine.construct 1337i64
let (U, S, V_T) = SVD.rsvd seed A 15

import "lib/github.com/diku-dk/linalg/linalg"
module L = mk_linalg f32

let A_reconstructed = L.matmul U S |> flip (L.matmul) V_T
```
The target rank should be `l <= n`.


## References

Significant code is based off of Tropp's Winter 2020 [*Randomized Algorithms for Matrix Computations*](https://tropp.caltech.edu/notes/Tro20-Randomized-Algorithms-LN.pdf).
