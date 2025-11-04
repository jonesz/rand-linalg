-- | Sketching operations which transform a matrix from the *ambient*
-- dimension to the *embedding* dimension.

import "../cbrng-fut/cbrng"
import "../cbrng-fut/distribution"
import "../../diku-dk/linalg/linalg"

module type sketch = {
  -- | The underlying scalar.
  type t

  -- | The underlying distribution the test matrix is sampled from.
  module dist: cbrng_distribution

  -- | A **l**eft **d**ense sketch forming `SA`.
  val sketch_ld [m] [n] : dist.engine.k -> (d: i64) -> [m][n]t -> [d][n]t

  -- | A **r**ight **d**ense sketch forming `AS`.
  val sketch_rd [m] [n] : dist.engine.k -> (d: i64) -> [m][n]t -> [m][d]t

  -- | A **l**eft **o**racle sketch forming `SA`.
  val sketch_lo [m] [n] : dist.engine.k -> (d: i64) -> ([m]t -> [n]t) -> [d][n]t

  -- | A **r**ight **o**racle sketch forming `AS`.
  val sketch_ro [m] [n] : dist.engine.k -> (d: i64) -> ([n]t -> [m]t) -> [m][d]t

  -- | Form the `S` test matrix with *ambient* dimension `n` and *embedding* dimension `d`.
  val sketch : dist.engine.k -> (d: i64) -> (n: i64) -> [d][n]t
}

local module mk_sketch
  (N: numeric)
  (D: cbrng_distribution with num.t = N.t) = {

  type t = N.t
  module dist = D

  -- https://futhark-lang.org/examples/matrix-multiplication.html
  local
  def matmul A B =
    map (\A_row ->
           map (\B_col ->
                  reduce (N.+) (N.i64 0) (map2 (N.*) A_row B_col))
               (transpose B))
        A

    -- | Sketch for the ambient dimensiion `n` and the embedding dimension `d`.
    def sketch dist seed d n : [d][n]t =
      tabulate (d * n) (D.rand seed dist) |> unflatten

    def sketch_ld [m] [n] dist seed d (A: [m][n]t) : [d][n]t =
      let S: [d][m]t = sketch dist seed d m
      in matmul S A

    def sketch_rd [m] [n] dist seed d (A: [m][n]t) : [m][d]t =
      let S: [d][n]t = sketch dist seed d n
      in matmul A (transpose S)

    def sketch_lo [m] [n] dist seed d (oracle: [m]t -> [n]t) : [d][n]t =
      let S: [d][m]t = sketch dist seed d m
      in map (oracle) S

    def sketch_ro [m] [n] dist seed d (oracle: [n]t -> [m]t) : [m][d]t =
      let S: [d][n]t = sketch dist seed d n
      in map (oracle) S |> transpose
}

-- | A sketch which utilizes a Gaussian test matrix sampling from  *X ~ N(0, d^-1)*.
module mk_gaussian_embedding (R: real) (T: integral) (E: cbrng_engine with t = T.t) : sketch with t = R.t = {
  type t = R.t
  module dist = gaussian_distribution R T E

  local module SK = mk_sketch R dist
  local def to_dist d = {mean = R.i64 0, stddev = R.i64 d |> flip (R.**) (R.i64 (-1i64))}

    def sketch_ld seed d = SK.sketch_ld (to_dist d) seed d
    def sketch_rd seed d = SK.sketch_rd (to_dist d) seed d
    def sketch_lo seed d = SK.sketch_lo (to_dist d) seed d
    def sketch_ro seed d = SK.sketch_ro (to_dist d) seed d
    def sketch seed d = SK.sketch (to_dist d) seed d
}
