-- | Sketching operations which transform a matrix from the *ambient*
-- dimension to the *embedding* dimension.

import "../cbrng-fut/cbrng"
import "../cbrng-fut/distribution"
import "../../diku-dk/linalg/linalg"

module type sketch_left_dense = {
  -- | The underlying distribution the test matrix samples from.
  module dist: cbrng_distribution

  -- | The underlying scalar type.
  type t

  -- | A sketch where a dense `A` is passed entirely in memory.
  val sketch [m] [n] : dist.engine.k -> (d: i64) -> [m][n]t -> [d][n]t
}

module type sketch_right_dense = {
  -- | The underlying distribution the test matrix samples from.
  module dist: cbrng_distribution

  -- | The underlying scalar type.
  type t

  -- | A sketch where a dense `A` is passed entirely in memory.
  val sketch [m] [n] : dist.engine.k -> (d: i64) -> [m][n]t -> [m][d]t
}

module type sketch_right_oracle = {
  -- | The underlying distribution the test matrix samples from.
  module dist: cbrng_distribution

  -- | The underlying scalar type.
  type t

  -- | A sketch where an oracle allows access to `A`.
  val sketch [m] [n] : dist.engine.k -> (d: i64) -> ([n]t -> [m]t) -> [m][d]t
}

module type sketch_left_oracle = {
  -- | The underlying distribution the test matrix samples from.
  module dist: cbrng_distribution

  -- | The underlying scalar type.
  type t

  -- | A sketch where an oracle allows access to `A`.
  val sketch [m] [n] : dist.engine.k -> (d: i64) -> ([m]t -> [n]t) -> [d][n]t
}

module type sketch = {
  module dense: {
    module left: {
      include sketch_left_dense
    }

    module right: {
      include sketch_right_dense
    }
  }

  module oracle: {
    module left: {
      include sketch_left_oracle
    }

    module right: {
      include sketch_right_oracle
    }
  }
}

local module mk_sketch (N: numeric) (D: cbrng_distribution with num.t = N.t) = {
  type t = N.t

  -- https://futhark-lang.org/examples/matrix-multiplication.html
  local def matmul A B =
    map (\A_row ->
           map (\B_col ->
                  reduce (N.+) (N.i64 0) (map2 (N.*) A_row B_col))
               (transpose B))
        A

  module dense = {
    module left = {
      def sketch [m] [n] (dist: D.distribution) (seed: D.engine.k) (d: i64) (A: [m][n]t) : [d][n]t =
        tabulate (d * m) (D.rand seed dist) |> unflatten |> flip (matmul) A
    }

    module right = {
      def sketch [m] [n] (dist: D.distribution) (seed: D.engine.k) (d: i64) (A: [m][n]t) : [m][d]t =
        tabulate (n * d) (D.rand seed dist) |> unflatten |> matmul A
    }
  }

  module oracle = {
    module left = {
      def sketch [m] [n] (dist: D.distribution) (seed: D.engine.k) (d: i64) (oracle: [m]t -> [n]t) : [d][n]t =
        tabulate (d * m) (D.rand seed dist) |> unflatten |> map (oracle)
    }

    module right = {
      def sketch [m] [n] (dist: D.distribution) (seed: D.engine.k) (d: i64) (oracle: [n]t -> [m]t) : [m][d]t =
        tabulate (n * d) (D.rand seed dist) |> unflatten |> transpose |> map (oracle) |> transpose
    }
  }
}

-- | A sketch which utilizes a Gaussian test matrix sampling from  *X ~ N(0, d^-1)*.
module mk_gaussian_embedding (R: real) (T: integral) (E: cbrng_engine with t = T.t) : sketch = {
  local module G = gaussian_distribution R T E
  local module S = mk_sketch R G

  local def to_dist d = {mean = R.i64 0, stddev = R.i64 d |> flip (R.**) (R.i64 (-1i64))}

  module dense = {
    module left = {
      module dist = G
      type t = R.t

      def sketch seed d A = S.dense.left.sketch (to_dist d) seed d A
    }

    module right = {
      module dist = G
      type t = R.t

      def sketch seed d A = S.dense.right.sketch (to_dist d) seed d A
    }
  }

  module oracle = {
    module left = {
      module dist = G
      type t = R.t

      def sketch seed d oracle = S.oracle.left.sketch (to_dist d) seed d oracle
    }

    module right = {
      module dist = G
      type t = R.t

      def sketch seed d oracle = S.oracle.right.sketch (to_dist d) seed d oracle
    }
  }
}
