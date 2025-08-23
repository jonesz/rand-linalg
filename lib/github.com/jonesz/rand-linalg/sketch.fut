import "../cbrng-fut/cbrng"
import "../cbrng-fut/distribution"

module mk_sketch (N: numeric) (D: cbrng_distribution with num.t = N.t) = {
  type t = N.t

  -- https://futhark-lang.org/examples/matrix-multiplication.html
  local def matmul A B =
    map (\A_row ->
           map (\B_col ->
                  reduce (N.+) (N.i64 0) (map2 (N.*) A_row B_col))
               (transpose B))
        A

  -- | Left sketches: SA
  module left = {
    -- | A left sketch where a dense `A` is passed entirely in memory.
    module dense = {
      def sketch [m] [n] (seed: D.engine.k) (cfg: D.distribution) (d: i64) (A: [m][n]t) : [d][n]t =
        tabulate (d * m) (D.rand seed cfg) |> unflatten |> flip (matmul) A
    }

    -- | A left sketch where a transpose matrix-vector product oracle is passed.
    module oracle = {
      def sketch [m] [n] (seed: D.engine.k) (cfg: D.distribution) (d: i64) (oracle: [m]t -> [n]t) : [d][n]t =
        tabulate (d * m) (D.rand seed cfg) |> unflatten |> map (oracle)
    }
  }

  -- | Right sketches: AS
  module right = {
    -- | A right sketch where a dense `A` is passed entirely in memory.
    module dense = {
      def sketch [m] [n] (seed: D.engine.k) (cfg: D.distribution) (d: i64) (A: [m][n]t) : [m][d]t =
        tabulate (n * d) (D.rand seed cfg) |> unflatten |> matmul A
    }

    -- | A right sketch where a matrix-vector product oracle is passed.
    module oracle = {
      def sketch [m] [n] (seed: D.engine.k) (cfg: D.distribution) (d: i64) (oracle: [n]t -> [m]t) : [m][d]t =
        tabulate (n * d) (D.rand seed cfg) |> unflatten |> transpose |> map (oracle) |> transpose
    }
  }
}
