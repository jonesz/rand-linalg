import "../cbrng-fut/cbrng"
import "../cbrng-fut/distribution"
import "../../diku-dk/linalg/linalg"

module type sketch_left_dense = {
  module dist: cbrng_distribution
  type t

  -- | A sketch where a dense `A` is passed entirely in memory.
  val sketch [m] [n] : dist.engine.k -> (d: i64) -> [m][n]t -> [d][n]t
}

module type sketch_right_dense = {
  module dist: cbrng_distribution
  type t

  -- | A sketch where a dense `A` is passed entirely in memory.
  val sketch [m] [n] : dist.engine.k -> (d: i64) -> [m][n]t -> [m][d]t
}

module type sketch_left_oracle = {
  module dist: cbrng_distribution
  type t

  -- | A sketch where an oracle allows access to `A`.
  val sketch [m] [n] : dist.engine.k -> (d: i64) -> ([m]t -> [n]t) -> [d][n]t
}

module type sketch_right_oracle = {
  module dist: cbrng_distribution
  type t

  -- | A sketch where an oracle allows access to `A`.
  val sketch [m] [n] : dist.engine.k -> (d: i64) -> ([n]t -> [m]t) -> [m][d]t
}

local module mk_sketch (N: numeric) (D: cbrng_distribution with num.t = N.t) = {
  -- https://futhark-lang.org/examples/matrix-multiplication.html
  local def matmul A B =
    map (\A_row ->
           map (\B_col ->
                  reduce (N.+) (N.i64 0) (map2 (N.*) A_row B_col))
               (transpose B))
        A

  -- | Left sketches: SA
  module left = {
    module dense = {
      module dist = D
      type t = N.t

      def sketch [m] [n] (seed: dist.engine.k) (cfg: dist.distribution) (d: i64) (A: [m][n]t) : [d][n]t =
        tabulate (d * m) (dist.rand seed cfg) |> unflatten |> flip (matmul) A
    }

    module oracle = {
      module dist = D
      type t = N.t

      def sketch [m] [n] (seed: dist.engine.k) (cfg: dist.distribution) (d: i64) (oracle: [m]t -> [n]t) : [d][n]t =
        tabulate (d * m) (dist.rand seed cfg) |> unflatten |> map (oracle)
    }
  }

  -- | Right sketches: AS
  module right = {
    module dense = {
      module dist = D
      type t = N.t

      def sketch [m] [n] (seed: dist.engine.k) (cfg: dist.distribution) (d: i64) (A: [m][n]t) : [m][d]t =
        tabulate (n * d) (dist.rand seed cfg) |> unflatten |> matmul A
    }

    module oracle = {
      module dist = D
      type t = N.t

      def sketch [m] [n] (seed: dist.engine.k) (cfg: dist.distribution) (d: i64) (oracle: [n]t -> [m]t) : [m][d]t =
        tabulate (n * d) (dist.rand seed cfg) |> unflatten |> transpose |> map (oracle) |> transpose
    }
  }
}

-- TODO: https://math.berkeley.edu/~mgu/MA273/RandSIv4a_paper.pdf
module naive_subspace_iteration (D: real) = {
  type t = D.t
  module L = mk_linalg D

  def subspace_iteration [m][n] q (A: [m][n]t) S =
    let AA_T = L.matmul A (transpose A)
    let AA_Tq = replicate q AA_T |> reduce_comm (L.matmul) (L.eye m)
    in L.matmul AA_Tq S
}
