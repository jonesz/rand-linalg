-- | Randomized Embeddings.

import "../cbrng-fut/cbrng"
import "../cbrng-fut/distribution"
import "sketch"

module random_sparse_signs_distribution
  (D: numeric)
  (T: integral)
  (E: cbrng_engine with t = T.t)
  : cbrng_distribution
    with num.t = D.t
    with engine.k = E.k
    with distribution = f32 = {
  module engine = E
  module num = D

  type seed = E.k
  type distribution = f32

  def rand seed p ctr =
    -- Determine the cutoffs for the range [E.min, E.max] which correspond to {0, 1, -1}.
    let z = (T.-) E.max E.min |> T.to_i64 |> f32.i64 |> (f32.*) (1f32 - p) |> T.f32
    let o = (T.-) E.max E.min |> T.to_i64 |> f32.i64 |> (f32.*) p |> (f32.*) 0.5_f32 |> T.f32 |> (T.+) z
    let r = E.rand seed ctr
    in if (T.<=) r z
       then D.i64 0
       else if (T.<=) r o
       then D.i64 1
       else D.i64 1 |> D.neg
}

module mk_gaussian_embedding (R: real) (T: integral) (E: cbrng_engine with t = T.t) = {
  type t = R.t
  module S = mk_sketch R (gaussian_distribution R T E)

  local def dist d = {mean = R.i64 0, stddev = (R.i64 d |> flip (R.**) (R.i64 1 |> R.neg))}

  module A = {
    def embed seed d A = S.left.A.sketch seed (dist d) d A
  }
  module B = {
    def embed seed d oracle = S.left.B.sketch seed (dist d) d oracle
  }
}
