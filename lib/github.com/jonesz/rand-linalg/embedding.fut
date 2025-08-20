-- | Randomized Embeddings.

import "../cbrng-fut/cbrng"
import "../cbrng-fut/distribution"
import "../../diku-dk/linalg/linalg"

module random_sparse_signs_distribution
  (D: numeric)
  (T: integral)
  (E: cbrng_engine with t = T.t)
  : cbrng_distribution
    with engine.k = E.k
    with num.t = D.t
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

module type embedding = {
  module engine: cbrng_engine
  type t

  val embed [m] [n] : engine.k -> (d: i64) -> [m][n]t -> [d][n]t
}

module guassian_embedding (R: real) (T: integral) (E: cbrng_engine with t = T.t)
  : embedding
    with engine.k = E.k
    with t = R.t = {
  module engine = E
  type t = R.t

  module L = mk_linalg R
  module G = gaussian_distribution R T E

  def embed [m] seed d A =
    let dist = {mean = R.i64 0, stddev = (R.i64 d |> flip (R.**) (R.i64 1 |> R.neg))}
    let S = tabulate (d * m) (G.rand seed dist) |> unflatten
    in L.matmul S A
}
