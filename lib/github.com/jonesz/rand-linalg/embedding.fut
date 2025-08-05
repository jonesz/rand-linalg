-- | Randomized Embeddings.

import "../cbrng-fut/cbrng"
import "../cbrng-fut/distribution"

module type randomized_embedding = {
  type t
  type param

  val mul [d] [m] [n] : param -> [m][n]t -> [d][n]t
}

module mk_gaussian_emb (R: real) (E: cbrng_engine) = {
  type t = R.t
  type param = ()
}

module mk_random_signs_emb
  (R: real)
  -- TODO: Remove these constraints on the generated value and the key.
  (E: cbrng_engine with t = u32 with k = i64)
  : randomized_embedding = {
  type t = R.t
  type param = {p: f32, seed: E.k}
  module U = uniform_real_distribution f32 u32 i64 E

  def mul (p: param) A =
    let sample k =
      let dist = (p.seed, {min_r = 0.0_f32, max_r = 1.0_f32})
      let z = U.rand dist k
      in if z <= (1f32 - p.p)
         then R.i64 0
         else if z <= (1f32 - p.p) + (p.p / 2f32)
           then R.i64 1 |> R.neg
           else R.i64 1
    in ???
}

module mk_srft_emb (R: real) (E: cbrng_engine) = {
  type t = R.t

  def mul seed A = ???
}
