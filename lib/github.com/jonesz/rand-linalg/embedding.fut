-- | Randomized Embeddings.

import "../cbrng-fut/cbrng"
import "../cbrng-fut/distribution"

module type randomized_embedding = {
  type t
  type param

  val mul [d] [m] [n] : param -> [m][n]t -> [d][n]t
}

module mk_gaussian_emb
  (R: real)
  -- TODO: Remove these constraints on the generated value and the key.
  (E: cbrng_engine with t = u32 with k = i64)
  : randomized_embedding = {
  type t = R.t
  type param = {d: i64, seed1: E.k, seed2: E.k}

  module G = gaussian_distribution R u32 i64 E

  def mul (p: param) A =
    let sample k =
      let dist =
        let mean = R.i64 0
        let stddev = R.i64 p.d |> (R./) (R.i64 1)
        in (p.seed1, p.seed2, {mean=mean, stddev=stddev})
      in G.rand dist k
    in ???
}

module mk_random_signs_emb
  (R: real)
  -- TODO: Remove these constraints on the generated value and the key.
  (E: cbrng_engine with t = u32 with k = i64)
  : randomized_embedding = {
  type t = R.t
  type param = {d: i64, p: f32, seed: E.k}

  module U = uniform_real_distribution f32 u32 i64 E

  def mul (p: param) A =
    let sample k =
      let dist = (p.seed, {min_r = 0.0_f32, max_r = 1.0_f32})
      let z = U.rand dist k

      -- \alpha^2 = 1 / dp, which implies E[S*S] = I_n.
      let alpha = f32.i64 p.d |> (f32.*) p.p |> (f32./) 1f32 |> f32.sqrt |> R.f32
      -- s_ij ~ { +1 w.p. p/2; -1 w.p. p/2; 0 w.p. 1 - p }
      let s = if z <= (1f32 - p.p)
         then R.i64 0
         else if z <= (1f32 - p.p) + (p.p / 2f32)
           then R.i64 1 |> R.neg
           else R.i64 1
      in (R.*) alpha s
    in ???
}

module mk_srft_emb (R: real) (E: cbrng_engine) = {
  type t = R.t

  def mul seed A = ???
}
