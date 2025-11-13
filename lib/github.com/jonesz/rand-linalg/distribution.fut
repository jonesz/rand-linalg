-- | Distributions useful for random sketches.
import "../cbrng-fut/cbrng"
import "../cbrng-fut/distribution"
import "sketch"

-- | A distribution where `s_ij ~ { +1 ~ p/2; -1 ~ p/2; 0 ~ 1 - p}`.`
module mk_random_sparse_signs
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

  -- | A number between [0, 1) (bernoulli `p`).
  type distribution = f32

  -- | Generate a random number given the seed, a bernoulli `p`, and a counter.
  def rand seed p ctr =
    -- Determine the cutoffs for the range `[E.min, E.max]` which correspond to `{0, 1, -1}`.
    let e_min_f32 = T.to_i64 E.min |> f32.i64
    let range = (T.-) E.max E.min |> T.to_i64 |> f32.i64

    -- We've determned the length of the range -- we need to move the left bound
    -- to `E.min`.
    let a_offset = (f32.*) range (1_f32 - p) |> (f32.+) e_min_f32
    let b_offset = (f32.*) range p |> (f32.*) 0.5_f32 |> (f32.+) a_offset

    let a_offset = T.f32 a_offset
    let b_offset = T.f32 b_offset

    let rnd = E.rand seed ctr
    in if (T.<) rnd a_offset
      then D.i64 0
      else if (T.<) rnd b_offset
        then D.i64 1
        else D.i64 1 |> D.neg
}
