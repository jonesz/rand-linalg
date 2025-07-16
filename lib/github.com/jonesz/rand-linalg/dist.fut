import "cbrng"

module type cbrng_distribution = {
  -- | The random number engine underlying this distribution.
  module engine: cbrng_engine

  -- | A module describing the type of values produced by this random
  -- distribution.
  module num: numeric

  -- | The dynamic configuration of the distribution.
  type distribution

  -- | A constructed configured distribution.
  type configuration

  -- | Construct a distribution given a key and a distribution.
  val construct : engine.k -> distribution -> configuration

  -- | Generate a random number given the the counter.
  val rand : configuration -> i64 -> num.t
}

module rademacher_distribution (D: numeric) (I: integral) (K: integral) (E: cbrng_engine with t = I.t with k = K.t)
  : cbrng_distribution
    with num.t = D.t
    with distribution = ()
    with configuration = K.t
    with engine.k = K.t = {
  module engine = E
  module num = D

  type distribution = ()
  type configuration = K.t

  def construct k _ = k

  def rand k ctr =
    -- If the underlying RNG is uniform, then the last bit will set 1/2 of the time.
    if E.rand k ctr |> I.get_bit 0 |> (==) 0i32
    then D.i32 1i32
    else D.i32 1i32 |> D.neg
}

-- | Normally distributed floats.
module normal_distribution
  (R: real)
  (I: integral)
  (K: integral)
  (E: cbrng_engine with t = I.t with k = K.t)
  : cbrng_distribution
    with distribution = {mean: R.t, stddev: R.t}
    with engine.k = K.t
    with num.t = R.t
    with configuration = (K.t, {mean: R.t, stddev: R.t}) = {
  def to_R (x: E.t) =
    R.u64 (u64.i64 (I.to_i64 x))

  module engine = E
  module num = R
  type distribution = {mean: num.t, stddev: num.t}
  type configuration = (K.t, distribution)

  def construct k dist = (k, dist)

  def normal (mean: num.t) (stddev: num.t) = {mean = mean, stddev = stddev}

  open R

  def rand (key: K.t, {mean, stddev}: distribution) ctr =
    let k1 = key
    let k2 = (K.+) key (K.i64 1337i64)
    -- Straight port from `diku-dk/cpprandom/random.fut`.
    -- Box-Muller where we only use one of the generated points.
    let u1 = E.rand k1 ctr
    let u2 = E.rand k2 ctr
    let u1 = (to_R u1 - to_R E.min) / (to_R E.max - to_R E.min)
    let u2 = (to_R u2 - to_R E.min) / (to_R E.max - to_R E.min)
    let r = sqrt (i32 (-2) * log u1)
    let theta = i32 2 * pi * u2
    in mean + stddev * (r * cos theta)
}

-- | This uniform integer distribution generates floats in a given
-- range with "equal" probability for each.
-- A straight port from `diku-dk/cpprandom`.
module uniform_real_distribution (R: real) (I: integral) (K: integral) (E: cbrng_engine with t = I.t with k = K.t)
  : cbrng_distribution
    with num.t = R.t
    with engine.k = K.t
    with distribution = (R.t, R.t)
    with configuration = (K.t, (R.t, R.t)) = {
  def to_R (x: E.t) =
    R.u64 (u64.i64 (I.to_i64 x))

  module engine = E
  module num = R

  type distribution = (R.t, R.t)
  type configuration = (K.t, (R.t, R.t))

  def construct k b = (k, b)

  def rand (k: K.t, (min_r: R.t, max_r: R.t)) ctr =
    let x = E.rand k ctr
    let x' = R.((to_R x - to_R E.min) / (to_R E.max - to_R E.min))
    in R.(min_r + x' * (max_r - min_r))
}

module uniform_int_distribution (D: numeric) (I: integral) (K: integral) (E: cbrng_engine with t = I.t with k = K.t)
  : cbrng_distribution
    with num.t = D.t
    with engine.k = K.t
    with distribution = (D.t, D.t)
    with configuration = (K.t, (D.t, D.t)) = {
  module engine = E
  module num = D
  module U = uniform_real_distribution f32 I K E

  type distribution = (D.t, D.t)
  type configuration = (K.t, (D.t, D.t))

  def construct k (a, b) = (k, (a, b))

  def rand (k: K.t, (D_min: D.t, D_max: D.t)) ctr =
    -- Generate a number between `[0, 1)`, scale it by the range, floor to `{0, ..., N-1}`,
    -- then shift by `D_min`.
    let x = U.rand (k, (0f32, 1f32)) ctr
    in (D.-) D_max D_min |> (D.+) (D.i64 1) |> D.to_i64 |> f32.i64 |> (f32.*) x |> f32.floor |> D.f32 |> (D.+) D_min
}
