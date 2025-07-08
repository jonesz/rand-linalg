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
    with configuration = K.t
    with distribution = ()
    with num.t = D.t
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

-- -- | Normally distributed floats.
-- module normal_distribution (R: real) (I: integral) (E: cbrng_engine with t = I.t): cbrng_distribution with = {
--   module engine = E
--   module num = R
--   type distribution = {mean: num.t, stddev: num.t}
--
--   def normal (mean: num.t) (stddev: num.t) = {mean = mean, stddev = stddev}
-- }
