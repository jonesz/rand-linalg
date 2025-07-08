import "cbrng"

module type cbrng_distribution = {
  -- | The random number engine underlying this distribution.
  module engine: cbrng_engine

  -- | A module describing the type of values produced by this random
  -- distribution.
  module num: numeric

  -- | The dynamic configuration of the distribution.
  type distribution

  -- | Generate a random number given the distributional parameters and the counter.
  val rand : distribution -> i64 -> num.t
}

module rademacher_distribution (D: numeric) (I: integral) (E: cbrng_engine with t = I.t)
  : cbrng_distribution
    with num.t = D.t
    with distribution = () = {
  module engine = E
  module num = D

  type distribution = ()

  def rand _ ctr =
    -- If the underlying RNG is uniform, then the last bit will set 1/2 of the time.
    if E.rand ctr |> I.get_bit 0 |> (==) 0i32
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
