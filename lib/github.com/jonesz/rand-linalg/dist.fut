import "../../diku-dk/cpprandom/random"

module type cbrng_distribution = {
  module engine: rng_engine
  module num: numeric

  type distribution

  val rand : distribution -> i32 -> num.t
}

module rademacher_distribution (D: real) (E: rng_engine) (S: {val seed : i32})
  : cbrng_distribution
    with num.t = D.t
    with engine.rng = E.rng
    with distribution = () = {
  module engine = E
  module num = D

  type distribution = ()

  def rand _ i =
    let rng = E.rng_from_seed [S.seed + i]
    let (_, x) = E.rand rng
    let cond = (E.int.get_bit 0 x) |> (==) 0i32
    in if cond
       then D.i32 1i32 |> D.neg
       else D.i32 1i32
}
