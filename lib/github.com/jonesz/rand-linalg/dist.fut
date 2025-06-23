module type distribution = {
  type t
  val rand : i64 -> t
}

module rademacher (R: real) (S: {val seed : i64}) : distribution with t = R.t = {
  type t = R.t

  def rand i =
    -- TODO: utilize an actual random algorithm.
    let cond = S.seed + i |> i64.get_bit 0 |> (==) 1i32
    in if cond
       then (R.i64 (-1i64))
       else (R.i64 1i64)
}
