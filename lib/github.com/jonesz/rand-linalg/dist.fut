module type distribution = {
  type t
  val rand : i64 -> t
}

module rademacher (R: real) (seed: i64) : distribution = {
  type t = R.t

  def rand i =
    seed + i |> R.i64
}
