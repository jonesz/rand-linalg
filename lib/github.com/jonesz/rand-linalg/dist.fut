module type distribution = {
  type t
  val rand : i64 -> t
}

module rademacher (R: real) (S: { val seed: i64 }) : distribution = {
  type t = R.t

  def rand i =
    S.seed + i |> R.i64
}
