--| Stochastic Trace Estimation.

module type ste = {
  type t
  val ste [n] : (m: i64) -> (rand: i64 -> t) -> (matvec: [n]t -> [n]t) -> t
}

-- The naive Girard-Hutchinson trace estimator.
module hutchinson (R: real) : ste = {
  type t = R.t

  def dotprod [n] (a: [n]t) (b: [n]t) : t =
    -- TODO: This should be a conjugate transpose...
    map2 (R.*) a b |> reduce (R.+) (R.i64 0i64)

  def ste [n] m rand matvec =
    let sample_n i = map (\n_i -> (+) n_i i |> rand) (iota n)
    in map (\m_i -> let x = sample_n (n * m_i) in matvec x |> dotprod x) (iota m) |> reduce (R.+) (R.i64 0i64) |> flip (R./) (R.i64 m)
}
