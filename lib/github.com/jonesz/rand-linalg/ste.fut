--| Stochastic Trace Estimation.

module type ste = {
  type t
  type dist
  val ste [n] : (m: i32) -> (rand: dist -> i32 -> t) -> (matvec: [n]t -> [n]t) -> t
}

-- The naive Girard-Hutchinson trace estimator.
module hutchinson (R: real)
  : ste
    with t = R.t
    with dist = () = {
  type t = R.t
  type dist = ()

  def dotprod [n] (a: [n]t) (b: [n]t) : t =
    -- TODO: This should be a conjugate transpose...
    map2 (R.*) a b |> reduce (R.+) (R.i64 0i64)

  def ste [n] (m: i32) rand matvec =
    let sample_n i = map (\n_i -> (+) n_i i |> i32.i64 |> rand ()) (iota n)
    in map (\m_i -> let x = sample_n (n * m_i) in matvec x |> dotprod x) (i64.i32 m |> iota) |> reduce (R.+) (R.i64 0i64) |> flip (R./) (R.i32 m)
}
