--| Stochastic Trace Estimation.

module type ste = {
  type t
  type dist

  val ste [n] : (m: i64) -> (rand: dist -> i64 -> t) -> (matvec: [n]t -> [n]t) -> t
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

  def ste [n] m rand matvec =
    let r_vec ctr = map (\i -> (+) i ctr |> rand ()) (iota n)
    in map (\m_i ->
              let w =
                (*) n m_i |> r_vec
              in matvec w |> dotprod w)
           (iota m)
       |> reduce (R.+) (R.i64 0i64)
       |> flip (R./) (R.i64 m)
}
