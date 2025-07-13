--| Stochastic Trace Estimation.

import "../../diku-dk/linalg/linalg"
import "../../diku-dk/linalg/qr"

module type ste = {
  type t
  val ste [n] : (m: i64) -> (rand: i64 -> t) -> (matvec: [n]t -> [n]t) -> t
}

-- Generate a random vector for the given `m` iteration.
local def r_vec 't (n: i64) (rand: i64 -> t) (m: i64) : [n]t =
  map (\i -> m * n + i |> rand) (iota n)

-- The naive Girard-Hutchinson trace estimator.
module hutchinson (R: real)
  : ste
    with t = R.t = {
  type t = R.t

  def dotprod [n] (a: [n]t) (b: [n]t) : t =
    -- TODO: This should be a conjugate transpose...
    map2 (R.*) a b |> reduce (R.+) (R.i64 0i64)

  def ste [n] m rand matvec =
    let r_vec ctr = map (\i -> (+) i ctr |> rand) (iota n)
    in map (\m_i ->
              let w = r_vec m_i
              in matvec w |> dotprod w)
           (iota m)
       |> reduce (R.+) (R.i64 0i64)
       |> flip (R./) (R.i64 m)
}

-- https://arxiv.org/pdf/2010.09649
module hutchplusplus (R: real)
  : ste
  with t = R.t = {
    type t = R.t

  module L = mk_linalg R
  module QR = mk_gram_schmidt R

  local def tr A =
    L.fromdiag A |> reduce (R.+) (R.i64 0)
     
  def ste [n] (m: i64) (rand: i64 -> t) (matvec: [n]t -> [n]t) =
    let mot = m / 3 -- "m_over_three"
    let z: [n + n][mot]t = map (r_vec mot rand) (iota (n + n))
    let S: [n][mot]t = take n z :> [n][mot]t
    let G: [n][mot]t = drop n z :> [n][mot]t

    -- A: [n][n]t; S: [n][mot]t; AS = [n][mot]t
    let AS: [n][mot]t = map (\s_i -> matvec s_i) S
    let (Q: [n][n]t, _) = QR.qr AS
    in
      let left = map (matvec) Q |> L.matmul (transpose Q) |> tr
      let right =
        let I_QQT = L.matmul Q (transpose Q) |> L.matsub (L.eye n)
        in L.matmul I_QQT G |> map (matvec) |> L.matmul I_QQT 
        |> L.matmul (transpose G) |> tr |> (R.*) (3_f32 / (f32.i64 m) |> R.f32)
      in (R.+) left right

}

-- module hutchplusplus_na (R: real)
--   : ste
--   with t = R.t = {
--     type t = R.t
-- 
--   module L = mk_linalg R
-- 
--   local def tr A =
--     L.fromdiag A |> reduce (R.+) (R.i64 0)
--      
--   def ste [n] (m: i64) (rand: i64 -> t) (matvec: [n]t -> [n]t) =
--     -- c1 < c2; c1 + c2 + c3 = 1
--     let (c1, c2, c3) = (0.25_f32, 0.50_f32, 0.25_f32)
--     let (m_c1, m_c2, m_c3) =
--       let m_f32 = f32.i64 m
--       in (c1 * m_f32 |> i64.f32, c2 * m_f32 |> i64.f32, c3 * m_f32 |> i64.f32)
-- 
--     let z = map (r_vec n rand) (iota m)
--     let S = take m_c1 z 
--     let R = drop m_c1 z |> take m_c2
--     let G = drop (m_c1 + m_c2) z
-- 
--     let Z = map matvec R
--     let W = map matvec S
-- 
--     -- TODO: Need to compute the pseudoinverse.
--     in ???
-- }
