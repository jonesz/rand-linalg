-- | ignore

import "../../diku-dk/linalg/linalg"
import "../cbrng-fut/cbrng"
import "../cbrng-fut/distribution"

module tm (T: real) = {
  type t = T.t

  module L = mk_linalg T
  module N = gaussian_distribution T u32 i64 squares32

  -- k_ij = exp(-||x_i - x_j||^2 / 2h)
  def rbfKernel 't [n] (p: [n]t) h =
    map (\i ->
           map (\j ->
                  ???)
               (iota n))
        (iota n)

  -- A = diag(1,...,1_R,0,...,0(n-R)) + e/(4n)GG^T; where G is G(0, I_n).
  def lowRankNoise seed R n e : [n][n]t =
    let r_elem = replicate R (T.i64 1)
    let p_elem = replicate (n - R) (T.i64 0)
    let d = L.todiag (r_elem ++ p_elem :> [n]t)
    let dist = ((squares32.construct seed), (squares32.construct (seed * 0xA)), {mean = (T.i64 0), stddev = (T.i64 1)})
    let G = map (\i -> map (\j -> N.rand dist (i * n + j)) (iota n)) (iota n)
    let noise = L.matmul G (transpose G)
    let scale = (T./) e (n * 4_i64 |> T.i64)
    in L.matscale scale noise |> L.matadd d

  -- A = diag(1,...,1_R,2^{-p},...,(n-R+1)^{-p})
  def polyDecay R n p : [n][n]t =
    let r_elem = replicate R (T.i64 1)
    let p_elem = map (\i -> T.i64 i |> flip (T.**) (T.neg p)) (2...n - R + 1)
    in r_elem ++ p_elem :> [n]t
    |> L.todiag

  -- A diag(1,...,1_R,10^{-q},...,10^{-(n-R)q})
  def expDecay R n q : [n][n]t =
    let r_elem = replicate R (T.i64 1)
    let p_elem = map (\i -> T.i64 i |> (T.*) q |> T.neg |> (T.**) (T.i64 10)) (1...n - R)
    in r_elem ++ p_elem :> [n]t
    |> L.todiag

  -- Examples from `Tro20-Randomized-Algorithms-LN.pdf`.
  def lowRankLowNoise seed R n = lowRankNoise seed R n (T.f32 0.005_f32)
  def lowRankMedNoise seed R n = lowRankNoise seed R n (T.f32 0.05_f32)
  def lowRankHiNoise seed R n = lowRankNoise seed R n (T.f32 0.5_f32)

  def polyDecaySlow R n = polyDecay R n (T.f32 0.5_f32)
  def polyDecayMed R n = polyDecay R n (T.f32 1.0_f32)
  def polyDecayFast R n = polyDecay R n (T.f32 2.0_f32)

  def expDecaySlow R n = expDecay R n (T.f32 0.01_f32)
  def expDecayMed R n = expDecay R n (T.f32 0.1_f32)
  def expDecayFast R n = expDecay R n (T.f32 0.5_f32)
}
