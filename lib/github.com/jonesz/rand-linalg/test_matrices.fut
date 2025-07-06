-- | ignore

import "../../diku-dk/linalg/linalg"
module L = mk_linalg f32

-- k_ij = exp(-||x_i - x_j||^2 / 2h)
def rbfKernel 't [n] (p: [n]t) h =
  map (\i ->
         map (\j ->
                ???)
             (iota n))
      (iota n)

-- A = diag(1,...,1_R,0,...,0(n-R)) + e/(4n)GG^T; where G is G(0, I_n).
def lowRankNoise rng R n e : [n][n]f32 =
  let r_elem = replicate R 1_f32
  let p_elem = replicate (R - n) 0_f32
  let d = L.todiag (r_elem ++ p_elem :> [n]f32)
  in ???

-- A = diag(1,...,1_R,2^{-p},...,(n-R+1)^{-p})
def polyDecay R n p : [n][n]f32 =
  let r_elem = replicate R 1f32
  let p_elem = map (\i -> f32.i64 i |> flip (**) (f32.neg p)) (2...n - R + 1)
  in r_elem ++ p_elem :> [n]f32
  |> L.todiag

-- A diag(1,...,1_R,10^{-q},...,10^{-(n-R)q})
def expDecay R n q : [n][n]f32 =
  let r_elem = replicate R 1f32
  let p_elem = map (\i -> f32.i64 i |> (*) q |> f32.neg |> (**) 10f32) (1...n - R)
  in r_elem ++ p_elem :> [n]f32
  |> L.todiag

-- Examples from `Tro20-Randomized-Algorithms-LN.pdf`.
def lowRankLowNoise rng R n = lowRankNoise rng R n 0.005_f32
def lowRankMedNoise rng R n = lowRankNoise rng R n 0.05_f32
def lowRankHiNoise rng R n = lowRankNoise rng R n 0.5_f32

def polyDecaySlow R n = polyDecay R n 0.5_f32
def polyDecayMed R n = polyDecay R n 1.0_f32
def polyDecayFast R n = polyDecay R n 2.0_f32

def expDecaySlow R n = expDecay R n 0.01_f32
def expDecayMed R n = expDecay R n 0.1_f32
def expDecayFast R n = expDecay R n 0.5_f32
