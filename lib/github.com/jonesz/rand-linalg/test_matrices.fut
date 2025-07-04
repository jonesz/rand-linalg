-- | ignore

import "../../diku-dk/linalg/linalg"
module L = mk_linalg f32

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
def polyDecaySlow R n = polyDecay R n 0.5_f32
def polyDecayMed R n = polyDecay R n 1.0_f32
def polyDecayFast R n = polyDecay R n 2.0_f32

def expDecaySlow R n = expDecay R n 0.01_f32
def expDecayMed R n = expDecay R n 0.1_f32
def expDecayFast R n = expDecay R n 0.5_f32
