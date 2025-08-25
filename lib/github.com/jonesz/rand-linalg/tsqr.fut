import "../../diku-dk/linalg/qr"
import "../../diku-dk/linalg/linalg"

-- TODO: Gram Schmidt is too unstable, replace this.
module QR = mk_gram_schmidt f32
module L = mk_linalg f32

-- | TSQR; we required that `k` be a power of 2; `a` should also be
-- divisible by k.
def tsqr [a] [b] (k: i64) (A: [a][b]f32) : ([a][a]f32, [a][b]f32) =
  -- Split into `k` different matrices.
  let A_split: [k][(a / k)][b]f32 = unflatten_3d (flatten A :> [k * (a / k) * b]f32)
  -- The underlying QR algorithm.
  let qr [m] [n] (Z: [m][n]f32): ([m][m]f32, [m][n]f32) = QR.qr Z

  let merge [d] [e] [f] (Rs: [d][e][f]f32): [d / 2][e * 2][f]f32 =
    let Rs: [d * e * f]f32 = flatten_3d Rs
    let Rs = Rs :> [(d / 2) * (e * 2) * f]f32
    in unflatten_3d Rs

  let bound = f32.i64 k |> f32.log2 |> i64.f32
  let R =
    loop (A_split) for _ in 0..<bound do
      let (_, R) = map (qr) A_split |> unzip
      in merge R

  let R = flatten R :> [a][b]f32

  -- TODO: Compute Q at the end.
  in (L.eye a, R)
