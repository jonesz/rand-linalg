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

-- | A tall and skinny thin QR factorization.
module TSQR (D: real)
  : {
      val qr [n] [m] : i64 -> (A: [n][m]D.t) -> ([n][m]D.t, [n][n]D.t)
    } = {

  -- The "6.1 Structured QR factorizations" thin QR algorithm.
  module local_qr : {
    val qr [n][m] : (A: [m][n]D.t) -> ([m][n]D.t, [n][n]D.t) 
  } = {
    def qr [n][m] (A: [m][n]D.t) : ([m][n]D.t, [n][n]D.t) =
      -- Compute a householder vector.

      let house [m] (x: [m]D.t) =
        let e_1 = [(D.i64 1)] ++ replicate (m-1) (D.i64 0) 
        let norm z = map2 (D.*) z z |> reduce (D.+) (D.i64 0) |> D.sqrt

        -- If using floating point arithmetic, we flip the sign bit dependent
        -- on the first element in the vector.
        let alpha = if (D.sgn x[0]) >= 0 then norm x |> D.neg else norm x 

        -- TODO: Can we not use the copy here?
        let u = copy x with [0] = x[0] - alpha
        let v = map (D./) u (norm u)
        in v

      let (Q, R) = ???
  }
  
  -- | TSQR; `k` must be a power of 2 and the number of rows in `A` (`[n]`) must be divisible by `k`.
  def qr [n] [m] k (A: [n][m]D.t) : ([n][n]D.t, [n][m]D.t) =

    -- Split the entire matrix `A` into `k` blocks.
    let A_blocked: [k][(n / k)][m]D.t =
      let A: [n * m]D.t = flatten A
      let A: [k * (n / k) * m]D.t = A :> [k * (n / k) * m]D.t
      in unflatten_3d A

    -- At each new depth of the tree, we concatenate each pair of the block matrices together.
    let merge [d] [e] [f] (X: [d][e][f]D.t): [d / 2][e * 2][f]D.t =
      let X: [d * e * f]D.t = flatten_3d X
      let X: [(d / 2) * (e * 2) * f]D.t = X :> [(d / 2) * (e * 2) * f]D.t
      in unflatten_3d X

    -- Depth of the tree.
    let bound = f32.i64 k |> f32.log2 |> i64.f32

    -- Underlying QR algorithm to call.
    let qr_econ [i] [j] (X: [i][j]D.t): ([i][i]D.t, [i][j]D.t) = ???

    let R =
      loop (A_blocked) for _ in 0..<bound do
        let (_, R) = map (qr_econ) A_blocked |> unzip
        in merge R

    -- `R` should be, at this point, `[1][n][m]D.t`.
    let R = flatten R :> [n][m]D.t

    -- TODO: Compute Q from R.
    let Q = ???
    in (Q, R)
}
