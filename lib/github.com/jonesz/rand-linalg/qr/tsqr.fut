import "../../../diku-dk/linalg/qr"
import "../../../diku-dk/linalg/linalg"

local def zero_below_main_diag [n] [m] 't (zero: t) (R: [m][n]t) : [m][n]t =
  map2 (\i -> map2 (\j x -> if j < i then zero else x) (iota n))
       (iota m)
       R

-- | An economical QR.
module econ_qr (D: real) = {
  module L = mk_linalg D

  def house [m] (x: [m]D.t) =
    -- If floating point, the sign bit should be flipped on `x[0]`.
    let a =
      let sgn = if (D.>=) x[0] (D.i64 0) then (D.i64 1i64) else (D.i64 (-1i64))
      in L.vecnorm x |> (D.*) sgn
    -- TODO: Remove the copy?
    let u = copy x with [0] = (D.+) x[0] a
    let u = map (\u_i -> (D./) u_i (L.vecnorm u)) u
    in u

  def qr [n] [m] (A: [m][n]D.t) : ([m][n]D.t, [n][n]D.t) =
    -- This algorithm's strategy is the normal QR: we compute the `R` matrix
    -- on the foreward pass, storing the `k`'th householder in the zero'ed `k-1`'th
    -- column. After forming `R`, we form `Q` with a backward pass.
    let (u, R) =
      -- Compute the initial householder transformation; the 0'th vector will not
      -- fit into `R`.
      let u: [m][1]D.t = [house A[0:, 0]] |> transpose
      -- TODO: Remove the copy?
      let R = copy A with [0:, 0:] = L.matmul (transpose u) A[0:, 0:] |> L.matmul u |> L.matscale (D.i64 2) |> L.matsub A[0:, 0:]
      let R =
        loop R for k in 1..<n do
          let u: [m - k][1]D.t = [house A[k:, k]] |> transpose
          -- TODO: Remove the copy?
          let R = copy R with [k:, k:] = L.matmul (transpose u) R[k:, k:] |> L.matmul u |> L.matscale (D.i64 2) |> L.matsub R[k:, k:]
          let R = R with [k:, (k - 1)] = flatten u
          in R
      in (u, R)
    -- Backward pass to form Q.
    let Q_wip: [m][n]D.t = (L.eye n) ++ (replicate (m - n) (replicate n (D.i64 0))) :> [m][n]D.t
    let (Q, R) =
      loop (Q, R) = (Q_wip, R)
      for k in (n - 1)..>0 do
        let u: [m - k][1]D.t = [R[k:, (k - 1)]] |> transpose
        let Q = copy Q with [k:, :] = L.matmul (transpose u) Q[k:, :] |> L.matmul u |> L.matscale (D.i64 2) |> L.matsub Q[k:, :]
        in (Q, R)
    -- Last `u` vector...
    let Q = copy Q with [0:, :] = L.matmul (transpose u) Q[0:, :] |> L.matmul u |> L.matscale (D.i64 2) |> L.matsub Q[0:, :]
    in (Q, (take n R) |> zero_below_main_diag (D.i64 0))
}

-- | A tall and skinny thin QR factorization.
module TSQR (D: real)
  : {
      val qr [n] [m] : i64 -> (A: [m][n]D.t) -> ([m][n]D.t, [n][n]D.t)
    } = {
  module QR_E = econ_qr D

  -- | TSQR; `k` must be a power of 2 and the number of rows in `A` (`[n]`) must be divisible by `k`.
  def qr [n] [m] k (A: [m][n]D.t) : ([m][n]D.t, [n][n]D.t) =
    -- Split the entire matrix `A` into `k` blocks.
    let A_blocked: [k][(n / k)][m]D.t =
      let A: [m * n]D.t = flatten A
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
    let qr_econ = QR_E.qr
    let R =
      loop (A_blocked) for _ in 0..<bound do
        let (_, R) = map (qr_econ) A_blocked |> unzip
        in merge R
    -- `R` should be, at this point, `[1][n][m]D.t`.
    let R = flatten R :> [n][m]D.t
    -- TODO: Compute Q from R.
    let Q = ???
    in ???
}
