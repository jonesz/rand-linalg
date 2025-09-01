-- | Operations for performing a QR decomposition. Specifically, these routines target
-- tall and skinny matrices (the result of a sketch) and are the thin
-- variant of QR; that is, their signature is `[m][n]t -> ([m][n]t, [n][n]t)`.
import "../../../diku-dk/linalg/linalg"

-- Taken from the `diku-dk/linalg` / Kasper et al. work.
local def zero_below_main_diag [n] [m] 't (zero: t) (R: [m][n]t) : [m][n]t =
  map2 (\i -> map2 (\j x -> if j < i then zero else x) (iota n))
       (iota m)
       R

-- | The thin QR.
module type thin_qr = {
  -- | The underlying scalar type.
  type t 
  -- | Some configuration variable for this QR (could be none).
  type cfg 

  -- | Perform a QR decomposition, returning `(Q, R)`.
  val qr [n][m] : cfg -> (A: [m][n]t) -> ([m][n]t, [n][n]t)
}

-- | A thin Householder QR.
module mk_householder_thin_qr (D: real) : thin_qr
  with t = D.t
  with cfg = () = {

  local module L = mk_linalg D
  type t = D.t
  type cfg = ()

  -- TODO: There's a construction with `(tau, v)` that should be investigated.`
  local def house [m] (x: [m]D.t) =
    -- If floating point, the sign bit should be flipped on `x[0]`.
    let a =
      let sgn = if (D.>=) x[0] (D.i64 0) then (D.i64 1i64) else (D.i64 (-1i64))
      in L.vecnorm x |> (D.*) sgn
    let u = copy x with [0] = (D.+) x[0] a
    let u = map (\u_i -> (D./) u_i (L.vecnorm u)) u
    in u

  def qr [n] [m] _ (A: [m][n]D.t) : ([m][n]D.t, [n][n]D.t) =
    -- The algorithm's strategy is the normal QR: we compute the `R` matrix
    -- on the forward pass, storing the `k`'th householder in the zero'ed `k-1`'th
    -- column. After forming `R`, we form `Q` with a backward pass.

    -- These update equations/indices are ripped from "3. The QR Decomposition" /
    -- https://acme.byu.edu/00000179-d4cb-d26e-a37b-fffb57800001/qr-decomposition-pdf

    -- Forward pass to form `R`.
    let (u, R) =
      -- Compute the initial householder transformation and store it; the 0'th vector will
      -- not fit into `R`.
      let u: [m][1]D.t = [house A[0:, 0]] |> transpose

      -- Perform the initial copy here. TODO: Place a uniqueness type on `A` `(A: *[m][n]D.t)`?
      -- TODO: I'm relatively sure that `[0:, 0:]` is just the entire matrix; the below throws
      -- some consuming error related to the below loop.
      -- let R = L.matmul (transpose u) A |> L.matmul u |> L.matscale (D.i64 2) |> L.matsub A
      let R = copy A with [0:, 0:] =
        L.matmul (transpose u) A[0:, 0:] |> L.matmul u |> L.matscale (D.i64 2) |> L.matsub A[0:, 0:]

      let R =
        loop R for k in 1..<n do
          let u: [m - k][1]D.t = [house R[k:, k]] |> transpose

          let R =
            -- TODO: Could these operations be chopped up further so that
            -- they're in place too? i.e. we form a new R `with` for each operation.
            let tmp = L.matmul (transpose u) R[k:, k:] |> L.matmul u |> L.matscale (D.i64 2)
            in R with [k:, k:] = L.matsub R[k:, k:] tmp

          -- Store the householder vector in the to-be-zeroed previous column.
          let R = R with [k:, (k - 1)] = flatten u
          in R
      in (u, R)

    -- Backward pass to form `Q`.
    let Q_wip: [m][n]D.t = (L.eye n) ++ (replicate (m - n) (replicate n (D.i64 0))) :> [m][n]D.t
    let Q =
      loop Q = Q_wip for k in (n - 1)..>0 do
        let u: [m - k][1]D.t = [R[k:, (k - 1)]] |> transpose
        let Q =
          let tmp = L.matmul (transpose u) Q[k:, :] |> L.matmul u |> L.matscale (D.i64 2)
          in Q with [k:, :] = L.matsub Q[k:, :] tmp
        in Q

    -- Retrieve the last `u` vector and form the matrix Q.
    let Q =
      -- TODO: I'm relatively sure that `[0:, 0:]` is just the entire matrix.
      -- let tmp = L.matmul (transpose u) Q[0:, :] |> L.matmul u |> L.matscale (D.i64 2)
      -- in Q with [0:, :]
      L.matmul (transpose u) Q |> L.matmul u |> L.matscale (D.i64 2) |> L.matsub Q
    in (Q, (take n R) |> zero_below_main_diag (D.i64 0))
}

-- | The tree-based tall and skinny QR decomposition scheme from [Demmel et al](https://bebop.cs.berkeley.edu/pubs/mhoemmen2008-tsqr-tech-report.pdf).
-- `A` **must be divisble by** `k` **and** `k` **must be a power of two.**
module mk_tsqr (D: real) (QR: { val qr [n] [m] : (A: [m][n]D.t) -> ([m][n]D.t, [n][n]D.t)})
  : thin_qr
  with t = D.t
  with cfg = i64 = {

  type t = D.t
  type cfg = i64
     
  -- | Perform a QR decomposition, splitting `A` into `k` leaves.
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

    let R =
      loop (A_blocked) for _ in 0..<bound do
        let (_, R) = map (QR.qr) A_blocked |> unzip
        in merge R

    -- `R` should be, at this point, `[1][n][m]D.t`.
    let R = flatten R :> [n][m]D.t

    -- TODO: Compute Q from R.
    let Q = ???
    in ???
}
