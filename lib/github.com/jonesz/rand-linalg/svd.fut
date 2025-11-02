-- | Randomized SVD.
import "../cbrng-fut/distribution"
import "../../diku-dk/linalg/linalg"
import "rangefinder"
import "qr/qr"

-- | Randomized SVD.
module type rsvd = {
  -- | The underlying scalar type.
  type t
  -- | The underlying distribution this sketches from.
  module dist: cbrng_distribution

  -- | Compute the rsvd for a matrix A with `m >> n`.
  -- the argument `l` is the target rank `k` + `p` where `p` could be {k, 5, 10}, etc.
  val rsvd [m] [n] : dist.engine.k -> [m][n]t -> (l: i64) -> ([m][l]t, [l][l]t, [l][n]t)
}

-- | Compute economical SVD via one-sided Jacobi iterations.
-- this module is serial; that is, the scheduling for columns is a simple `foldl`.
module mk_one_sided_jacobi_serial (R: real) : {
  type t
  val svd [l] : [l][l]t -> ([l][l]t, [l][l]t, [l][l]t)
} with t = R.t = {
  type t = R.t

  local module TQ = mk_householder_thin_qr R
  local module LA = mk_linalg R

  -- [Algorithm 6](https://icl.utk.edu/files/publications/2018/icl-utk-1341-2018.pdf)
  def is_orth bij bii bjj =
    let eps = R.f32 0.01f32 -- TODO: What should eps be?
    in (R.<) (R.abs bij) ((R.*) bii bjj |> R.sqrt |> (R.*) eps)

  -- TODO: Work out a parallel schedule!
 def serial_schedule n =
    -- TODO: This filter operation sucks; we know the number of indices per sweep too.
    map (\j -> map (\k -> (j, k)) (iota n)) (iota n) |> flatten |> filter (\(j, k) -> j > k)

  def sweep [l] (S: [l][l]t) (V: [l][l]t) : ([l][l]t, [l][l]t, bool) =
    -- Compute `(c, s, orth)` for the Jacobi rotation of `S` for columns `i, j`.
    let core S i j: (R.t, R.t, bool) =
      let col_i = S[0:, i]
      let col_j = S[0:, j]

      let d_ii = LA.dotprod col_i col_i
      let d_jj = LA.dotprod col_j col_j
      let d_ij = LA.dotprod col_i col_j

      in if (is_orth d_ij d_ii d_jj |> not)

        -- The columns aren't orthogonal; orthogonalize...
        then
          let was_orth = false
          -- [4.2.2 1JAC](www.irisa.fr/sage/bernard/publis/SVD-Chapter06.pdf)
          let alpha = (R.*) (R.i64 2) d_ij
          let beta  = d_jj
          let gamma = (R.+) ((R.*) alpha alpha) ((R.*) beta beta) |> R.sqrt

          in if (R.>) beta (R.i64 0) 

            then
              let c = (R.*) (R.i64 2) gamma |> (R./) ((R.+) beta gamma) |> R.sqrt 
              let s = (R.*) (R.i64 2) gamma |> (R.*) c |> (R./) alpha
              in (c, s, was_orth)

            else 
              let s = (R.*) (R.i64 2) gamma |> (R./) ((R.-) gamma beta) |> R.sqrt
              let c = (R.*) (R.i64 2) gamma |> (R.*) s |> (R./) alpha
              in (c, s, was_orth)

        -- The columns are orthogonal; pass...
        else
          let was_orth = true
          let c = R.i64 0
          let s = R.i64 0
          in (c, s, was_orth)
      
    in foldl (\(S_k, V_k, all_orth) (col_i, col_j) ->
      let (c, s, was_orth) = core S_k col_i col_j

      in if was_orth

        -- If the columns were orthogonal, there's no change to the matrices.
        then (S_k, V_k, and [was_orth, all_orth])

        else
          let S_k_col_i: [l]R.t = S_k[0:, col_i]
          let S_k_col_j: [l]R.t = S_k[0:, col_j]

          let S_k_col_i_new = map2 (\a b -> (R.+) ((R.*) c a) ((R.*) s b)) S_k_col_i S_k_col_j
          let S_k_col_j_new = map2 (\a b -> (R.+) ((R.*) s a |> R.neg) ((R.*) c b)) S_k_col_i S_k_col_j

          -- TODO: Dump the copy somehow?
          let S_k =
              copy S_k with [0:, col_i] = S_k_col_i_new
                       with [0:, col_j] = S_k_col_j_new

          let V_k_col_i: [l]R.t = V_k[0:, col_i]
          let V_k_col_j: [l]R.t = V_k[0:, col_j]

          let V_k_col_i_new = map2 (\a b -> (R.+) ((R.*) c a) ((R.*) s b)) V_k_col_i V_k_col_j
          let V_k_col_j_new = map2 (\a b -> (R.+) ((R.*) s a |> R.neg) ((R.*) c b)) V_k_col_i V_k_col_j

          -- TODO: Dump the copy somehow?
          let V_k =
              copy V_k with [0:, col_i] = V_k_col_i_new
                       with [0:, col_j] = V_k_col_j_new

                in (S_k, V_k, and [all_orth, was_orth])
      ) (S, V, true) (serial_schedule l)

  -- | Return (U, S, V^T); Assume that m == n.
  def svd [n] (A: [n][n]t) : ([n][n]t, [n][n]t, [n][n]t) =
    -- TODO: We should sweep based on some evaluation of the off-diagonal elements
    -- and their convergence to zero... instead we'll sweep ten times.
    let (US, V, _) =
      loop (US, V, all_orth) = (A, LA.eye n, true) while all_orth do
        sweep US V

    let S_diag = transpose US |> map (LA.vecnorm)
    let U = transpose US |> map2 (\sigma_i us_i -> LA.vecscale ((R./) (R.i64 1) sigma_i) us_i) S_diag |> transpose
    let S = LA.todiag S_diag

    -- TODO: Change all negative singular values to positive.
    -- TODO: Sort the singular values.
    in (U, S, (transpose V))
}

module mk_rsvd (R: real) (T: rangefinder with t = R.t) : rsvd with t = R.t = {
  module dist = T.dist
  type t = R.t

  module LA = mk_linalg R
  module TQ = mk_householder_thin_qr R
  module SVD = mk_one_sided_jacobi_serial R

  -- Retun (U, S, V^T); Assuming that m >> n.
  def rsvd [m] [n] seed (A: [m][n]t) (l: i64) =

    -- Compute `Q` with the rangefinder.
    let Q: [m][l]t = T.rangefinder seed l A
    -- Form `C = Q*A`, the approximation; note: l <= n.
    let C: [l][n]t = LA.matmul (transpose Q) A

    -- Following the randomized rangefinder problem `l << n`...
    -- See [Dongarra](https://icl.utk.edu/files/publications/2018/icl-utk-1341-2018.pdf);
    -- "if m << n, it is more efficient to perform an LQ factorization of A".
    -- LQ can be achieved through our `thin_qr` and a transpose.
    let (Q_qr, R_qr) : ([n][l]t, [l][l]t) = transpose C |> TQ.qr ()
    let (Q_qr, L_qr) : ([l][n]t, [l][l]t) = (transpose Q_qr, transpose R_qr)

    let (U_hat, S, V_T) = SVD.svd L_qr

    -- Revert the LQ transformation.
    let V = LA.matmul V_T Q_qr
    -- Form the SVD..
    let U = LA.matmul Q U_hat

    in (U, S, V)
}

import "../cbrng-fut/cbrng"
import "sketch"

-- | RSVD with Gaussian embeddings and the default rangefinder.
module mk_rsvd_default (R: real) : rsvd with t = R.t = {
  local module EMB = mk_gaussian_embedding R u32 squares32
  local module RSVD = mk_rsvd R (mk_rangefinder_dense_default R EMB.dense.right)

  type t = R.t
  module dist = RSVD.dist

  def rsvd = RSVD.rsvd
}

