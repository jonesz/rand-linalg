-- | Randomized SVD.

import "rangefinder"
import "../cbrng-fut/distribution"
import "../../diku-dk/linalg/linalg"
import "qr/qr"

module one_sided_jacobi_serial (R: real) = {
  type t = R.t

  module TQ = mk_householder_thin_qr R
  module LA = mk_linalg R

  -- [Algorithm 6](https://icl.utk.edu/files/publications/2018/icl-utk-1341-2018.pdf)
  def is_orth bij bii bjj =
    let eps = R.f32 0.01f32 -- TODO: What should eps be?
    in (R.<) (R.abs bij) ((R.*) bii bjj |> R.sqrt |> (R.*) eps)

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
        -- ORTHOGONALIZE
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

        -- WAS ALREADY ORTHOGONAL.
        else
          let was_orth = true
          let c = R.i64 0
          let s = R.i64 0
          in (c, s, was_orth)
      
    in foldl (\(S_k, V_k, all_orth) (col_i, col_j) ->
      let (c, s, was_orth) = core S_k col_i col_j

      in if was_orth
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

  -- Return (U, S, V^T)
  def svd_econ [l] (A: [l][l]t) : ([l][l]t, [l][l]t, [l][l]t) =
    -- Following the randomized rangefinder problem `l << n`; we need a square
    -- matrix before the two-sided jacobi iteration: compute the LQ factorization.
    -- See [Dongarra](https://icl.utk.edu/files/publications/2018/icl-utk-1341-2018.pdf);
    -- "if m << n, it is more efficient to perform an LQ factorization of A".
    -- LQ can be achieved through our `thin_qr` and a transpose.

    -- TODO: With a one-sided Jacobi, we work on A^T * A without explicitly forming it;
    -- computing the initial QR probably isn't necessary?

    -- let (Q, R) = transpose A |> TQ.qr ()
    -- `thin_qr([n][l]t)``
    -- let (Q, L) = (transpose Q, transpose R)
    -- At this point, we compute the SVD on the matrix `L`, which is `[l][l]t` (square)
    -- and suitable for a two-sided jacobi.

    -- TODO: We should sweep based on some evaluation of the off-diagonal elements
    -- and their convergence to zero... instead we'll sweep ten times.
    let (S, V, _) =
      loop (S, V, all_orth) = (A, LA.eye l, true) while all_orth do
        sweep S V

    let U = LA.eye l
    -- TODO: Sort the singular values?
    -- TODO: Change all negative singular values to positive.
    -- TODO: Pad the matrix with zero columns.
    -- TODO: Determine whether this construction of V is correct.
    -- let V = LA.matmul (transpose V) Q
    in (U, S, V)
}

module type svd = {
  module dist: cbrng_distribution
  type t
  val svd_econ [m] [n] : dist.engine.k -> [m][n]t -> (l: i64) -> ([m][l]t, [l][l]t, [l][n]t)
}

-- module randomized_svd (R: real) (T: rangefinder with t = R.t) : svd with t = R.t = {
--   module dist = T.dist
--   type t = R.t
-- 
--   module L = mk_linalg R
--   module SVD = one_sided_jacobi_serial R
-- 
--   -- Assuming that m >> n.
--   def svd_econ [m] [n] seed (A: [m][n]t) (l: i64) =
--     -- Form `Q = qr_econ(AQ)`.
--     let Q: [m][l]t = T.rangefinder seed l A
--     -- Form `C = Q*A`, the approximation; note: l <= n.
--     let C: [l][n]t = L.matmul (transpose Q) A
--     let (U_hat, S, V): ([l][l]t, [l][l]t, [l][n]t) = SVD.svd_econ C
--     let U = L.matmul Q U_hat
--     in (U, S, V)
-- }
