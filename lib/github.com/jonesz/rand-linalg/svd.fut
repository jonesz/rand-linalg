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

-- | Compute the SVD via one-sided Jacobi iterations.
module mk_one_sided_jacobi_slow (R: real) : {
  type t
  val svd [l] : [l][l]t -> ([l][l]t, [l][l]t, [l][l]t)
} with t = R.t = {
	type t = R.t
	module LA = mk_linalg R

	-- Compute a 2x2 Schur Decomposition (8.5.2).
	-- https://math.ecnu.edu.cn/~jypan/Teaching/books/2013%20Matrix%20Computations%204th.pdf
	let schurr_decomp A_k (i: i64, j: i64) : (t, t) =
		let col_i = A_k[0:, i]
		let col_j = A_k[0:, j]

		-- form b_ij, b_jj, b_ii
		let (b_ij, b_jj, b_ii) =
			let tmp = map2 (LA.dotprod) [col_i, col_j, col_i] [col_j, col_j, col_i]
			in (tmp[0], tmp[1], tmp[2])

		-- TODO: Stopping criteria should be based on this condition, not number of sweeps.
		-- This is the convergence formula out of the `utk.edu` Dongarra PDF.
		-- let cond =
		-- 	let eps = R.f32 1e-4f32
		-- 	in (R.>=) (R.abs b_ij) ((R.*) b_jj b_ii |> R.sqrt |> (R.*) eps)

		-- TODO: See the note above about stopping criteria.
		-- in if cond
		-- then
			-- See this [algorithm](https://en.wikipedia.org/wiki/Jacobi_rotation#Numerically_stable_computation)
			let tau = (R.-) b_ii b_jj |> flip (R./) ((R.*) (R.i64 2) b_ij)
			let t   = (R.*) tau tau |> (R.+) (R.i64 1) |> R.sqrt |> (R.+) (R.abs tau) |> (R./) (R.sgn tau)
			let cs  = (R.*) t t |> (R.+) (R.i64 1) |> R.sqrt |> (R./) (R.i64 1)
			let sn  = (R.*) cs t
			in (cs, sn)

		-- TODO: See the note above about stopping criteria.
		-- else
		--   (R.i64 1, R.i64 0)

	-- Given the Jacobi rotation parameters, compute the new columns; returns these new columns along
	-- with the entire set of indices to be used in the later `scatter_3d`.
	let new_cols_indices [l] (X: [l][l]t) (x_i: i64) (cs: t) (sn: t) (i: i64) (j: i64): ([l+l]t, [2*l](i64, i64, i64)) =
		-- Compute the new columns.
		let X_col_i_new = map2 (\a b -> (R.+) ((R.*) cs a) ((R.*) sn b)) X[0:, i] X[0:, j]
		let X_col_j_new = map2 (\a b -> (R.+) ((R.*) sn a |> R.neg) ((R.*) cs b)) X[0:, i] X[0:, j]

		let col_pairs =
			map (\col -> tabulate l (\row -> (x_i, row, col))) [i, j] |> flatten 

		in (X_col_i_new ++ X_col_j_new, col_pairs)

	-- A chess round-robin scheduling algorithm; fix the head player, rotate, and then
	-- pair each player at the ends.
	let round_robin [l] (xs: [l]i64) : [l]i64 = [0i64] ++ rotate 1 (tail xs) :> [l]i64
	let pairs [l] (xs: [l]i64) : [l/2i64](i64, i64) =
		tabulate (l / 2i64) (\i ->
			-- `0` should pair with `l`, `1` should pair with `l-1`, so forth...
			let a = drop i xs |> head
			let b = reverse xs |> drop i |> head
			in (a, b)
		)

	-- Compute a single parallel sweep over the entirety of a schedule.
	let parallel_sweep [l] (AV: *[2][l][l]t) =
		let xs = iota l
		let (AV, _) = 
			loop (A_k_V_k, xs) = (AV, xs) for _cycle < l do
				let col_pairs: [l/2](i64, i64) = pairs xs
				-- Compute the rotations Schurr Decomposition on `A_k`.
				let rotations: [l/2](t, t)     = map (schurr_decomp A_k_V_k[0]) col_pairs

				-- We now have the rotation for each column pair -- determine the new
				-- values for the columns in both `A_k` and `V_k`.

				let (A_k_V_k_values, A_k_V_k_indices) =
					let (tmp_0, tmp_1) = map (\i -> map2 (\(cs_i, sn_i) (col_i, col_j) ->
						new_cols_indices A_k_V_k[i] i cs_i sn_i col_i col_j) rotations col_pairs |> unzip) (iota 2) |> unzip
					in (flatten_3d tmp_0 :> [l * 2 * l]t, flatten_3d tmp_1 :> [l * 2 * l](i64, i64, i64))

			 	let A_k_V_k = scatter_3d A_k_V_k A_k_V_k_indices A_k_V_k_values
				in (A_k_V_k, round_robin xs)

		in AV

    -- | Compute (U, S, V_T).
	def svd [l] (A: [l][l]t) : ([l][l]t, [l][l]t, [l][l]t) =
		let A_hat_V =
			-- TODO: See the note about convergence above. 10~ sweeps should be enough.
			loop A_k_V_k = [copy A, LA.eye l] for _iter < 10 do
				parallel_sweep A_k_V_k

		let A_hat = A_hat_V[0]
		let V = A_hat_V[1]

        -- The singular values are the euclidean norms of each column.
		let S = transpose A_hat |> map (LA.vecnorm) |> LA.todiag

        -- Compute U; (A_hat = US) == (A_hat * S_inv = U)
		let U =
			let S_inv = LA.fromdiag S |> map (\x -> (R./) (R.i64 1) x) |> LA.todiag
			in LA.matmul A_hat S_inv

		in (U, S, (transpose V))
}

module mk_rsvd (R: real) (T: rangefinder with t = R.t) : rsvd with t = R.t = {
  module dist = T.dist
  type t = R.t

  module LA = mk_linalg R
  module TQ = mk_householder_thin_qr R
  module SVD = mk_one_sided_jacobi_slow R

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
  local module RSVD = mk_rsvd R (mk_rangefinder_dense_default R EMB)

  type t = R.t
  module dist = RSVD.dist

  def rsvd = RSVD.rsvd
}

