-- | Randomized SVD.
import "../cbrng-fut/distribution"
import "../../diku-dk/linalg/linalg"
import "rangefinder"
import "qr/qr"
import "sketch"

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

-- | A randomized turnstile SVD.
module type rsvd_turnstile = {
	-- | The underlying scalar type.
	type t
	-- | The underlying distribution this sketches from.
	module dist: cbrng_distribution

	type seeds = (dist.engine.k, dist.engine.k, dist.engine.k, dist.engine.k)

	-- | Initiliaze `X, Y, Z` with sketch size parameter `k` where `r <= k << min(m, n)`.
	val initialize : (m: i64) -> (n: i64) -> (k: i64) -> (seed: i64)
		-> ([k][n]t, [m][k]t, [2 * k][2 * k]t, seeds)

	val linear_update_col [m] [n] [k] : seeds -> (X: *[k][n]t) -> (Y: *[m][k]t) -> (Z: *[2 * k][2 * k]t) -> (col: i64) -> (H: [m]t)
		-> ([k][n]t, [m][k]t, [2 * k][2 * k]t)

	val linear_update_row [m] [n] [k] : seeds -> (X: *[k][n]t) -> (Y: *[m][k]t) -> (Z: *[2 * k][2 * k]t) -> (row: i64) -> (H: [n]t)
		-> ([k][n]t, [m][k]t, [2 * k][2 * k]t)

	val linear_update_entry [m] [n] [k] : seeds -> (X: *[k][n]t) -> (Y: *[m][k]t) -> (Z: *[2 * k][2 * k]t) -> (idx: (i64, i64)) -> (H: t)
		-> ([k][n]t, [m][k]t, [2 * k][2 * k]t)

	-- TODO: `r` should be the target rank, which should be a little bit below `k`.
	val sketchy_svd [m] [n] [k] : seeds -> (X: [k][n]t) -> (Y: [m][k]t) -> (Z: [2 * k][2 * k]t) -> (_r: i64)
		-> ([m][k]t, [k][k]t, [k][n]t)
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
	-- with the entire set of indices to be used in the later `scatter_2d`.
	let new_cols_indices [l] (X: [l][l]t) (cs: t) (sn: t) (i: i64) (j: i64): ([l+l]t, [l+l](i64, i64)) =
		-- Compute the new columns.
		let X_col_i_new = map2 (\a b -> (R.+) ((R.*) cs a) ((R.*) sn b)) X[0:, i] X[0:, j]
		let X_col_j_new = map2 (\a b -> (R.+) ((R.*) sn a |> R.neg) ((R.*) cs b)) X[0:, i] X[0:, j]

		let col_pairs =
			map (\col -> tabulate l (\row -> (row, col))) [i, j] |> flatten :> [l+l](i64, i64)

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
	let parallel_sweep [l] (A_k: *[l][l]t) (V_k: *[l][l]t) =
		let xs = iota l
		let (A_k, V_k, _) = 
			loop (A_k, V_k, xs) = (A_k, V_k, xs) for _cycle < l do
				let col_pairs: [l/2](i64, i64) = pairs xs
				let rotations: [l/2](t, t)     = map (schurr_decomp A_k) col_pairs

				-- We now have the rotation for each column pair -- determine the new
				-- values for the columns in both `A_k` and `V_k`.

				let (A_k_col_values, A_k_col_indices, V_k_col_values, V_k_col_indices) =
					let tmp =
						map (\X -> map2 (\(cs_i, sn_i) (col_i, col_j) ->
							new_cols_indices X cs_i sn_i col_i col_j) rotations col_pairs |> unzip) [A_k, V_k]
					in (flatten tmp[0].0, flatten tmp[0].1, flatten tmp[1].0, flatten tmp[1].1)


				-- NOTE: It might be tempting to attempt a `scatter_3d` here, combining
				-- `A_k` and `V_k` into a single `[2][l][l]t` vector; well I tried it
				-- and it turned out to be slower (see bookmark `osj_scatter_3d`).

				let A_k = scatter_2d A_k A_k_col_indices A_k_col_values
				let V_k = scatter_2d V_k V_k_col_indices V_k_col_values

				in (A_k, V_k, round_robin xs)

		in (A_k, V_k)

    -- | Compute (U, S, V_T).
	def svd [l] (A: [l][l]t) : ([l][l]t, [l][l]t, [l][l]t) =
		let (A_hat, V) =
			-- TODO: See the note about convergence above. 10~ sweeps should be enough.
			loop (A_k, V_k) = (copy A, LA.eye l) for _iter < 10 do
				parallel_sweep A_k V_k

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

module mk_sketchy_svd(R: real) (SK: sketch with t = R.t): rsvd_turnstile with t = R.t = {
	type t = R.t

	module dist = SK.dist
	type seeds = (dist.engine.k, dist.engine.k, dist.engine.k, dist.engine.k)

	module LA = mk_linalg R
  	module TQ = mk_householder_thin_qr R
  	module SVD = mk_one_sided_jacobi_slow R

	def initialize m n k seed =
		let X = LA.matzeros k n
		let Y = LA.matzeros m k
		let Z = LA.matzeros (2i64 * k) (2i64 * k)
		let s = (dist.engine.construct seed, dist.engine.construct (seed + 1),
				 dist.engine.construct (seed + 2), dist.engine.construct (seed + 3))
		in (X, Y, Z, s)

	-- [Slide 18](https://tropp.caltech.edu/slides/Tro21-SketchySVD-slides.pdf)
	def linear_update_col [m] [n] [k] (seed: seeds) (X: *[k][n]t) (Y: *[m][k]t) (Z: *[2 * k][2 * k]t) (col: i64) (H: [m]t) =
		-- X <- X + SEED1 H
		-- Y <- Y + H SEED2
		-- Z <- Z + SEED3 H SEED4

		-- X <- X + SEED1 (k x m) * H (m x 1)
		let X =
			let omega = SK.sketch seed.0 k m
			let X_new_col = LA.matvecmul_row omega H |> map2 (R.+) (copy X[0:, col])
			in X with [0:, col] = X_new_col

		-- Y <- Y + H (m x 1) * SEED2 (n x k)
		-- NOTE: H is actually [m x n], but we have a single slice of it.
		let Y =
			let omega = SK.sketch seed.1 k n |> transpose
			let Y_new_mat =
				tabulate_2d m k (\m_i k_i -> (R.*) H[m_i] omega[col, k_i])
			in LA.matadd Y Y_new_mat

		--      [2 * k][2 * k]   [2*k][m]  [m x 1]t    [n][2*k]
		-- Z <- Z              + SEED3    * H        * SEED4
		-- NOTE: H is actually [2 * k][n] after the multiplication. TODO: Does this present an issue?
		let Z =
			let omega = SK.sketch seed.2 (2 * k) m 
			let rho   = SK.sketch seed.3 (2 * k) n |> transpose
			let Z_new_mat =
				let left = LA.matvecmul_row omega H
				in tabulate_2d (2*k) (2*k) (\m_i k_i -> (R.*) left[m_i] rho[col, k_i])
			in LA.matadd Z Z_new_mat

		in (X, Y, Z)

	def linear_update_row seeds X Y Z idx H =
		???

	def linear_update_entry seeds X Y Z (idx_row, idx_column) H = 
		???

	def sketchy_svd [k] [m] [n] (seed: seeds) (X: [k][n]t) (Y: [m][k]t) (Z: [2 * k][2 * k]t) r =
		let Q: [m][k]t =
			let (Q, _) = TQ.qr () Y
			in Q

		let P: [n][k]t =
			let (P, _) = transpose X |> TQ.qr ()
			in P

		let C =
			-- TODO: If we perform the LQ decomposition on the right side,
			-- we can perform the SVDs on both at the same time.
			let left =
				let tmp: [2 * k][k]t = SK.sketch_ld seed.2 (2 * k) Q
				let (tmp_Q, tmp_R) = TQ.qr () tmp
				let (U, S, V_T) = SVD.svd tmp_R
				let S_inv = LA.fromdiag S |>
					map (\x ->
						if (R.==) x (R.i64 0)
							then (R.i64 0)
							else (R./) (R.i64 1) x
					)
				|> LA.todiag
				let R_inv = LA.matmul (transpose V_T) S_inv |> flip (LA.matmul) (transpose U)
				in LA.matmul R_inv (transpose tmp_Q)

			let right =
				let tmp: [k][2 * k]t = SK.sketch_rd seed.3 (2 * k) (transpose P)
				let (tmp_Q, tmp_R) = transpose tmp |> TQ.qr ()
				let tmp_L = transpose tmp_R

				let (U, S, V_T) = SVD.svd tmp_L
				let S_inv = LA.fromdiag S |>
					map (\x ->
						if (R.==) x (R.i64 0)
						then (R.i64 0)
						else (R./) (R.i64 1) x
					)
				|> LA.todiag
				let R_inv = LA.matmul (transpose V_T) S_inv |> flip (LA.matmul) (transpose U)
				in LA.matmul tmp_Q R_inv

			in LA.matmul left Z |> flip (LA.matmul) right

		let (U_hat, S_hat, V_T_hat) = SVD.svd C
		let U = LA.matmul Q U_hat
		let V_T = LA.matmul P (transpose V_T_hat)
		in (U, S_hat, (transpose V_T))
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

-- | Turnstile RSVD with Gaussian embeddings.
module mk_sketchy_svd_default (R: real) : rsvd_turnstile with t = R.t = {
	type t = R.t
	local module EMB = mk_gaussian_embedding R u32 squares32
	local module SKETCHY = mk_sketchy_svd R EMB

	module dist = SKETCHY.dist

	open SKETCHY
}
