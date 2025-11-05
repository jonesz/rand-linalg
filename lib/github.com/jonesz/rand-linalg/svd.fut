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
		-> (*[k][n]t, *[m][k]t, *[2 * k][2 * k]t)

	val linear_update_row [m] [n] [k] : seeds -> (X: *[k][n]t) -> (Y: *[m][k]t) -> (Z: *[2 * k][2 * k]t) -> (row: i64) -> (H: [n]t)
		-> ([k][n]t, [m][k]t, [2 * k][2 * k]t)

	val linear_update_entry [m] [n] [k] : seeds -> (X: *[k][n]t) -> (Y: *[m][k]t) -> (Z: *[2 * k][2 * k]t) -> (idx: (i64, i64)) -> (H: t)
		-> ([k][n]t, [m][k]t, [2 * k][2 * k]t)

	-- TODO: `r` should be the target rank, which should be a little bit below `k`.
	val sketchy_svd [m] [n] [k] : seeds -> (X: [k][n]t) -> (Y: [m][k]t) -> (Z: [2 * k][2 * k]t) -> (_r: i64)
		-> ([m][k]t, [k][k]t, [k][n]t)
}

-- | Compute an economical SVD via one-sided Jacobi iterations.
-- this module is serial; that is, the scheduling for columns is a simple `foldl`.
module mk_one_sided_jacobi_slow (R: real) : {
  type t
  val svd [l] : [l][l]t -> ([l][l]t, [l][l]t, [l][l]t)
} with t = R.t = {
	type t = R.t
	module LA = mk_linalg R

    -- Build the Jacobi rotation matrix.
	let jacobi (l: i64) (cs: t) (sn: t) (p: i64) (q: i64) =
		(LA.eye l) with [p, p] = cs
				   with [p, q] = sn
				   with [q, p] = sn |> R.neg
				   with [q, q] = cs

	let serial_schedule l =
		-- TODO: The `filter` sucks but we should try to move off the serial schedule.
		map (\j -> map (\k -> (j, k)) (iota l)) (iota l) |> flatten |> filter (\(a, b) -> a < b)

    -- | Compute (U, S, V_T).
	def svd [l] (A: [l][l]t) : ([l][l]t, [l][l]t, [l][l]t) = 
		-- [Algorithm 6](https://icl.utk.edu/files/publications/2018/icl-utk-1341-2018.pdf)
		let (A_hat, V) =
          -- TODO: Check for actual convergence rather than doing X amount of sweeps.
		  loop (A_k, V_k) = (A, LA.eye l) for _i < 10 do
			foldl (\(A_k, V_k) (i, j) ->
				let col_i = A_k[0:, i]
				let col_j = A_k[0:, j]

				let b_ij = LA.dotprod col_i col_j
				let b_jj = LA.dotprod col_j col_j 
				let b_ii = LA.dotprod col_i col_i 

				let cond =
					let eps = R.f32 1e-4f32
					in (R.>=) (R.abs b_ij) ((R.*) b_jj b_ii |> R.sqrt |> (R.*) eps)

				in if cond
				then
					-- Compute a 2x2 Schur Decomposition (8.5.2).
					-- https://math.ecnu.edu.cn/~jypan/Teaching/books/2013%20Matrix%20Computations%204th.pdf
					let tau = (R.-) b_jj b_ii |> flip (R./) ((R.*) (R.i64 2) b_ij)
					let t =
						if (R.>=) tau (R.i64 0)
						then
							(R.*) tau tau |> (R.+) (R.i64 1) |> R.sqrt |> (R.+) tau |> (R./) (R.i64 1) 
						else
							(R.*) tau tau |> (R.+) (R.i64 1) |> R.sqrt |> (R.-) tau |> (R./) (R.i64 1)

					let cs = (R.*) t t |> (R.+) (R.i64 1) |> R.sqrt |> (R./) (R.i64 1)
					let sn = (R.*) cs t

                    -- TODO: Update specific columns rather than computing the full
                    -- matmul.

					let J = jacobi l cs sn i j
					let A_k = LA.matmul A_k J
					let V_k = LA.matmul V_k J

					in (A_k, V_k)
				else
					(A_k, V_k)

			) (A_k, V_k) (serial_schedule l)

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

module mk_sketchy_svd(R: real) (SK: sketch with t = R.t): rsvd_turnstile = {
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
module mk_sketchy_svd_default (R: real) : rsvd_turnstile = {
	type t = R.t
	local module EMB = mk_gaussian_embedding R u32 squares32
	local module SKETCHY = mk_sketchy_svd R EMB

	module dist = SKETCHY.dist


	open SKETCHY
}
