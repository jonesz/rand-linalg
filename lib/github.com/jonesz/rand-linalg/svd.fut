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

